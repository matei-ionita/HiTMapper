#' @useDynLib HiTMapper, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise filter mutate
#' @import igraph
#' @import ggraph
#' @importFrom mclust Mclust mclustBIC
#' @importFrom glmnet glmnet
#' @importFrom KernSmooth bkde
#' @import leidenAlg
NULL


#' @title HiTMapper
#' @description Wrapper for the core Mapper functionality.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param total_nodes Approximate number of nodes for the Mapper graph.
#' @param method Clustering method to use for each bin. Currently implemented
#' options are "som" (recommended) or "kmeans".
#' @param min_node_size Minimum number of cells in a node. No need to change
#' the default, unless there are very few total cells.
#' @param grid_size A vector of integers, specifying the number of
#' overlapping bins along each principal component. Its length determines
#' the number of principal components used.
#' @param overlap Fraction of overlap between neighboring bins.
#' @param n_passes Number of passes over all data points when allocating
#' cells to nodes via the "som" method.
#' @param scale Logical, whether to scale the data before PCA.
#' @param filter A matrix with few columns and same number of rows as data.
#' To be used as filter, instead of the default principal components.
#' @param resolution Used in community detection. Larger values lead to more
#' communities.
#' @param k Number of nearest neighbors to use when constructing graph.
#' @param force_sep Character vector, must be a subset of colnames(data).
#' Weakens graph edges between nodes that have different modality for any
#' of the given markers. This is a kludge and should be used sparingly.
#' @export
HiTMapper <- function(data, total_nodes,
                      method = "som",
                      min_node_size=50,
                      grid_size=c(10,10), overlap=0,
                      n_passes=10,
                      scale = FALSE, filter=NULL,
                      resolution=8, defs=NULL,
                      k=8, force_sep=c()) {

  if(!is.matrix(data))
    stop("Please enter your data in matrix format.")

  if (is.null(filter)) {
    message("Computing filter function...")
    rank <- length(grid_size)
    filter <- get_filter(data, rank, scale)
  }

  message("Computing and clustering level sets...")
  level_sets <- get_level_sets(filter = filter, grid_size=grid_size, overlap=overlap)
  bins <- apply_level_sets(filter = filter, boundaries = level_sets)

  nodes <- cluster_level_sets(data, bins, total_nodes=total_nodes,
                              min_node_size=min_node_size,
                              n_passes=n_passes,
                              verbose=verbose,method=method)
  node_stats <- get_stats(data, nodes)
  thresholds <- get_modality_thresholds_gmm(data, nodes)

  message("Constructing graph...")
  gr <- get_graph_nn(node_stats$q50, k=k, thresholds, force_sep)
  mapper <- list(nodes=nodes, gr=gr, node_stats=node_stats,
                 thresholds=thresholds)

  message("Detecting communities...")
  mapper <- detect_communities(mapper, data, resolution, defs)

  message(paste("Done! Graph has", length(nodes), "nodes,", length(E(gr)),
                "edges,", length(levels(mapper$community)),"communities."))

  return(mapper)
}


#' @title detect_communities
#' @description Use Leiden algorithm to identify communities in the
#' HiTMapper network. A community is a group of similar nodes,
#' which approximates the concept of "cell type".
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param mapper Existing mapper object.
#' @param resolution passed to the Leiden algorithm, controls
#' the number of communities. Increase for more, decrease for fewer.
#' @param defs A data frame of cell population definitions.
#' @export
detect_communities <- function(mapper, data, resolution=8, defs=NULL) {
  community <- leiden_clustering(mapper$gr, resolution)
  community_medians <- get_community_medians(community, mapper$nodes, data)

  diff <- sweep(community_medians, 2, mapper$thresholds)
  modality <- diff
  modality[which(diff<=0)] <- "lo"
  modality[which(diff>0)] <- "hi"
  
  mapper$community <- community
  mapper$community_medians <- community_medians
  mapper$modality <- modality
  mapper$clustering <- assign_cells(data, mapper, community)
  
  if (!is.null(defs)) {
    mapper <- label_communities(mapper, defs)
    
    mapper <- distinguish_communities(mapper, data)
    mapper$community_medians <- get_community_medians(mapper$community, 
                                                      mapper$nodes, data)
    row.names(mapper$community_medians) <- levels(mapper$community)
  }
  
  return(mapper)
}


#' @title extract_features
#' @description Extract features (cell type proportions across samples)
#' node compositions (node proportions across samples).
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param sample_mapping A vector giving the sample membership of
#' each data point.
#' @param mapper Existing mapper object.
#' @export
extract_features <- function(data, sample_mapping, mapper) {
  mapper$features <- get_contingency_table(mapper$clustering,
                                           sample_mapping)
  mapper$node_composition <- node_composition(mapper, sample_mapping)
  mapper$sample_names <- unique(sample_mapping)
  return(mapper)
}


#' @title label_communities
#' @description Wrapper for community detection, labeling,
#' and extracting features (cell type percentages).
#' @param mapper Existing mapper object.
#' @param defs Data frame of phenotype definitions.
#' @param additional Named list, specifying attributes to be
#' appended to phenotypes. E.g. list(CD38="activated", Ki67="proliferating").
#' @export
label_communities <- function(mapper, defs, additional=list()) {
  matches <- match_defs(defs, mapper$modality)
  phenos <- c(defs$Phenotype, "Other")
  labels <- phenos[matches]
  labels <- append_additional(additional, mapper$modality, labels, matches)
  labels <- make.unique(labels)

  levels(mapper$community) <- labels
  row.names(mapper$community_medians) <- labels
  if (!is.null(mapper$features))
    colnames(mapper$features) <- c(labels, "Unassigned")
  if (!is.null(mapper$clustering))
    levels(mapper$clustering) <- c(levels(mapper$community), "Unassigned")

  return(mapper)
}

