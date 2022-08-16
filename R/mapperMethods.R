#' @useDynLib HiTMapper, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise
#' @import igraph
#' @import ggraph
#' @import leidenAlg
NULL


#' @title HiTMapper
#' @description Wrapper for the core Mapper functionality.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param total_nodes Approximate number of nodes for the Mapper graph.
#' @param outlier_cutoff Discard nodes containing fewer cells than this.
#' @param grid_size A vector of integers, specifying the number of
#' overlapping bins along each principal component. Its length determines
#' the number of principal components used.
#' @param overlap Fraction of overlap between neighboring bins.
#' @param scale Logical, whether to scale the data before PCA.
#' @param verbose Logical, whether to output detailed info.
#' @export
HiTMapper <- function(data, total_nodes,
                      method = "som", graph = "knn",
                      outlier_cutoff=nrow(data)/1e4,
                      grid_size=c(10,10), overlap=0.15,
                      n_passes=10,
                      scale = FALSE, verbose=FALSE, filter=NULL,
                      merge=TRUE, split = c(), resolution=2,
                      k=8) {

  if(!is.matrix(data))
    stop("Please enter your data in matrix format.")

  if (is.null(filter)) {
    message("Computing filter function...")
    rank <- length(grid_size)
    filter <- get_filter(data, rank, scale)
  }

  message("Computing and clustering level sets...")
  if (graph == "knn")
    overlap <- 0

  level_sets <- get_level_sets(filter = filter, grid_size=grid_size, overlap=overlap)
  bins <- apply_level_sets(filter = filter, boundaries = level_sets)

  nodes <- cluster_level_sets(data, bins, total_nodes=total_nodes,
                              outlier_cutoff=outlier_cutoff,
                              n_passes=n_passes,
                              verbose=verbose,method=method)
  node_stats <- get_stats(data, nodes)
  thresholds <- get_modality_thresholds(node_stats$q50)

  message("Constructing Mapper graph...")

  if (graph == "knn") {
    gr <- get_graph_nn(node_stats$q50, k=k)
  } else {
    gr <- get_graph_c(nodes, node_stats)

  }

  mapper <- list(nodes=nodes, gr=gr, node_stats=node_stats,
                 thresholds=thresholds)

  message("Detecting communities...")
  mapper <- detect_communities(mapper, data, resolution, merge, split)

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
#' @param merge Logical, whether to merge communities which have the same
#' modality across all markers in the data. It is recommended to use a
#' large resolution to obtain many communities, some of which are then merged.
#' @export
detect_communities <- function(mapper, data, resolution=2, merge=TRUE,
                               split = c(), min_split=0) {
  community <- leiden_clustering(mapper$gr, resolution)

  for (split_marker in split) {
    message(paste("Splitting", split_marker))
    community <- split_communities(community, mapper, split_marker, min_split)
  }

  community_medians <- get_community_medians(community, mapper$nodes, data)

  diff <- sweep(community_medians, 2, mapper$thresholds)
  modality <- diff
  modality[which(diff<=0)] <- "lo"
  modality[which(diff>0)] <- "hi"

  if (merge) {
    group_index <- merge_communities(modality)
    community <- group_index[community] %>% as.factor()
    community_medians <- get_community_medians(community, mapper$nodes, data)

    dup <- duplicated(group_index)
    modality <- modality[which(!dup),]
  }

  mapper$community <- community
  mapper$community_medians <- community_medians
  mapper$modality <- modality
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
  mapper$clustering <- assign_cells(data, mapper, mapper$community)
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
  # labels <- make_unique_modality(labels, mapper$modality)

  levels(mapper$community) <- labels
  row.names(mapper$community_medians) <- labels
  if (!is.null(mapper$features))
    colnames(mapper$features) <- c(labels, "Unassigned")
  if (!is.null(mapper$clustering))
    levels(mapper$clustering) <- c(levels(mapper$community), "Unassigned")

  return(mapper)
}


#' @title getFilter
#' @description Get the filter function (PCA projection) used by HiTMapper.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param scale Logical, whether to scale the data before PCA.
#' @param plot_path Generate plots of the filter color-coded by
#' each of the variables in the data matrix and save them at this path.
#' No plots generated unless path is provided.
#' @export
get_filter <- function(data, rank=2, scale=FALSE, plot_path=NULL) {
  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")

  if (scale) {
    cov <- cor(data)
  } else {
    cov <- cov(data)
  }
  pr <- princomp(covmat=cov, scores=FALSE)
  pr$center <- apply(data, 2, mean)
  filter <- apply_PCA(data, pr, rank)

  # To do: show the overlapping bins on the PCA plot,
  # to help visualize what Mapper is doing
  if(!is.null(plot_path))
    plotFilter(data,filter,plot_path)

  return(filter)
}


