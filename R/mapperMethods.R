#' @useDynLib HiTMapper, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise filter select mutate if_else pull
#' @importFrom stringr str_split fixed
#' @importFrom tibble tibble as_tibble
#' @import igraph
#' @import ggraph
#' @importFrom leidenAlg leiden.community
NULL


#' @title HiTMapper
#' @description Wrapper for the core Mapper functionality.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param total_nodes Approximate number of nodes for the Mapper graph.
#' @param resolution Used in community detection. Larger values lead to more
#' communities. We recommend running HiTMapper with default resolution,
#' and later customizing if necessary using the detect_communities function,
#' at very little computational cost.
#' @param defs A data frame of cell population definitions, to be used for
#' labeling communities.
#' @param grid_size Array of integers, the number of bins along 
#' each principal component.
#' @param min_node_size Minimum number of cells in a node.
#' @param n_passes Number of passes over all data points when clustering.
#' @export
HiTMapper <- function(data, total_nodes=1000,
                      resolution=1, defs=NULL,
                      grid_size=c(10,10),
                      min_node_size=50, n_passes=10) {

  if(!is.matrix(data))
    stop("Please enter your data in matrix format.")

  message("Clustering and creating graph...")
  centroids <- clustering_main(data, cov(data), grid_size, 
                               total_nodes, min_node_size, n_passes)
  message("Mapping...")
  l <- assign_datapoints(data, centroids)
  mapping <- l$mapping
  sim <- l$sim
  
  tab <- tabulate(mapping)
  zer <- which(tab==0)

  # remove nodes which contain no data points
  if(length(zer) > 0) {
    centroids <- centroids[-zer,]
    sim <- sim[-zer,-zer]

    for (z in sort(zer, decreasing = TRUE)) {
      mapping[which(mapping >= z)] <- mapping[which(mapping >= z)]-1
    }
    tab <- tabulate(mapping)
  }

  sim <- get_weights(tab, sim, min(min_node_size, ceiling(nrow(data)/1000)))
  colnames(centroids) <- colnames(data)
  gr <- get_graph_sim(sim)
  mapper <- list(gr=gr, centroids=centroids,
                 mapping=mapping)
  if(!is.null(defs))
    mapper$defs <- defs
  
  message("Detecting communities...")
  mapper <- detect_communities(mapper, data, resolution)

  return(mapper)
}


#' @title Use existing mapper model to predict cluster membership of new data.
#' @param mapper Existing mapper object.
#' @param data New data matrix, must have same columns as that used for
#' mapper model.
#' @export
HiTMapper_predict <- function(mapper, data) {
  mapper$mapping <- predict_datapoints(data, mapper$centroids)
  mapper$clustering <- mapper$community[mapper$mapping]
  return(mapper)
}



#' @title Perform community detection on the existing graph
#'  with a new resolution, then re-label the communities.
#' @param mapper Existing mapper object.
#' @param data The data matrix used for obtaining mapper object.
#' @param resolution Numeric, resolution for the Leiden algorithm.
#' @export
detect_communities <- function(mapper, data, resolution) {
  mapper$community <- leiden_clustering(mapper$gr, resolution)
  mapper$clustering <- mapper$community[mapper$mapping]
  mapper <- parse_communities(mapper, data)

  message(paste("Done! Graph has", length(V(mapper$gr)), "nodes,",
                length(E(mapper$gr)), "edges,",
                length(levels(mapper$community)),"communities."))
  return(mapper)
}


#' @title label_communities
#' @description Wrapper for community detection, labeling,
#' and extracting features (cell type percentages).
#' @param mapper Existing mapper object.
#' @param defs Data frame of phenotype definitions.
#' @param thresholds Data frame of user-supplied thresholds,
#' for a subset of markers,
#' to override those estimated by the algorithm.
#' @export
label_communities <- function(mapper, defs, thresholds=NULL) {
  mapper$thresholds <- get_modality_thresholds_hier(mapper$centroids)
  
  if(!is.null(thresholds))
    mapper$thresholds[thresholds$channel] <- thresholds$value
  
  diff <- sweep(mapper$community_medians, 2, mapper$thresholds)
  mapper$modality <- ifelse(diff<=0, "lo", "hi")
  
  matches <- match_defs(defs, mapper$modality)
  phenos <- c(defs$Phenotype, "Other")
  labels <- phenos[matches]
  labels <- make.unique(labels)

  levels(mapper$community) <- labels
  row.names(mapper$community_medians) <- labels
  levels(mapper$clustering) <- labels

  if (!is.null(mapper$features))
    colnames(mapper$features) <- labels

  return(mapper)
}


#' @title plot_mapper
#' @description Generates plots of the network color-coded by
#' each of the specified markers.
#' @param mapper Existing mapper object.
#' @param markers A subset of the marker names.
#' @param path The path where plots are saved.
#' @param device The device used for image encoding. Supports
#' png for bitmaps, pdf for vector graphics.
#' @export
plot_mapper <- function(mapper, markers=colnames(mapper$centroids),
                        path = NULL, device="png") {
  
  if (device!="png" & device!="pdf")
    stop("Supported devices are png and pdf.")
  
  size <- as.integer(table(mapper$mapping))
  layout <- create_layout(mapper$gr, layout="fr")
  
  comm <- data.frame(x=layout$x,y=layout$y,
                     c=mapper$community) %>%
    group_by(c) %>%
    summarise(x=median(x), y=median(y))
  
  for (marker in markers) {
    g <- ggraph(layout) +
      geom_edge_link(aes(alpha = weight)) +
      geom_node_point(shape=21, aes(size=size, fill = mapper$centroids[,marker])) +
      scale_fill_gradient2(low = "white", mid="white",
                           high = "red",name = marker) +
      scale_size(range=c(1,6), name="count") +
      scale_edge_alpha(guide="none") +
      theme_graph(base_family = "sans") +
      geom_label(data=comm, aes(x=x,y=y,label=c),
                 alpha=0.5, inherit.aes = FALSE)
    
    if (is.null(path)) {
      plot(g)
    } else {
      ggsave(plot=g, filename=paste0(path,marker, ".", device))
    }
  }
  
  if (!is.null(mapper$community)) {
    g <- ggraph(layout) +
      geom_edge_link(aes(alpha = weight)) +
      geom_node_point(aes(color=mapper$community, size=size)) +
      scale_color_discrete(name="community") +
      scale_edge_alpha(guide="none") +
      scale_size(range=c(1,6), name="count") +
      theme_graph(base_family = "sans") +
      theme(text=element_text(size = 18))
    
    if (is.null(path)) {
      plot(g)
    } else {
      ggsave(plot=g, filename=paste0(path,"community.", device))
    }
  }
}
