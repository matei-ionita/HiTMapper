#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @import ggraph
#' @import leidenAlg
#' @import e1071
#' @import rdist
NULL


#' @title flowMapper
#' @description Wrapper for the core Mapper functionality.
#' @param data A data matrix.
#' @param nx Number of bins (level sets) for the first
#' filter dimension.
#' @param ny Number of bins (level sets) for the second
#' filter dimension.
#' @param overlap Fractional overlap between neighboring bins.
#' @param nodemax In level set clustering stage, split level
#' sets with larger cardinality.
#' @param intersectionSize Minimum number of shared datapoints
#' required for an edge.
#' @export
flowMapper <- function(data, scale = TRUE,
                       nx=10, ny=10, overlap=0.3,
                       nodemax=1000, kmax = 20,
                       intersectionSize=NULL,
                       iou=0.03,
                       clusterMethod="kmeans",
                       outlierCutoff=25,
                       verbose=FALSE) {

  message("Computing filter function and level sets...")
  pr <- prcomp(data, rank. = 2, scale. = scale)
  levelSets <- getLevelSets(filter = pr$x, nx=nx, ny=ny, overlap=overlap)
  bins <- applyLevelSets(filter = pr$x, levelSets = levelSets)

  message("Clustering the level sets...")
  cluster <- clusterFibers(data, bins, method=clusterMethod, nodemax=nodemax,
                           kmax=kmax, verbose=verbose, outlierCutoff=outlierCutoff)

  message("Comstructing Mapper graph...")
  gr <- getGraph(cluster, M=intersectionSize, iou=iou)

  cluster$nodes <- mapToCenters(data, bins, cluster$centers)
  nodeMedians <- getMedians(data, cluster$nodes)

  message(paste("Done! Graph has", length(cluster$nodes), "nodes,",
                length(E(gr)), "edges."))
  pr$x <- NULL
  mapper <- list(bins=bins, nodes=cluster$nodes, gr=gr, nodeMedians=nodeMedians,
                 centers=cluster$centers, pr=pr, levelSets=levelSets)
  return(mapper)
}



resetEdges <- function(mapper, intersectionSize=NULL, iou=NULL) {
  cluster <- list()
  cluster$nodes <- mapper$nodes
  cluster$edges <- matrix(nrow=0,ncol=2)

  gr <- getGraph(cluster, M=intersectionSize, iou=iou)
  mapper$gr <- gr
  return(mapper)
}


#' @title applyMapper
#' @description Map new data to existing Mapper network.
#' @param data A data matrix.
#' @param mapper Existing mapper network
#' @export
applyMapper <- function(data, mapper) {
  filter <- applyPCA(data,mapper$pr)
  bins <- applyLevelSets(filter = filter, levelSets = mapper$levelSets)
  nodes <- mapToCenters(data, bins, mapper$centers)
  sizes <- sapply(nodes, length)

  return(list(bins = bins, nodes = nodes, sizes = sizes))
}

applyPCA <- function(data, pr) {

  scaled <- data %>%
    as.matrix() %>%
    sweep(MARGIN = 2, STATS = pr$center) %>%
    sweep(MARGIN = 2, STATS = pr$scale, FUN = "/")

  filter <- scaled %*% pr$rotation
  return(filter)
}

getMapping <- function(data, nodes, community) {
  mapping <- character(nrow(data))

  for (cl in levels(community)) {
    mapping[do.call(nodes[which(community == cl)], what=c)] <- cl
  }

  return(mapping)
}


#' @title louvainClustering
#' @description Meta-clustering of Mapper nodes using Leiden method.
#' @param graph An igraph object, for example, from the output of flowMapper.
#' @param resolution Numeric value controlling the number of meta-clusters:
#' increase for more, decrease for fewer.
#' @export
leidenClustering <- function(graph, resolution=1) {
  leid <- leiden.community(graph, resolution=resolution)
  clust <- as.factor(leid$membership)
  return(clust)
}


#' @title Basic features.
#' @description Extract features from multiple samples, as a contingency table
#' of sample vs node membership.
#' @param nodes A list of integer vectors, each containing the data
#' points assigned to each node.
#' @param sampleMapping An integer vector, containing the sample of
#' origin for each data point.
#' @param normalized Logical, whether the contingency table should
#' be normalized across features. (Normalization across samples takes
#' place either way.
#' @export
getSampleFeatures <- function(nodes, sampleMapping, cl, normalized=TRUE) {
  pctg <- sapply(nodes, function(node) sampleMapping[node] %>%
                   tabulate(nbins = 112)) %>% t()

  unq <- as.integer(as.character(unique(cl)))
  tm <- sapply(unq, function(i) apply(pctg[which(cl==i),,drop=FALSE], 2, sum))

  pctg <- apply(tm, 1, function(row) row/sum(row)) %>% t()

  return(data.frame(pctg))

  # pctg <- apply(pctg, 2, function(col) col/sum(col))
  # pctgNorm <- apply(pctg, 1, function(row) row/mean(row))
  #
  # if(normalized)
  #   return(pctgNorm)
  # else
  #   return(pctg)
}

#' @title Moving average features.
#' @description Extract features from multiple samples, as a contingency table
#' of sample vs node membership. Smoothed over neighboring nodes in the network
#' for de-noising purposes, analogous to a moving average.
#' @param nodes A list of integer vectors, each containing the data
#' points assigned to each node.
#' @param sampleMapping An integer vector, containing the sample of
#' origin for each data point.
#' @param normalized Logical, whether the contingency table should
#' be normalized across features. (Normalization across samples takes
#' place either way.
#' @export
getSampleFeaturesSmooth <- function(mapper, sampleMapping, norm=FALSE) {
  pctg <- sapply(mapper$nodes, function(node) sampleMapping[node] %>%
                   tabulate(nbins = 112))

  pctg <- apply(pctg, 1, function(row) row / sum(row)) %>% t()

  adj <- as_adj_list(mapper$gr)
  pctg <- sapply(seq_along(adj), function(i) {
    apply(pctg[,c(i,adj[[i]]),drop=FALSE], 1, sum)
  })

  if(norm) {
    pctg <- apply(pctg, 1, function(row) row / sum(row)) %>% t()
    pctg <- apply(pctg, 2, function(col) col / sum(col))
  }

  return(data.frame(pctg))
}



