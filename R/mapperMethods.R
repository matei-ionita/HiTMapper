#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @import ggraph
#' @import leidenAlg
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
                       nx=10, ny=10, overlap=0.1,
                       nodemax=1000, kmax = 20,
                       intersectionSize=10) {
  pr <- prcomp(data, rank. = 2, scale. = scale)
  bins <- getBins(filter = pr$x, nx=nx, ny=ny, overlap=overlap)
  nodes <- clusterFibers(data, bins, nodemax=nodemax, kmax=kmax)
  gr <- getGraph(nodes, data, M=intersectionSize)
  nodeMedians <- getMedians(data, nodes)

  mapper <- list(bins=bins, nodes=nodes, gr=gr, nodeMedians=nodeMedians)
  return(mapper)
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


#' @title getSampleFeatures
#' @description Create a contingency table of sample vs node membership.
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



