#' @import tidyverse
#' @import igraph
#' @import ggraph
#' @importFrom NetworkToolbox louvain
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
flowMapper <- function(data, nx=10, ny=10, overlap=0.1, nodemax=1000, intersectionSize=10) {
  pr <- prcomp(data, rank. = 2, scale. = TRUE)
  bins <- getBins(filter = pr$x, nx=nx, ny=ny, overlap=overlap)
  nodes <- clusterFibers(data, bins, nodemax=nodemax)
  gr <- getGraph(nodes, data, M=intersectionSize)
  nodeMedians <- getMedians(data, nodes)

  mapper <- list(bins=bins, nodes=nodes, gr=gr, nodeMedians=nodeMedians)
  return(mapper)
}

#' @title louvainClustering
#' @description Meta-clustering of Mapper nodes using Louvain method.
#' @param graph An igraph object, for example, from the output of flowMapper.
#' @param gamma Numeric value controlling the number of meta-clusters:
#' increase for more, decrease for fewer.
#' @export
louvainClustering <- function(graph, gamma=1) {
  adj <- as_adjacency_matrix(graph, type = "both")
  louv <- louvain(as.matrix(adj), gamma=gamma)
  clust <- as.factor(louv$community)
  return(clust)
}
