#' @import ggplot2
#' @import magrittr
#' @import igraph
#' @import ggraph
#' @import leidenAlg
#' @import e1071
#' @import rdist
NULL


#' @title HiTMapper
#' @description Wrapper for the core Mapper functionality.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param kNodes Approximate number of nodes for the Mapper graph.
#' @param outlierCutoff Discard nodes containing fewer cells than this.
#' @param nx Number of bins for the first filter dimension.
#' @param ny Number of bins for the second filter dimension.
#' @param overlap Fraction of overlap between neighboring bins.
#' @param scale Logical, whether to scale the data before PCA.
#' @param verbose Logical, whether to output detailed info.
#' @export
HiTMapper <- function(data, kNodes, outlierCutoff=25,
                      nx=10, ny=10, overlap=0.3,
                      scale = FALSE, verbose=FALSE) {

  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")

  message("Computing filter function...")
  if (scale) {
    cov <- cor(data)
  } else {
    cov <- cov(data)
  }
  pr <- princomp(covmat=cov, scores=FALSE)
  pr$center <- apply(data, 2, mean)
  filter <- applyPCA(data, pr)

  message("Computing level sets...")
  levelSets <- getLevelSets(filter = filter, nx=nx, ny=ny, overlap=overlap)
  bins <- applyLevelSets(filter = filter, levelSets = levelSets)
  rm(filter)
  gc()

  message("Clustering the level sets...")
  cluster <- clusterLevelSets(data, bins, kNodes=kNodes,
                              outlierCutoff=outlierCutoff,
                              verbose=verbose)

  message("Constructing Mapper graph...")
  inters <- getInters(cluster, mode = "iou")
  gr <- getGraph(inters)

  cluster$nodesNew <- mapToCenters(data, bins, cluster$centers)
  nodeStats <- getStats(data, cluster$nodes)

  message(paste("Done! Graph has", length(cluster$nodes), "nodes,",
                length(E(gr)), "edges."))

  mapper <- list(bins=bins, nodes=cluster$nodes, gr=gr, nodeStats=nodeStats,
                 nodesNew=cluster$nodesNew, inters=inters,
                 centers=cluster$centers, pr=pr, levelSets=levelSets)
  return(mapper)
}



#' @title pruneEdges
#' @description Prune edges with low weight, to obtain a sparser network
#' and minimize spurious connections.
#' @param mapper Existing mapper object.
#' @param cutoff Prune edges whose weight is smaller than this fraction of
#' the maximum weight.
#' @export
pruneEdges <- function(mapper, cutoff=0.01) {
  if (cutoff < 0 | cutoff > 1)
    stop("Cutoff must be given as fraction of the max weight; please enter
         a value between 0 and 1.")

  inters <- mapper$inters / max(mapper$inters)
  inters[which(inters<cutoff)] <- 0

  n0 <- length(E(mapper$gr))
  mapper$gr <- getGraph(inters)
  n1 <- length(E(mapper$gr))

  message(paste("Pruned", n0, "edges down to", n1, "."))
  return(mapper)
}



#' @title leidenClustering
#' @description Community detection using Leiden method.
#' @param graph An igraph object, for example, from the output of HiTMapper.
#' @param resolution Numeric value controlling the number of communities:
#' increase for more, decrease for fewer.
#' @export
leidenClustering <- function(graph, resolution=1) {
  leid <- leiden.community(graph, resolution=resolution)
  clust <- as.factor(leid$membership)
  return(clust)
}


#' @title getFilter
#' @description Get the filter function (PCA projection) used by HiTMapper.
#' @param data A data matrix or data frame, with observations
#' along rows and variables along columns.
#' @param scale Logical, whether to scale the data before PCA.
#' @param plotPath Generate plots of the filter color-coded by
#' each of the variables in the data matrix and save them at this path.
#' No plots generated unless path is provided.
#' @export
getFilter <- function(data, scale=FALSE, plotPath=NULL) {
  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")

  if (scale) {
    cov <- cor(data)
  } else {
    cov <- cov(data)
  }
  pr <- princomp(covmat=cov, scores=FALSE)
  pr$center <- apply(data, 2, mean)
  filter <- applyPCA(data, pr)

  # To do: show the overlapping bins on the PCA plot,
  # to help visualize what Mapper is doing
  if(!is.null(plotPath))
    plotFilter(data,filter,plotPath)

  return(filter)
}



#' @title applyMapper
#' @description Map new data to existing Mapper network.
#' @param data A data matrix or data frame.
#' @param mapper Existing mapper object.
#' @export
applyMapper <- function(data, mapper) {
  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")

  filter <- applyPCA(data,mapper$pr)
  bins <- applyLevelSets(filter = filter, levelSets = mapper$levelSets)
  nodes <- mapToCenters(data, bins, mapper$centers)
  sizes <- sapply(nodes, length)

  return(list(bins = bins, nodes = nodes, sizes = sizes))
}



#' @title assignCells
#' @description Map individual data points to network nodes or communities.
#' @param data A data matrix or data frame.
#' @param mapper Existing mapper object.
#' @param community A factor giving community membership for the nodes.
#' If present, outputs mapping of data points to communities. If absent,
#' outputs mapping of data points to network nodes.
#' @export
assignCells <- function(data, mapper, community=NULL) {
  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")

  mapping <- integer(nrow(data))
  distances <- numeric(nrow(data)) + Inf

  for (i in seq_along(mapper$nodes)) {
    node <- mapper$nodes[[i]]
    m <- mapper$nodeStats$q50[i,]

    distNew <- (t(data[node,])-m)^2 %>% apply(2, sum)
    closer <- which(distNew < distances[node])

    distances[node[closer]] <- distNew[closer]
    mapping[node[closer]] <- i
  }
  distances <- sqrt(distances)

  if (!is.null(community)) {
    # Some cells are unassigned
    n <- length(mapper$nodes)
    mapping[which(mapping==0)] <- n+1
    communityAugmented <- c(as.character(community), "Unassigned") %>%
      as.factor()

    mapping <- communityAugmented[mapping]
  }

  df <- data.frame(mapping=mapping,dist=distances)
  return(df)
}



#' @title nodeComposition
#' @description Contingency table of sample vs node membership.
#' @param mapper Existing mapper object.
#' @param samples An integer vector, containing the sample of
#' origin for each data point.
#' @param scale Logical, whether the contingency table should
#' be scaled across features. (Scaling across samples takes
#' place either way.
#' @export
nodeComposition <- function(mapper, samples, scale=FALSE) {
  if (is.factor(samples)) {
    unq <- levels(samples)
  } else {
    unq <- unique(samples)
  }

  counts <- lapply(mapper$nodes, function(node) {samples[node] %>% table}) %>%
    sapply(function(v) {
      w <- v[unq]
      w[which(is.na(names(w)))] <- 0
      return(unname(w))
    })

  pctg <- apply(counts, 1, function(row) {
    if (sum(row) == 0)
      return(row)
    return(row / sum(row))
  }) %>% t()

  if(scale)
    pctg <- apply(pctg, 2, function(col) {
      if (sum(col) == 0)
        return(col)
      return(col / sum(col))
    })

  return(data.frame(pctg))
}






