#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @import ggraph
#' @import leidenAlg
#' @import e1071
#' @import rdist
NULL


#' @title HiTMapper
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
HiTMapper <- function(data, scale = FALSE,
                       nx=10, ny=10, overlap=0.3,
                       nodemax=NULL, kmax = NULL,
                       kNodes=500,
                       intersectionSize=NULL,
                       iou=0.03,
                       clusterMethod="kmeans",
                       outlierCutoff=25,
                       verbose=FALSE) {

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
  cluster <- clusterFibers(data, bins, method=clusterMethod, nodemax=nodemax,
                           kmax=kmax, kNodes=kNodes,
                           verbose=verbose, outlierCutoff=outlierCutoff)

  message("Constructing Mapper graph...")
  inters <- getInters(cluster, mode = "iou")
  gr <- getGraph(inters, M=intersectionSize, iou=iou)

  cluster$nodesNew <- mapToCenters(data, bins, cluster$centers)
  nodeStats <- getStats(data, cluster$nodes)

  message(paste("Done! Graph has", length(cluster$nodes), "nodes,",
                length(E(gr)), "edges."))

  mapper <- list(bins=bins, nodes=cluster$nodes, gr=gr, nodeStats=nodeStats,
                 nodesNew=cluster$nodesNew, inters=inters,
                 centers=cluster$centers, pr=pr, levelSets=levelSets)
  return(mapper)
}


getFilter <- function(data, scale=FALSE) {
  if (scale) {
    cov <- cor(data)
  } else {
    cov <- cov(data)
  }
  pr <- princomp(covmat=cov, scores=FALSE)
  pr$center <- apply(data, 2, mean)
  filter <- flowMapper:::applyPCA(data, pr)

  return(filter)
}

plotFilter <- function(data, filter, base="") {
  sel <- sample(nrow(data), 1e5)
  df <- cbind(filter[sel,], data[sel,])

  for(param in colnames(data)) {
    message(param)
    g <- ggplot(df, aes(x=Comp.1, y=Comp.2)) +
      geom_point(aes_(color=as.name(param))) +
      scale_color_gradient(low = "black", high = "orange", name = param) +
      theme_bw()

    path <- paste0(base, "/", param, ".png")
    png(path, width=1200, height=900)
    plot(g)
    dev.off()
  }
}


#' @title pruneEdges
#' @description Prune edges with low weight, to obtain a sparser network
#' and minimize spurious connections.
#' @param mapper Existing mapper network.
#' @param cutoff Prune edges whose weight is smaller than this fraction of
#' the maximum weight.
#' @export
pruneEdges <- function(mapper, cutoff=0.1) {
  if (cutoff < 0 | cutoff > 1)
    stop("Cutoff must be given as fraction of the max weight; please enter
         a value between 0 and 1.")

  inters <- mapper$inters / max(mapper$inters)
  inters[which(inters<cutoff)] <- 0
  mapper$gr <- getGraph(inters)
  return(mapper)
}


#' @title applyMapper
#' @description Map new data to existing Mapper network.
#' @param data A data matrix.
#' @param mapper Existing mapper network
#' @export
applyMapper <- function(data, mapper) {
  filter <- applyPCA(data,mapper$pr)
  # filter <- predict(mapper$pr, data)
  # filter <- filter$scores[,c(1,2)]
  bins <- applyLevelSets(filter = filter, levelSets = mapper$levelSets)
  nodes <- mapToCenters(data, bins, mapper$centers)
  sizes <- sapply(nodes, length)

  return(list(bins = bins, nodes = nodes, sizes = sizes))
}

applyPCA <- function(data, pr, rank=2) {

  # scaled <- data %>%
  #   as.matrix() %>%
  #   sweep(MARGIN = 2, STATS = pr$center) %>%
  #   sweep(MARGIN = 2, STATS = pr$scale, FUN = "/")

  filter <- sweep(as.matrix(data), MARGIN=2, STATS=pr$center) %*% pr$loadings
  return(filter[,seq(rank)])
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



