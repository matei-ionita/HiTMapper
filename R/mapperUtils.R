getLevelSets <- function(filter, nx = 10, ny = 10, overlap = 0.1) {
  # Split the dimensionally reduced space into overlapping bins
  # Assuming 2D. Setting ny=1 is effectively 1D.

  minx <- min(filter[,1])
  maxx <- max(filter[,1])

  miny <- min(filter[,2])
  maxy <- max(filter[,2])

  sizex <- (maxx-minx) / (nx*(overlap+1)+overlap)
  sizey <- (maxy-miny) / (ny*(overlap+1)+overlap)

  x0 <- minx + seq(0,9)*sizex*(1+overlap)
  x1 <- x0 + sizex*(1+2*overlap)

  y0 <- miny + seq(0,9)*sizey*(1+overlap)
  y1 <- y0 + sizey*(1+2*overlap)

  levelSets <- list(x0=x0, x1=x1, y0=y0, y1=y1)
  return(levelSets)
}


applyLevelSets <- function(filter, levelSets, outlierCutoff = 0) {
  nx <- length(levelSets$x0)
  ny <- length(levelSets$y0)

  binsx <- lapply(seq(nx), function(i) which(filter[,1] >= levelSets$x0[i] &
                                             filter[,1] <= levelSets$x1[i]))
  binsy <- lapply(seq(ny), function(j) which(filter[,2] >= levelSets$y0[j] &
                                             filter[,2] <= levelSets$y1[j]))
  bins <- list()
  for (i in seq(nx))
    for (j in seq(ny)) {
      binid <- j + (i-1)*ny
      bins[[binid]] <- intersect(binsx[[i]], binsy[[j]])
    }

  counts <- sapply(bins, length)
  bins <- bins[which(counts >= outlierCutoff)]

  return(bins)
}


clusterFibers <- function(data, bins, method = "kmeans", outlierCutoff = 50,
                          kNodes=NULL,
                          nodemax, kmax, verbose=verbose) {
  if (method == "kmeans")
    return(clusterFibersKmeans(data, bins, outlierCutoff, nodemax, kmax, verbose=verbose,
                               kNodes))

  if (method == "fuzzy")
    return(clusterFibersFuzzy(data, bins, outlierCutoff, nodemax, kmax, verbose=verbose,
                              kNodes))

  # Enter more methods here

  stop("Invalid method chosen.")
}


clusterFibersFuzzy <- function(data, bins, outlierCutoff, nodemax, kmax, verbose,
                               kNodes=NULL) {
  nodes <- list()
  edges <- matrix(nrow = 0, ncol = 2)

  for (bin in bins) {
    if (length(bin) < outlierCutoff)
      next

    k <- ceiling(length(bin)/nodemax)
    k <- min(k,kmax)

    cat("Splitting", length(bin), "into", k, "nodes.\n")

    if (k == 1) {
      nodes <- c(nodes, list(bin))
      next
    }

    fuzzy <- cmeans(data[bin,], centers=k, m=1.5, iter.max=50)
    newNodes <- fuzzy$cluster %>%
      nodesFromMapping(bin = bin)

    keep <- which(sapply(newNodes, length) >= outlierCutoff)

    if (length(keep) > 1) {
      newEdges <- cosineDistance(fuzzy$membership[,keep])
      newEdges <- newEdges + length(nodes)
      edges <- rbind(edges, newEdges)
    }

    nodes <- c(nodes, newNodes[keep])
  }

  return(list(nodes=nodes, edges=edges))
}


cosineDistance <- function(X) {
  inn <- t(X) %*% X
  norms <- diag(inn)
  ans <- diag(1/sqrt(norms)) %*% inn %*% diag(1/sqrt(norms))
  edges <- which(ans > 0.05, arr.ind = TRUE)
  edges <- edges[which( edges[,"row"]< edges[,"col"] ),,drop=FALSE]

  return(edges)
}


getK <- function(data, bins, outlierCutoff, kNodes,
                 fracSumSq = 0.7) {
  binLength <- sapply(bins, length)
  notOut <- which(binLength>= outlierCutoff)

  binSumSq <- sapply(bins, function(bin) {
    if (length(bin) < 2)
      return(0)

    v <- var(data[bin,])
    return(sum(diag(v)))
  }) %>% sqrt()
  kSumSq <- fracSumSq * binSumSq * kNodes / sum(binSumSq[notOut])

  # binDist <- sapply(bins, function(bin) {
  #   if (length(bin) < 2)
  #     return(0)
  #
  #   m <- apply(data[bin,],2,mean)
  #   d <- abs(t(data[bin,])-m) %>%
  #     apply(2,sum) %>%
  #     mean()
  #
  #   return(d)
  # })
  # kDist <- fracSumSq * binDist * kNodes / sum(binDist[notOut])

  binLog <- pmax(log2(binLength),0)
  kLog <- (1-fracSumSq) * binLog * kNodes / sum(binLog[notOut])

  binK <- ceiling(kSumSq + kLog)
  allK <- cbind(binLength, kLog, kSumSq, binK)
  # print(allK[notOut,])
  return(data.frame(allK))
}


clusterFibersKmeans <- function(data, bins, outlierCutoff, nodemax, kmax, verbose,
                                kNodes=NULL) {
  nodes <- list()
  centers <- list()

  allK <- getK(data, bins, outlierCutoff, kNodes)

  for (i in seq_along(bins)) {
    bin <- bins[[i]]
    if (length(bin) < outlierCutoff)
      next

    # k <- ceiling(length(bin)/nodemax)
    # k <- min(k,kmax)
    k <- allK$binK[i]

    if (k == 1) {
      nodes <- c(nodes, list(bin))
      centers[[i]] <- data[1,,drop=FALSE]-data[1,,drop=FALSE]
      next
    }

    if(verbose)
      message(paste("Bin size", length(bin), "k", k))

    km <- suppressWarnings(kmeans(data[bin,], centers=k, nstart = 10))

    newNodes <- km$cluster %>%
      nodesFromMapping(bin = bin, nCent = k)

    keep <- which(sapply(newNodes, length) >= outlierCutoff)
    centers[[i]] <- km$centers[keep,,drop=FALSE]
    nodes <- c(nodes, newNodes[keep])
  }

  return(list(centers = centers, nodes = nodes,
              edges = matrix(nrow=0,ncol=2)))
}

mapToCenters <- function(data, bins, centers) {
  nodes <- list()

  for (i in seq_along(bins)) {
    bin <- bins[[i]]
    if (i > length(centers) || is.null(centers[[i]]))
      next



    d <- cdist(data[bin,],centers[[i]])
    cluster <- apply(d, 1, which.min)
    newNodes <- nodesFromMapping(mapping=cluster, bin=bin, nCent=nrow(centers[[i]]))
    nodes <- c(nodes, newNodes)
  }

  return(nodes)
}


nodesFromMapping <- function(mapping, bin, nCent) {
  nodes <- list()

  for (i in seq_len(nCent)) {
    mapi <- which(mapping == i)
    nodes[[i]] <- bin[mapi]
  }
  # unq <- unique(mapping)
  # nodes <- lapply(unq, function(i) bin[which(mapping==i)])

  return(nodes)
}



getInters <- function(cluster, mode="iou") {
  nodes <- cluster$nodes

  n <- length(nodes)
  sizes <- sapply(nodes, length)
  inters <- matrix(0, nrow=n, ncol=n)

  for (i in seq(1,n-1)) {
    node <- nodes[[i]]
    for (j in seq(i+1,n)) {
      neighbor <- nodes[[j]]
      int <- intersect(node, neighbor)

      if (mode == "iou") {
        uni <- union(node, neighbor)
        inters[i,j] <- length(int) / length(uni)
        inters[j,i] <- length(int) / length(uni)
      } else {
        inters[i,j] <- length(int)
        inters[j,i] <- length(int)
      }
    }
  }

  return(inters)
}

getAdjList <- function(inters, M, iou, mode="iou") {
  n <- nrow(inters)
  adjList <- vector(mode = "list", length = n)

  if (mode == "iou") {
    for (i in seq(n)) {
      adjList[[i]] <- which(inters[i,] > iou)
    }
  } else {
    for (i in seq(n)) {
      adjList[[i]] <- which(inters[i,] > M)
    }
  }

  return(adjList)
}


getGraph <- function(inters, M = 10, iou=NULL) {
  # adjList <- getAdjList(inters, M, iou)
  # gr <- graph_from_adj_list(adjList, mode = "all")
  gr <- graph_from_adjacency_matrix(inters/max(inters),
                                    weighted = TRUE, mode="undirected")

  # E(gr)$weight <- rep(1, length(E(gr))) # for Leiden
  # print(nrow(inters) / length(E(gr)))

  return(gr)
}


getStats <- function(data, nodes) {
  q25 <- sapply(nodes, function(node) apply(data[node,], 2, function(x) quantile(x,0.25))) %>% t()
  q50 <- sapply(nodes, function(node) apply(data[node,], 2, median)) %>% t()
  q75 <- sapply(nodes, function(node) apply(data[node,], 2, function(x) quantile(x,0.75))) %>% t()
  return(list(q25=q25,q50=q50,q75=q75))
}


