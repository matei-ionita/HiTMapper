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
                          nodemax, kmax) {
  if (method == "kmeans")
    return(clusterFibersKmeans(data, bins, outlierCutoff, nodemax, kmax))

  if (method == "fuzzy")
    return(clusterFibersFuzzy(data, bins, outlierCutoff, nodemax, kmax))

  # Enter more methods here

  stop("Invalid method chosen.")
}


clusterFibersFuzzy <- function(data, bins, outlierCutoff, nodemax, kmax) {
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


clusterFibersKmeans <- function(data, bins, outlierCutoff, nodemax, kmax) {
  nodes <- list()
  centers <- list()

  for (i in seq_along(bins)) {
    bin <- bins[[i]]
    if (length(bin) < outlierCutoff)
      next

    k <- ceiling(length(bin)/nodemax)
    k <- min(k,kmax)

    if (k == 1) {
      nodes <- c(nodes, list(bin))
      centers[[i]] <- data[1,,drop=FALSE]-data[1,,drop=FALSE]
      next
    }

    km <- kmeans(data[bin,], centers=k)

    newNodes <- km$cluster %>%
      nodesFromMapping(bin = bin)

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

    # n <- length(nodes)
    d <- cdist(data[bin,],centers[[i]])
    cluster <- apply(d, 1, which.min)
    newNodes <- nodesFromMapping(mapping=cluster, bin=bin)
    nodes <- c(nodes, newNodes)
  }

  return(nodes)
}


nodesFromMapping <- function(mapping, bin) {
  n <- length(unique(mapping))
  nodes <- list()

  for (i in seq(n)) {
    nodes[[i]] <- bin[which(mapping == i)]
  }

  return(nodes)
}


getAdjList <- function(cluster, M) {
  nodes <- cluster$nodes
  edges <- cluster$edges

  n <- length(nodes)
  adjList <- vector(mode = "list", length = n)
  sizes <- sapply(nodes, length)

  for (i in seq(1,n-1)) {
    node <- nodes[[i]]
    for (j in seq(i+1,n)) {
      neighbor <- nodes[[j]]
      int <- intersect(node, neighbor)

      if (length(int) > M) {
        adjList[[i]] <- c(adjList[[i]], j)
        adjList[[j]] <- c(adjList[[j]], i)
      }
    }
  }

  for (k in seq_len(nrow(edges))) {
    i <- edges[k,1]
    j <- edges[k,2]

    adjList[[i]] <- c(adjList[[i]], j)
    adjList[[j]] <- c(adjList[[j]], i)
  }

  return(adjList)
}


getGraph <- function(cluster, data, M = 10) {
  adjList <- getAdjList(cluster, M)
  gr <- graph_from_adj_list(adjList, mode = "all")

  E(gr)$weight <- rep(1, length(E(gr))) # for Leiden
  print(nrow(cluster$edges) / length(E(gr)))

  return(gr)
}


getMedians <- function(data, nodes) {
  medians <- sapply(nodes, function(node) apply(data[node,], 2, median)) %>% t()
  return(medians)
}


