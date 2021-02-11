getBins <- function(filter, nx = 10, ny = 10, overlap = 0.1, outlierCutoff = 10) {
  # Split the dimensionally reduced space into overlapping bins
  # Assuming 2D. Setting ny=1 is effectively 1D.

  minx <- min(filter[,1])
  maxx <- max(filter[,1])

  miny <- min(filter[,2])
  maxy <- max(filter[,2])

  sizex <- (maxx-minx) / (nx*(overlap+1)+overlap)
  sizey <- (maxy-miny) / (ny*(overlap+1)+overlap)

  startx <- minx + seq(0,9)*sizex*(1+overlap)
  stopx <- startx + sizex*(1+2*overlap)

  starty <- miny + seq(0,9)*sizey*(1+overlap)
  stopy <- starty + sizey*(1+2*overlap)

  binsx <- lapply(seq(nx), function(i) which(filter[,1] >= startx[i] & filter[,1] <= stopx[i]))
  binsy <- lapply(seq(ny), function(j) which(filter[,2] >= starty[j] & filter[,2] <= stopy[j]))

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

clusterFibers <- function(data, bins, method = "kmeans", outlierCutoff = 50, nodemax = 1000) {
  if (method == "kmeans")
    return(clusterFibersKmeans(data, bins, outlierCutoff, nodemax))

  # Enter more methods here

  stop("Invalid method chosen.")
}

clusterFibersKmeans <- function(data, bins, outlierCutoff, nodemax) {
  nodes <- list()

  for (bin in bins) {
    if (length(bin) < outlierCutoff)
      next

    k <- ceiling(length(bin)/nodemax)
    k <- min(k,20)

    if (k == 1) {
      nodes <- c(nodes, list(bin))
      next
    }

    km <- kmeans(data[bin,], centers=k)
    newNodes <- km$cluster %>%
      nodesFromMapping(bin = bin)

    keep <- which(sapply(newNodes, length) >= outlierCutoff)
    nodes <- c(nodes, newNodes[keep])
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

getAdjList <- function(nodes, M) {
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

  return(adjList)
}


getGraph <- function(bins, data, M = 10) {
  adjList <- getAdjList(bins, M)
  gr <- graph_from_adj_list(adjList, mode = "all")
  return(gr)
}


getMedians <- function(data, nodes) {
  medians <- sapply(nodes, function(node) apply(data[node,], 2, median)) %>% t()
  return(medians)
}


