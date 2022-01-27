get_level_sets <- function(filter, nx = 10, ny = 10, overlap = 0.1) {
  # Split the dimensionally reduced space into overlapping bins
  # Assuming 2D. Setting ny=1 is effectively 1D.

  minx <- min(filter[,1])
  maxx <- max(filter[,1])

  miny <- min(filter[,2])
  maxy <- max(filter[,2])

  sizex <- (maxx-minx) / (nx*(overlap+1)+overlap)
  sizey <- (maxy-miny) / (ny*(overlap+1)+overlap)

  x0 <- minx + seq(0,nx-1)*sizex*(1+overlap)
  x1 <- x0 + sizex*(1+2*overlap)

  y0 <- miny + seq(0,ny-1)*sizey*(1+overlap)
  y1 <- y0 + sizey*(1+2*overlap)

  level_sets <- list(x0=x0, x1=x1, y0=y0, y1=y1)
  return(level_sets)
}

apply_level_sets <- function(filter, level_sets, outlier_cutoff = 0) {
  nx <- length(level_sets$x0)
  ny <- length(level_sets$y0)

  binsx <- lapply(seq(nx), function(i) which(filter[,1] >= level_sets$x0[i] &
                                             filter[,1] <= level_sets$x1[i]))
  binsy <- lapply(seq(ny), function(j) which(filter[,2] >= level_sets$y0[j] &
                                             filter[,2] <= level_sets$y1[j]))
  bins <- list()
  for (i in seq(nx))
    for (j in seq(ny)) {
      binid <- j + (i-1)*ny
      bins[[binid]] <- intersect(binsx[[i]], binsy[[j]])
    }

  counts <- sapply(bins, length)
  bins <- bins[which(counts >= outlier_cutoff)]

  return(bins)
}


cluster_level_sets <- function(data, bins, total_nodes, outlier_cutoff,
                             verbose, npc, method="kmeans") {
  if (method == "kmeans")
    return(cluster_kmeans(data, bins, total_nodes, outlier_cutoff))

  # Enter more methods here?

  stop("Invalid method chosen.")
}


cluster_kmeans <- function(data, bins, total_nodes, outlier_cutoff) {
  all_k <- get_k(data, bins, outlier_cutoff, total_nodes)
  nodes <- Map(function(bin, k) cluster_bin(bin, k, data, outlier_cutoff),
               bins, all_k$binK) %>% do.call(what=c)
  return(nodes)
}

cluster_bin <- function(bin, k, data, outlier_cutoff) {
  if(length(bin) < outlier_cutoff)
    return(list())

  if(k < 2)
    return(list(bin))

  km <- suppressWarnings(kmeans(data[bin,], centers=k, nstart = 10))
  new_nodes <- km$cluster %>%
    nodes_from_mapping(bin = bin, k = k)
  keep <- which(vapply(new_nodes, length, integer(1)) >= outlier_cutoff)

  return(new_nodes[keep])
}


get_k <- function(data, bins, outlier_cutoff, total_nodes,
                 fracSumSq = 0.75) {
  # Allocate nodes to bins, based on a combination of
  # bin size and bin variance

  binLength <- sapply(bins, length)
  notOut <- which(binLength>= outlier_cutoff)

  binLog <- pmax(log2(binLength/outlier_cutoff),0)
  kLog <- binLog * total_nodes / sum(binLog[notOut])

  binSumSq <- sapply(bins, function(bin) {
    if (length(bin) < 2)
      return(0)
    v <- var(data[bin,])
    return(sum(diag(v)))
  }) %>% sqrt()

  kSumSq <- binSumSq * total_nodes / sum(binSumSq[notOut])
  kSumSq <- pmin(kSumSq, floor(kLog))
  kSumSq <- kSumSq * total_nodes/sum(kSumSq)

  binK <- ceiling(fracSumSq*kSumSq + (1-fracSumSq)*kLog)
  binK <- pmin(binK, floor(binLength/outlier_cutoff))

  allK <- cbind(binLength, binLog, kLog, kSumSq, binK)
  return(data.frame(allK))
}


nodes_from_mapping <- function(mapping, bin, k) {
  lapply(seq_len(k), function(i) bin[which(mapping==i)])
}


get_graph <- function(nodes, node_stats, w, cutoff=0.01) {
  inters <- get_inters(nodes)
  inters <- inters/max(inters)
  inters[which(inters<cutoff)] <- 0

  gr <- graph_from_adjacency_matrix(inters,
                                    weighted = TRUE,
                                    mode="undirected")
  gr <- compute_weights(gr, node_stats$q50, w)
  return(gr)
}

get_inters <- function(nodes) {
  n <- length(nodes)
  inters <- matrix(0, nrow=n, ncol=n)

  for (i in seq(1,n-1)) {
    for (j in seq(i+1,n)) {
      val <- iou(nodes[[i]], nodes[[j]])
      inters[i,j] <- val
      inters[j,i] <- val
    }
  }
  return(inters)
}

iou <- function(node, neighbor) {
  int <- intersect(node,neighbor)
  uni <- union(node,neighbor)
  return(length(int)/length(uni))
}

compute_weights <- function(gr, medians, w, decay=4) {
  medians <- sweep(medians, MARGIN=2,STATS=w,FUN="*")
  d <- dist(medians) %>% as.matrix()

  edges <- as_edgelist(gr) %>% data.frame()
  names(edges) <- c("s", "t")
  d_rest <- d[as.matrix(edges[,c("s","t")])]
  new_weights <- exp(-decay*d_rest/mean(d_rest))

  E(gr)$weight <- new_weights
  return(gr)
}

get_stats <- function(data, nodes) {
  q25 <- sapply(nodes, function(node) apply(data[node,], 2, function(x) quantile(x,0.25))) %>% t()
  q50 <- sapply(nodes, function(node) apply(data[node,], 2, median)) %>% t()
  q75 <- sapply(nodes, function(node) apply(data[node,], 2, function(x) quantile(x,0.75))) %>% t()
  return(list(q25=q25,q50=q50,q75=q75))
}


plot_filter <- function(data, filter, path="") {
  sel <- sample(nrow(data), min(1e5, nrow(data)))
  df <- cbind(filter[sel,], data[sel,])

  for(param in colnames(data)) {
    message(param)
    g <- ggplot(df, aes(x=Comp.1, y=Comp.2)) +
      geom_point(aes_(color=as.name(param))) +
      scale_color_gradient(low = "black", high = "orange", name = param) +
      theme_bw()

    path_full <- paste0(path, "/", param, ".png")
    png(path_full, width=1200, height=900)
    plot(g)
    dev.off()
  }
}

apply_PCA <- function(data, pr, rank=2) {
  if(is.matrix(data)) {
    filter <- sweep(data, MARGIN=2, STATS=pr$center) %*% pr$loadings[,seq(rank)]
  } else {
    filter <- sweep(as.matrix(data), MARGIN=2, STATS=pr$center) %*% pr$loadings[,seq(rank)]
  }

  return(filter)
}
