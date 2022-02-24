get_level_sets <- function(filter, grid_size=c(10,10), overlap = 0.1) {
  # Split the dimensionally reduced space into overlapping bins

  n <- length(grid_size)
  extrema <- lapply(seq_len(n), function(i) c(min(filter[,i]),max(filter[,i])))

  widths <- Map(function(extrem, dim_size) {
    (extrem[2]-extrem[1]) / (dim_size*(overlap+1)+overlap)
  }, extrema, grid_size)

  boundaries <- Map(function(extrem, dim_size, width) {
    lower <- extrem[1] + seq(0,dim_size-1)*width*(1+overlap)
    upper <- lower + width * (1+2*overlap)
    return(list(lower=lower,upper=upper))
  }, extrema, grid_size, widths)

  return(boundaries)
}

apply_level_sets <- function(filter, boundaries,
                               outlier_cutoff=0, max_bin_size=nrow(filter)/10) {
  grid_size <- vapply(boundaries, function(b) length(b$lower), integer(1))
  bins_1D <- Map(function(dim, boundaries_dim) get_bins_1D(filter, dim, boundaries_dim),
                 seq_along(grid_size), boundaries)

  bins <- Reduce(intersect_bins, bins_1D) %>%
    lapply(split_bins, max_bin_size=max_bin_size) %>%
    do.call(what=c)

  counts <- vapply(bins, length, integer(1))
  bins <- bins[which(counts >= outlier_cutoff)]
  return(bins)
}

split_bins <- function(bin, max_bin_size) {
  if (length(bin) < max_bin_size)
    return(list(bin))

  n <- ceiling(length(bin)/max_bin_size)
  ind <- sample(seq_len(n), length(bin), replace=TRUE)

  split_bin <- lapply(seq_len(n), function(i) bin[which(ind==i)])
  return(split_bin)
}

get_bins_1D <- function(filter, dim, boundaries_dim) {
  dim_size <- length(boundaries_dim$lower)
  lapply(seq(dim_size), function(i) which(filter[,dim] >= boundaries_dim$lower[i] &
                                          filter[,dim] <= boundaries_dim$upper[i]))
}

intersect_bins <- function(bins_1, bins_2) {
  lapply(bins_1, function(bin_1) {
    lapply(bins_2, function(bin_2) intersect(bin_1,bin_2))
  }) %>% do.call(what=c)
}



cluster_level_sets <- function(data, bins, total_nodes, outlier_cutoff,
                             verbose, method="kmeans") {
  if (method == "kmeans")
    return(cluster_kmeans(data, bins, total_nodes, outlier_cutoff))

  # Enter more methods here?

  stop("Invalid method chosen.")
}


cluster_kmeans <- function(data, bins, total_nodes, outlier_cutoff) {
  all_k <- get_k(data, bins, outlier_cutoff, total_nodes)
  nodes <- Map(function(bin, k) cluster_bin_kmeans(bin, k, data, outlier_cutoff),
               bins, all_k) %>% do.call(what=c)
  return(nodes)
}


cluster_bin_kmeans <- function(bin, k, data, outlier_cutoff) {
  if(length(bin) < outlier_cutoff)
    return(list())

  if(k < 2)
    return(list(bin))

  km <- suppressWarnings(kmeans(data[bin,], centers=k, nstart = 10))
  new_nodes <- km$cluster %>%
    nodes_from_mapping(bin = bin, lev=seq_len(k))
  keep <- which(vapply(new_nodes, length, integer(1)) >= outlier_cutoff)

  return(new_nodes[keep])
}


get_k <- function(data, bins, outlier_cutoff, total_nodes,
                 fracSumSq = 0.75) {
  # Allocate nodes to bins, based on a combination of
  # bin size and bin variance

  binLength <- sapply(bins, length)

  binLog <- pmax(log2(binLength),0)
  kSize <- binLog * total_nodes / sum(binLog)

  binSumSq <- sapply(bins, function(bin) {
    if (length(bin) < 2)
      return(0)
    v <- var(data[bin,])
    return(sum(diag(v)))
  }) %>% sqrt()

  kSumSq <- binSumSq * total_nodes / sum(binSumSq)
  kSumSq <- pmin(kSumSq, floor(kSize))
  kSumSq <- kSumSq * total_nodes/sum(kSumSq)

  k <- ceiling(fracSumSq*kSumSq + (1-fracSumSq)*kSize)
  k <- pmin(k, floor(binLength/outlier_cutoff))
  return(k)
}


nodes_from_mapping <- function(mapping, bin, lev) {
  lapply(lev, function(i) bin[which(mapping==i)])
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
