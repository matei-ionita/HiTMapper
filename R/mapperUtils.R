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
                               min_node_size=0, max_bin_size=nrow(filter)/10) {
  grid_size <- vapply(boundaries, function(b) length(b$lower), integer(1))
  bins_1D <- Map(function(dim, boundaries_dim) get_bins_1D(filter, dim, boundaries_dim),
                 seq_along(grid_size), boundaries)

  bins <- Reduce(intersect_bins, bins_1D) %>%
    lapply(split_bins, max_bin_size=max_bin_size) %>%
    do.call(what=c)

  counts <- vapply(bins, length, integer(1))
  bins <- bins[which(counts >= min_node_size)]
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


cluster_level_sets <- function(data, bins, total_nodes, min_node_size,
                               n_passes, verbose, method) {
  if (method == "kmeans")
    return(cluster_kmeans(data, bins, total_nodes, min_node_size))

  if (method == "som")
    return(cluster_som(data, bins, total_nodes, min_node_size, n_passes))

  stop("Invalid method chosen.")
}


cluster_kmeans <- function(data, bins, total_nodes, min_node_size) {
  all_k <- get_k(data, bins, min_node_size, total_nodes)
  nodes <- Map(function(bin, k) cluster_bin_kmeans(bin, k, data, min_node_size),
               bins, all_k) %>% do.call(what=c)
  return(nodes)
}

cluster_som <- function(data, bins, total_nodes, min_node_size, n_passes) {
  all_k <- get_k(data, bins, min_node_size, total_nodes)
  nodes <- Map(function(bin, k) cluster_bin_som(bin, k, data, min_node_size, n_passes),
               bins, all_k) %>% do.call(what=c)
  return(nodes)
}


cluster_bin_kmeans <- function(bin, k, data, min_node_size) {
  if(length(bin) < min_node_size)
    return(list())

  if(k < 2)
    return(list(bin))

  km <- suppressWarnings(kmeans(data[bin,], centers=k, nstart = 10))
  new_nodes <- km$cluster %>%
    nodes_from_mapping(bin = bin, lev=seq_len(k))
  keep <- which(vapply(new_nodes, length, integer(1)) >= min_node_size)

  return(new_nodes[keep])
}


cluster_bin_som <- function(bin, k, data, min_node_size, n_passes) {
  
  if(length(bin) < 5)
    return(list())

  if(length(bin) < min_node_size)
    return(list(bin))

  if(k < 2)
    return(list(bin))

  nx <- floor(sqrt(k))
  ny <- ceiling(sqrt(k))
  som <- som(data[bin,], nx, ny, n_passes)
  new_nodes <- som$mapping %>%
    nodes_from_mapping(bin = bin, lev=seq_len(nx*ny))
  keep <- which(vapply(new_nodes, length, integer(1)) >= 5)

  return(new_nodes[keep])
}


get_k <- function(data, bins, min_node_size, total_nodes,
                         frac_sum_sq = 0.75) {
  # Allocate nodes to bins, based on a combination of
  # bin size and bin variance
  
  bin_length <- sapply(bins, length)
  bin_log_size <- pmax(log2(bin_length),0)
  weight_size <- bin_log_size / sum(bin_log_size)
  
  bin_sum_sq <- sapply(bins, function(bin) {
    if (length(bin) < 2)
      return(0)
    v <- var(data[bin,])
    return(sum(diag(v)))
  }) %>% sqrt()
  weight_sum_sq <- bin_sum_sq / sum(bin_sum_sq)
  
  weight <- frac_sum_sq*weight_sum_sq + (1-frac_sum_sq)*weight_size
  max_k_bin <- ceiling(bin_length/min_node_size)
  
  k <- allocate_with_constraints(total_nodes, max_k_bin, weight, 
                                 alloc=rep(0, length(bins)))
  return(k)
}

allocate_with_constraints <- function(total_nodes, max_k_bin, weight, alloc,
                                      depth=1) {
  if (sum(alloc)>=total_nodes | depth > 10)
    return(round(alloc))
  
  to_alloc <- total_nodes - sum(alloc)
  not_at_max <- which(alloc!=max_k_bin)
  renorm <- sum(weight[not_at_max])
  
  alloc_new <- pmin(weight*to_alloc/renorm, max_k_bin-alloc)

  return(allocate_with_constraints(total_nodes, max_k_bin, weight, 
                                   alloc+alloc_new, depth=depth+1))
}


nodes_from_mapping <- function(mapping, bin, lev) {
  lapply(lev, function(i) bin[which(mapping==i)])
}


get_graph_nn <- function (medians, k, thresholds, force_sep)
{
  d <- dist(medians) %>% as.matrix()
  nn <- lapply(seq(nrow(medians)), function(i) {
    setdiff(order(d[i,])[seq(k+1)], i)
  })

  jac <- lapply(nn, function(node) {
    lapply(node, function(nbr) {
      inter <- intersect(node, nn[[nbr]])
      denom <- length(node) + length(nn[[nbr]]) - length(inter)
      numer <- length(inter)
      return(numer / denom)
    }) %>% do.call(what=c)
  })
  
  jac <- force_smaller_weights_on_thresholds(jac, nn, medians, k,
                                             thresholds, force_sep)
  
  gr <- graph_from_adj_list(nn)
  w <- do.call(what=c, args=jac)
  ## most recent version of ggraph only works with
  ## positive weights in FR algorithm
  w[which(w==0)] <- 1e-5
  E(gr)$weight <- w

  return(gr)
}

force_smaller_weights_on_thresholds <- function(jac, nn, medians, k,
                                                thresholds, force_sep) {
  for (i in seq_along(nn)) {
    for (j in seq(k)) {
      nbr <- nn[[i]][j]
      flag <- FALSE
      for (marker in force_sep) {
        if ( (medians[i,marker] - thresholds[marker]) *
             (medians[nbr,marker] - thresholds[marker]) < 0) {
          flag <- TRUE
        }
      }
      if (flag) {
        jac[[i]][j] <- jac[[i]][j]/10
      }
    }
  }
  
  return(jac)
}

get_stats <- function(data, nodes) {
  q50 <- sapply(nodes, function(node) apply(data[node,], 2, median)) %>% t()
  return(list(q50=q50))
}


get_filter <- function(data, rank=2, scale=FALSE) {
  if(!is.matrix(data) & !is.data.frame(data))
    stop("Please enter your data in matrix or data.frame format.")
  
  if (scale) {
    cov <- cor(data)
  } else {
    cov <- cov(data)
  }
  pr <- princomp(covmat=cov, scores=FALSE)
  pr$center <- apply(data, 2, mean)
  filter <- apply_PCA(data, pr, rank)
  
  return(filter)
}


apply_PCA <- function(data, pr, rank=2) {
  if(is.matrix(data)) {
    filter <- sweep(data, MARGIN=2, STATS=pr$center) %*% pr$loadings[,seq(rank)]
  } else {
    filter <- sweep(as.matrix(data), MARGIN=2, STATS=pr$center) %*% pr$loadings[,seq(rank)]
  }

  return(filter)
}
