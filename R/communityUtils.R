leiden_clustering <- function(gr, resolution, out) {
  leid <- leiden.community(gr, resolution=resolution)
  clust <- as.factor(leid$membership)
  # change 0-based to 1-based
  levels(clust) <- as.numeric(levels(clust)) + 1
  
  return(as.factor(clust))
  
}


get_weights <- function(tab, sim, min_n) {
  union_size <- outer(tab, tab, FUN="+")
  sim <- sim/union_size
  
  outliers <- which(tab <= min_n)
  sim[outliers,] <- sim[outliers,] / 10
  sim[setdiff(seq(nrow(sim)),outliers),outliers] <- 
    sim[setdiff(seq(nrow(sim)),outliers),outliers] / 10
  return(sim)
}


get_graph_sim <- function(sim) {
  gr <- graph_from_adjacency_matrix(sim, mode="undirected", weighted = TRUE)
  # remove edges with tiny weights, to help plotting
  gr <- delete_edges(gr, which(E(gr)$weight < 1e-3))
  return(gr)
}


parse_communities <- function(mapper, data) {
  # compute community centroids and use them for labeling
  mapper$community_medians <- compute_centroids(data, mapper$clustering,
                                                length(levels(mapper$community)))
  colnames(mapper$community_medians) <- colnames(data)
  
  if (!is.null(mapper$defs)) {
    mapper <- label_communities(mapper, mapper$defs)
  }
  
  return(mapper)
}


get_contingency_table <- function(mapping, samples) {
  # tabulate cluster percentages for all samples
  base <- unname(table(samples)) %>% as.numeric()
  tab <- table(mapping, samples) %>% as.matrix()
  tab <- apply(tab, 1, function(row) row/base)
  if(length(base)==1)
    tab <- as.matrix(tab) %>% t()
  return(tab)
}


get_modality_pos <- function(marker) {
  pos <- which(marker>=0)
  diana <- cluster::diana(marker[pos])
  diana_cut <- cutree(as.hclust(diana), k=2)

  if (mean(marker[pos][which(diana_cut==1)]) < 
      mean(marker[pos][which(diana_cut==2)])) {
    mod <- c("lo", "hi")
  } else {
    mod <- c("hi", "lo")
  }
  
  result <- rep("lo", length(marker))
  result[pos] <- mod[diana_cut]
  return(result)
}


get_modality_thresholds_hier <- function(centroids) {
  node_modality <- apply(centroids, 2, get_modality_pos)
  markers <- colnames(centroids)
  
  thresholds <- vapply(markers, function(marker) {
    neg <- which(node_modality[,marker]=="lo")
    pos <- which(node_modality[,marker]=="hi")
    lower <- max(centroids[neg,marker])
    upper <- min(centroids[pos,marker])
    return((lower+upper)/2)
  }, numeric(1))
  return(thresholds)
}


match_defs <- function(defs, modality) {
  ind <- integer(nrow(modality))

  for (i in seq(nrow(defs))) {
    markers <- names(defs)[which(defs[i,] %in% c("hi", "lo"))]
    match <- apply(modality[,markers,drop=FALSE], 1, function(x) {
      def <- as.character(defs[i,markers])
      chx <- as.character(x)
      all.equal(chx,def)
    })
    sel <- which(match == "TRUE")
    ind[sel] <- i
  }

  ind[which(ind==0)] <- nrow(defs)+1
  return(ind)
}

