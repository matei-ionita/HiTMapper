leiden_clustering <- function(gr, resolution=1) {
  leid <- leiden.community(gr, resolution=resolution)
  clust <- as.factor(leid$membership)
  # change 0-based to 1-based
  levels(clust) <- as.numeric(levels(clust)) + 1
  return(clust)
}


build_graph <- function(data, mapping, m=10, k=20) {
  # sel_list <- lapply(unique(mapping), function(i) {
  #   mi <- which(mapping==i)
  #   sample(mi, min(m, length(mi)))
  # })
  # sel <- unlist(sel_list)
  chosen <- sample_cells(mapping, unique(mapping), m)
  dat <- data[chosen,]
  node_of_origin <- mapping[chosen]
  # node_of_origin <- rep(unique(mapping), sapply(sel_list, length))

  knn <- hnsw_knn(dat, k=k)
  edgelist <- get_edgelist(knn$idx)
  
  df_edgelist <- tibble(ind1 = edgelist[,1],
                        ind2 = edgelist[,2],
                        weight = edgelist[,3]) %>%
    mutate(node1 = node_of_origin[ind1],
           node2 = node_of_origin[ind2]) %>%
    group_by(node1, node2) %>%
    summarise(weight=sum(weight)) %>%
    dplyr::arrange(node2) %>%
    tidyr::pivot_wider(names_from="node2", values_from="weight", values_fill=0) %>%
    dplyr::arrange(node1)
  sim <- as.matrix(df_edgelist)[,-1]
  dia <- sqrt(diag(sim))
  sim <- t(t(sim/dia)/dia)
  sim <- sim + t(sim) - diag(1, nrow=nrow(sim))

  gr <- graph_from_adjacency_matrix(sim, mode="undirected", weighted=TRUE)
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

