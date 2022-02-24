leiden_clustering <- function(graph, resolution) {
  leid <- leiden.community(graph, resolution=resolution)
  clust <- as.factor(leid$membership)
  return(clust)
}


get_community_medians <- function(community, nodes, data) {
  lev <- levels(community)
  comm_events <- lapply(lev, function(l) Reduce(union, nodes[which(community==l)]))
  medians <- sapply(comm_events, function(comm) apply(data[comm,], 2, median)) %>% t()

  return(medians)
}

get_contingency_table <- function(mapping, samples) {
  base <- unname(table(samples)) %>% as.numeric()
  tab <- table(mapping, samples) %>% as.matrix()
  tab <- apply(tab, 1, function(row) row/base)
  if(length(base)==1)
    tab <- as.matrix(tab) %>% t()
  return(tab)
}


merge_communities <- function(modality) {
  mod_collapse <- apply(modality,1, function(x) paste(x,collapse=""))
  vals <- unique(mod_collapse)
  inv <- seq_along(vals)
  names(inv) <- vals
  group_index <- unname(inv[mod_collapse])
  return(group_index)
}

get_modality <- function(marker) {
  diana <- cluster::diana(marker)
  # diana <- cluster::agnes(marker)
  diana_cut <- cutree(as.hclust(diana), k=2)

  if (mean(marker[which(diana_cut==1)]) < mean(marker[which(diana_cut==2)])) {
    mod <- c("lo", "hi")
  } else {
    mod <- c("hi", "lo")
  }
  return(mod[diana_cut])
}

get_modality_thresholds <- function(medians) {
  node_modality <- apply(medians, 2, get_modality)
  markers <- colnames(medians)

  thresholds <- vapply(markers, function(marker) {
    neg <- which(node_modality[,marker]=="lo")
    pos <- which(node_modality[,marker]=="hi")
    lower <- max(medians[neg,marker])
    upper <- min(medians[pos,marker])
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

append_additional <- function(additional, modality, labels, matches) {
  for (marker in names(additional)) {
    pheno <- additional[[marker]]
    sel <- which(modality[,marker] == "hi")
    if (marker %in% names(defs)) {
      sel <- setdiff(sel,which(defs[matches,marker]== "hi"))
    }

    labels[sel] <- paste(labels[sel], pheno)
  }
  return(labels)
}


node_composition <- function(mapper, samples, scale=FALSE) {
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

  if (length(unq) == 1)
    counts <- matrix(counts, nrow=1)

  pctg <- apply(counts, 1, scale_vec) %>% t()

  if(scale)
    pctg <- apply(pctg, 2, scale_vec)

  return(data.frame(pctg))
}

scale_vec <- function(vec) {
  if (sum(vec) == 0)
    return(vec)
  return(vec/sum(vec))
}


make_unique_modality <- function(labels, modality) {
  lab <- unique(labels)
  unique_labels <- lapply(lab, function(label) {
    sel <- which(labels==label)
    if (length(sel)<=1)
      return(label)

    mod <- modality[sel,]
    diff <- which(!apply(mod,2,function(col) all(col==col[1])))
    mod <- mod[,diff,drop=FALSE]

    to_append <- apply(mod, 1, function(row) {
      paste(colnames(mod), unname(row), collapse=" ")
    })

    return(paste(label, to_append))
  }) %>% do.call(what=c)

  return(unique_labels)
}


assign_cells <- function(data, mapper, community=NULL) {
  mapping <- integer(nrow(data))
  distances <- numeric(nrow(data)) + Inf

  for (i in seq_along(mapper$nodes)) {
    node <- mapper$nodes[[i]]
    m <- mapper$node_stats$q50[i,]

    dist_new <- (t(data[node,])-m)^2 %>% apply(2, sum)
    closer <- which(dist_new < distances[node])

    distances[node[closer]] <- dist_new[closer]
    mapping[node[closer]] <- i
  }
  distances <- sqrt(distances)

  if (!is.null(community)) {
    # Some cells are unassigned
    n <- length(mapper$nodes)
    mapping[which(mapping==0)] <- n+1
    community_augmented <- c(as.character(community), "Unassigned")
    community_augmented <- factor(community_augmented,
                                  levels = c(seq_along(levels(community)), "Unassigned"))
    mapping <- community_augmented[mapping]
  }

  return(mapping)
}

