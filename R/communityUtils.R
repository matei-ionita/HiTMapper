calibrate_weights <- function(mapper) {
  medians <- mapper$nodeStats$q50
  d <- dist(medians) %>% as.matrix()

  edges <- as_edgelist(mapper$gr) %>% data.frame()
  names(edges) <- c("s", "t")
  edges$weight <- E(mapper$gr)$weight
  d_rest <- d[as.matrix(edges[,c("s","t")])]
  new_weights <- exp(-3*d_rest/mean(d_rest))

  E(mapper$gr)$weight <- new_weights
  return(mapper)
}


get_labels <- function(mapper, defs, additional=list()) {
  pos <- apply(mapper$community_medians, 2, function(v) {
    dv <- cluster::diana(v)
    dv2 <- cutree(as.hclust(dv), k=2)

    if (mean(v[which(dv2==1)]) < mean(v[which(dv2==2)])) {
      lab <- c("lo", "hi")
    } else {
      lab <- c("hi", "lo")
    }
    return(lab[dv2])
  })
  # pos[,"CD159c"] <- "lo"

  labels <- rep("Other", nrow(pos))
  ind <- integer(nrow(pos))
  defs <- mutate_if(defs, is.character, str_replace_all, pattern="mid", replacement="hi")
  for (i in seq(nrow(defs))) {
    markers <- names(defs)[which(defs[i,] %in% c("hi", "lo"))]
    sel <- which(apply(pos[,markers,drop=FALSE], 1, function(x) {
      def <- as.character(defs[i,markers])
      chx <- as.character(x)
      all.equal(chx,def)
    }) == "TRUE")
    labels[sel] <- defs$Phenotype[i]
    ind[sel] <- i
  }

  for (marker in names(additional)) {
    pheno <- additional[[marker]]
    sel <- which(pos[,marker] == "hi")
    if (marker %in% names(defs)) {
      sel <- intersect(sel,which(defs[ind,marker]!= "hi"))
    }

    labels[sel] <- paste(labels[sel], pheno)
  }

  posvec <- apply(pos,1, function(x) paste(x,collapse=""))
  for (un in unique(labels)) {
    sel <- which(labels==un)
    if(length(unique(posvec[sel]))>1) {
      df <- tibble(labels=labels[sel],
                   phenos=posvec[sel]) %>%
        group_by(phenos)
      labels[sel] <- paste(labels[sel], group_indices(df))
    }
  }

  lab_node <- labels[mapper$community]
  return(as.factor(lab_node))
}
