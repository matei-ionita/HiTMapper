leiden_clustering <- function(graph, resolution) {
  leid <- leiden.community(graph, resolution=resolution)
  clust <- as.factor(leid$membership)
  levels(clust) <- as.numeric(levels(clust)) + 1
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


get_modality_thresholds_gmm <- function(data, nodes, max_n=1e5) {
  markers <- colnames(data)

  max_n_node <- floor(min(max_n, nrow(data)/5)/length(nodes))
  sel <- lapply(nodes, function(node) {
    sample(node, min(length(node), max_n_node))
  }) %>% do.call(what = c)

  thresholds <- vapply(markers, get_modality_gmm, numeric(1), 
                       data=data, sel=sel)
  return(thresholds)
}


get_modality_gmm <- function(data, sel, marker) {
  # Mixture model with two components
  dat <- data[sel,marker]
  n_zero <- length(which(dat==0))
  dat <- dat[which(dat!=0)]
  
  fit <- Mclust(dat, G=2, modelNames = "V", verbose = FALSE)
  par <- fit$parameters
  
  m <- find_local_minimum_btn_means(par, data, marker)
  return(m)
}

find_local_minimum_btn_means <- function(par, data, marker) {
  kde <- bkde(data[,marker])
  
  mean1 <- min(par$mean[1], par$mean[2])
  mean2 <- max(par$mean[1], par$mean[2])
  idx_btn <- which(kde$x > mean1 & kde$x < mean2)
  
  idx_min <- idx_btn[which.min(kde$y[idx_btn])]
  
  return(kde$x[idx_min])
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

    labels[sel] <- paste(labels[sel], pheno, sep="_")
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
    mapping <- factor(community_augmented[mapping],
                      levels=c(levels(community), "Unassigned"))
  }
  
  return(mapping)
}


distinguish_communities <- function(mapper, data) {
  
  communities <- levels(mapper$community)
  
  unique_names <- communities %>% 
    str_split(fixed(".")) %>% sapply(function(x) x[1]) %>%
    unique() %>% sort()
  
  community_details <- lapply(unique_names, distinguish_name,
                              clustering=mapper$clustering, data=data,
                              communities=communities) %>% 
    do.call(what=rbind)
  
  new_community <- get_new_community(community_details)
  mapper$community <- new_community[as.character(mapper$community)] %>%
    unname() %>%
    as.factor()
  mapper$clustering <- new_community[as.character(mapper$clustering)] %>%
    unname() %>%
    as.factor()
  
  return(mapper)
}

get_new_community <- function(community_details) {
  new_community <- c(community_details$name, "Unassigned")
  names(new_community) <- c(community_details$community, "Unassigned")
  
  community_specific <- community_details %>%
    dplyr::filter(sensitivity > 0.85 & specificity > 0.9) %>%
    mutate(pheno=paste(name,pheno, sep="_"))
  
  new_community[community_specific$community] <- community_specific$pheno
  return(new_community)
}

distinguish_name <- function(name, clustering, data, communities) {
  communities_simple <- communities %>% 
    str_split(fixed(".")) %>% sapply(function(x) x[1])
  comm_name <- communities[which(communities_simple==name)]
  
  if (length(comm_name) == 1)
    return(tibble(name=name, community=name, sensitivity=0, specificity=0, 
                  pheno="", desc="Only 1 community with this name."))
  
  lapply(comm_name, distinguish_community, name=name, data=data, 
         clustering=clustering) %>% 
    do.call(what=rbind)
}

distinguish_community <- function(community, name, clustering, data) {
  model <- select_features(data, clustering,
                           pop1=name,
                           pop2=paste0(community, "$"),
                           max_markers = 3)
  
  if(length(model$coeff) < 2)
    return(tibble(name=name, 
                  community=community, 
                  sensitivity=0, specificity=0,
                  pheno="", desc="Found no separating markers."))
  
  coeff_name <- names(model$coeff[-1])
  coeff_vals <- unname(model$coeff[-1]) %>% round(2)
  ord <- order(abs(coeff_vals), decreasing = TRUE)
  
  tab <- table(model$pred$class, model$pred$prediction) %>% as.matrix()
  
  coeff <- matrix(c(coeff_name[ord], coeff_vals[ord]),ncol=2)
  desc <- coeff %>% t() %>% as.character() %>% paste(collapse=" ")
  
  coeff_sign <- matrix(c(coeff_name[ord], 
                         if_else(coeff_vals[ord]>0, "hi", "lo")), ncol=2)
  pheno <- coeff_sign %>% t() %>% as.character() %>% paste(collapse="_")
  
  tib <- tibble(name=name, community=community,
                sensitivity = tab[2,2] / (tab[2,2] + tab[2,1]),
                specificity = tab[1,1] / (tab[1,1] + tab[1,2]),
                pheno = pheno, desc=desc)
  return(tib)
}



select_features <- function(data, labels, pop1, pop2,
                            subsample=1e6, max_markers=5,
                            nlambda=50) {
  sel <- grep(paste0(pop1,"|",pop2), labels)
  if (!is.null(subsample) & subsample < length(sel))
    sel <- sample(sel, subsample)
  
  x <- data[sel,] %>% scale()
  y <- if_else(grepl(pop2, labels[sel]), pop2, pop1)
  w <- if_else(y==pop1, 1, length(which(y==pop1))/length(which(y==pop2)))
  
  model <- glmnet(x, y, weights=w, family="binomial", alpha=1, 
                  nlambda=nlambda, pmax = max_markers, 
                  lambda.min.ratio = 1e-3) %>% suppressWarnings()
  
  n <- ncol(coef(model))
  coeff <- coef(model)[,n]
  coeff <- coeff[coeff!=0]
  
  if (length(coeff) < 2) {
    return(list(coeff=coeff))
  }
  
  prob <- predict(model, newx=x %>% scale(), type="response")[,n]
  if (mean(prob[which(y==pop1)]) > mean(prob[which(y==pop2)]) ){
    pophi <- pop1
    poplo <- pop2
  } else {
    pophi <- pop2
    poplo <- pop1
  }
  
  thresh <- get_optimal_threshold(class=y, prob=prob)
  pred <- if_else(prob > thresh, pophi, poplo)
  
  tib <- tibble(class=y,
                prob=prob,
                prediction=pred)
  
  return(list(pred = tib, coeff = coeff, sel = sel, thresh=thresh))
}


get_optimal_threshold <- function(class, prob) {
  un <- unique(class)
  class1 <- which(class==un[2])
  class2 <- which(class==un[1])
  
  kde1 <- bkde(prob[class1], range.x = c(0.01,0.99))
  kde2 <- bkde(prob[class2], range.x = c(0.01,0.99))
  
  diff <- kde1$y-kde2$y * length(class2)/length(class1)
  signs <- sign(kde1$y-kde2$y * length(class2)/length(class1))
  
  M1 <- which.max(kde1$y)
  M2 <- which.max(kde2$y)
  
  
  if (M1 < M2) {
    interv <- seq(M1,M2)
    thresh <- kde1$x[interv[which.min(abs(diff[interv]))]]
  } else {
    interv <- seq(M2,M1)
    thresh <- kde1$x[interv[which.min(abs(diff[interv]))]]
  }
  
  return(thresh)
}







