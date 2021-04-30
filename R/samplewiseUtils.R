
get_mappers <- function(path, files, ignore, kNodes, nx, ny, overlap) {
  lapply(files, function(file) {
    message(file)

    inpath <- paste0(path, file)
    ff <- read.FCS(inpath, transformation=FALSE, truncate_max_range=FALSE)
    keep <- setdiff(flowCore::colnames(ff), ignore)
    data <- exprs(ff)[,keep]

    mapper <- HiTMapper(data, kNodes=kNodes,
                        nx=nx,ny=ny, overlap=overlap)
    mapper <- pruneEdges(mapper, cutoff = 0.01)

    file_base <- str_remove(file, ".fcs")
    n <- length(V(mapper$gr))
    vertex_names <- paste(file_base, seq(n), sep="_")
    V(mapper$gr)$name <- vertex_names

    return(mapper)
  })
}



get_meta_gr <- function(mappers, a=0.25, b=0.1, w=0.5) {
  edge_list <- matrix(nrow=0,ncol=2)
  weights <- c()
  n <- length(mappers)

  for (i in seq(n)) {
    for (j in seq(i,n)) {
      if (i==j) {
        gr <- mappers[[i]]$gr
        edge_list <- rbind(edge_list, get.edgelist(gr))
        weights <- c(weights, E(gr)$weight)
      } else {
        name_i <- V(mappers[[i]]$gr)$name
        name_j <- V(mappers[[j]]$gr)$name

        mi <- mappers[[i]]$nodeStats$q50
        mj <- mappers[[j]]$nodeStats$q50
        dis <- rdist::cdist(mi, mj)

        pairs1 <- cbind(seq(nrow(mi)), apply(dis,1,which.min))
        pairs2 <- cbind(apply(dis,2,which.min), seq(nrow(mj)))
        pairs <- rbind(pairs1,pairs2) %>% unique()
        # ker <- exp(-a*dis^2)
        #
        # pairs <- which(ker > b, arr.ind = TRUE)
        el <- cbind(name_i[pairs[,1]], name_j[pairs[,2]])
        edge_list <- rbind(edge_list, el)
        weights <- c(weights, rep(w, nrow(pairs)))
        # weights <- c(weights, ker[pairs])
      }
    }
  }
  gr <- graph_from_edgelist(edge_list, directed=FALSE)
  E(gr)$weight <- weights

  perm <- get_perm(mappers,gr)
  gr <- permute(gr,perm)

  sam <- V(gr)$name %>%
    str_split(pattern="_") %>%
    sapply(function(x) x[1])
  V(gr)$sample <- sam

  return(gr)
}



get_perm <- function(mappers, gr) {
  n <- length(mappers)

  nms_old <- lapply(seq(n), function(i) V(mappers[[i]]$gr)$name) %>% do.call(what=c)
  nms <- V(gr)$name
  o <- order(nms)
  perm <- order(nms_old)[order(o)]

  return(perm)
}



get_mapping <- function(path, files, mappers, ignore) {
  mapping <- lapply(seq_along(files), function(i) {
    inpath <- paste0(path, files[i])
    ff <- read.FCS(inpath, transformation=FALSE, truncate_max_range=FALSE)
    keep <- setdiff(colnames(ff), ignore)
    data <- exprs(ff)[,keep]

    mapping <- assignCells(data, mappers[[i]], V(mappers[[i]]$gr)$name)

    return(as.character(mapping$mapping))
  }) %>% do.call(what=c)
}


get_sample_mapping <- function(gr, Nevents) {
  fns <- V(gr)$name %>%
    str_split(pattern="_") %>%
    sapply(function(x) x[1]) %>%
    unique()
  return(rep.int(fns, times=Nevents))
}


get_contingency_table <- function(mapping, samples) {
  base <- unname(table(samples))
  tab <- table(mapping, samples) %>% as.matrix()
  tab <- apply(tab, 1, function(row) row/base)

  return(tab)
}


get_community_medians <- function(community, mappers) {
  size <- lapply(mappers, function(mapper) sapply(mapper$nodes, length)) %>%
    do.call(what=c)

  node_medians <- lapply(mappers, function(mapper) mapper$nodeStats$q50) %>%
    do.call(what=rbind)

  community_medians <- sapply(levels(community), function(lev) {
    comm <- which(community == lev)
    med <- apply(node_medians[comm,], 2, function(x) weightedMedian(x,w=size[comm]))
  }) %>% t()

  return(community_medians)
}



get_n_events <- function(path, files) {
  sapply(files, function(file) {
    inpath <- paste0(path, file)
    ff <- read.FCS(inpath, transformation=FALSE, truncate_max_range=FALSE)
    return(nrow(ff))
  })
}

pool_data <- function(path, files, ignore) {
  lapply(files, function(file) {
    inpath <- paste0(path, file)
    ff <- read.FCS(inpath, transformation=FALSE, truncate_max_range=FALSE)
    keep <- setdiff(flowCore::colnames(ff), ignore)
    data <- exprs(ff)[,keep]
    return(data)
  }) %>% do.call(what=rbind)
}



get_sample_mapping_pooled <- function(path, files) {
  lapply(files, function(file) {
    inpath <- paste0(path, file)
    basename <- str_remove(file, ".fcs")
    ff <- read.FCS(inpath, transformation=FALSE, truncate_max_range=FALSE)
    return(rep(basename, nrow(ff)))
  }) %>% do.call(what=c)
}


get_community_medians_pooled <- function(community, mapper, data) {
  lev <- levels(community)
  commEvents <- lapply(lev, function(l) Reduce(union, mapper$nodes[which(community==l)]))
  commStats <- getStats(data, commEvents)

  return(commStats$q50)
}


