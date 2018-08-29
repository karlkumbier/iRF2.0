generalizedRIT <- function(rf=NULL, x=NULL, rforest=NULL,
                           weights=rep(1, nrow(x)),
                           varnames.grp=colnames(x),
                           rit.param=list(depth=5,
                                          ntree=500,
                                          nchild=2,
                                          class.id=1,
                                          min.nd=1,
                                          class.cut=NULL),
                           get.prevalence=FALSE,
                           int.sign=TRUE,
                           int.subs=TRUE,
                           out.file=NULL,
                           n.core=1) {

  out <- list()
  if (!is.null(rf) & !is.null(x)| !is.null(rforest))
    stop('Supply random forest or read forest output')
  if (is.null(rit.param$class.cut) & is.null(rit.param$class.id))
    stop('Supply class.id (classification) or class.cut (regression)')

  # Set feature names for grouping interactions
  if (is.null(varnames.grp)) varnames.grp <- as.character(1:ncol(x))
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)
  
  # Read RF object to extract decision path metadata
  if (is.null(rforest)) {
    rforest <- readForest(rf, x=x,
                          return.node.feature=TRUE,
                          return.node.obs=TRUE,
                          varnames.grp=varnames.grp,
                          get.split=TRUE,
                          n.core=n.core)
  }

  # Collapse node feature matrix for unsigned iRF
  if (!int.sign) {
    rforest$node.feature <- rforest$node.feature[,1:p] +
    rforest$node.feature[,(p + 1):(2 * p)]
  }

  if (!is.null(out.file)) save(file=out.file, rforest)

  # Select class specific leaf nodes
  rfpred <- rforest$tree.info$prediction
  if (is.null(rit.param$class.cut))
    select.id <- rfpred == rit.param$class.id + 1
  else
    select.id <- rfpred > rit.param$class.cut

  # Run RIT on leaf nodes of selected class to find
  # prevalent interactions
  if (sum(select.id) < 2) {
    return(character(0))
  } else {
    
    wt <- Matrix::colSums(t(rforest$node.obs) * weights)
    out$int <- runRIT(subsetReadForest(rforest, select.id),
                      weights=wt[select.id],
                      rit.param=rit.param,
                      int.subs=int.subs,
                      n.core=n.core)
    
    if (is.null(out$int)) return(NULL)
    int.names <- nameInts(out$int, varnames.unq, signed=int.sign)

    if (get.prevalence) {
      out$prev$i1 <- mclapply(out$int, prevalence,
                              nf=rforest$node.feature[select.id,],
                              wt=wt[select.id],
                              mc.cores=n.core)
      out$prev$i1 <- unlist(out$prev$i1)

      out$prev$i0 <- mclapply(out$int, prevalence,
                              nf=rforest$node.feature[!select.id,],
                              wt=wt[!select.id],
                              mc.cores=n.core)
      out$prev$i0 <- unlist(out$prev$i0)

      names(out$prev$i1) <- int.names
      names(out$prev$i0) <- int.names
    }

    out$int <- int.names
  }
  return(out)
}

runRIT <- function(rforest, weights, rit.param, int.subs=TRUE, n.core=1) {

  # remove nodes below specified size threshold
  id.rm <- rforest$tree.info$size.node < rit.param$min.nd
  if (all(id.rm)) {
    warning(paste('No nodes with greater than ', rit.param$min.nd,
                  'observations. Using all nodes.'))
    id.rm <- rep(FALSE, length(id.rm))
  }
  
  rforest <- subsetReadForest(rforest, !id.rm)
  wt <- weights[!id.rm]
  interactions <- RIT(rforest$node.feature[wt > 0,],
                      weights=wt[wt > 0],
                      depth=rit.param$depth,
                      n_trees=rit.param$ntree,
                      branch=rit.param$nchild,
                      output_list=TRUE,
                      n_cores=n.core)$Interaction

  if (int.subs) {
    interactions <- lapply(interactions, intSubsets)
    interactions <- unlist(interactions, recursive=FALSE)
    interactions <- unique(interactions)
  }

  return(interactions)
}

intSubsets <- function(int) {
  # Generate all lower order subsets of specified interaction
  int.order <- length(int)
  int.subs <- lapply(1:int.order, combn, x=int, simplify=FALSE)
  int.subs <- unlist(int.subs, recursive=FALSE)
  return(int.subs)
}


nameInts <- function(ints, varnames, signed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by
  # variable names. Naming convention for an interaction is:
  #   <variable1(sign)>_ <variable2(sign)>_...
  
  varnames <- unique(varnames)
  p <- length(varnames)
  if (signed) 
    ints.signs <- lapply(ints, function(z) ifelse(z > p, '+', '-'))
  else 
    ints.signs <- ''
  
  # Adjust indexing to match varnames
  ints <- lapply(ints, function(z) z %% p + p * (z == p | z == 2 * p))
  ints.name <- mapply(function(i, s) paste0(varnames[i], s),
                      ints, ints.signs, SIMPLIFY=FALSE)
  ints.name <- sapply(ints.name, function(z) paste(sort(z), collapse='_'))
  return(ints.name)
}

summarizePrev <- function(prev) {
  # Summarize interaction prevalence across bootstrap samples
  require(data.table)
  nbs <- length(prev)
  prev0 <- unlist(lapply(prev, function(z) z$i0))
  prev1 <- unlist(lapply(prev, function(z) z$i1))

  prev <- data.table(int=names(prev1), prev1=prev1, prev0=prev0) %>%
    group_by(int) %>%
    summarize(prev1=mean(prev1), prev0=mean(prev0), n=n()/nbs) %>%
    mutate(diff=(prev1-prev0)) %>%
    arrange(desc(diff))
  return(prev)
}

prevalence <- function(int, nf, wt=rep(1, ncol(nf))) {
  # Calculate the decision path prevalence of an interaction
  id.rm <- wt == 0
  wt <- wt[!id.rm]

  intord <- length(int)
  if (intord == 1)
    int.id <- nf[!id.rm, int] != 0
  else
    int.id <- Matrix::rowSums(nf[!id.rm, int] != 0) == intord

  prev <- sum(wt[int.id]) / sum(wt)
  return(prev)
}

subsetReadForest <- function(rforest, subset.idcs) {
  # Subset nodes from readforest output 
  if (!is.null(rforest$node.feature))
    rforest$node.feature <- rforest$node.feature[subset.idcs,]

  if (!is.null(rforest$tree.info))
    rforest$tree.info <- rforest$tree.info[subset.idcs,]

  if (!is.null(rforest$node.obs))
    rforest$node.obs <- rforest$node.obs[subset.idcs,]

  return(rforest)
}

