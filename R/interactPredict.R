interactPredict <- function(x, int, read.forest, varnames.grp=1:ncol(x), 
                            hard.region=FALSE, qcut=0.5, nrule=1000, 
                            min.node=1, is.split=FALSE) {
  # Generate RF predictions for given interactions, using only information
  # from interacting features.
  p <- ncol(x)
  stopifnot(p == length(varnames.grp))
  stopifnot(p == ncol(read.forest$node.feature) / 2)
  
  tree.info <- read.forest$tree.info[tree.info$size.node >= min.node,]
  nf <- read.forest$node.feature[tree.info$size.node >= min.node,]
  
  # Get feature indices for interaction terms 
  if (!is.split) int <- strsplit(int, '_')[[1]]
  int.adj <- isPositive(int)
  int.unsgn <- intUnsign(int) 
  id <- mapply(function(i, a) {
    intId(int=i, varnames=varnames.grp, adj=a)
  }, int.unsgn, int.adj, SIMPLIFY=TRUE)
  id.raw <- id %% p  + p * (id == p | id == 2 * p) 
  id.pos <- id > p
  
  # if classification, subset to class 1 leaf nodes
  if (all(tree.info$prediction %in% 1:2)) {
    nf <- nf[tree.info$prediction == 2,]
    tree.info <- tree.info[tree.info$prediction == 2,]
    tree.info$prediction <- tree.info$prediction - 1
  }
  
  # Subset node feature matrix to and data matrix to interacting features
  nf <- nf[,id]
  x <- x[,id.raw]
  if (is.null(dim(nf))) {
    int.nds <- nf != 0
    nf <- as.matrix(nf[int.nds])
    x <- as.matrix(x)
  } else {
    int.nds <- Matrix::rowSums(nf != 0) > 0
    nf <- nf[int.nds,]
  }
  tree.info <- tree.info[int.nds,]
  
  # Set response values for each region proportional to node size
  if (hard.region) {
    nf <- matrix(apply(nf, MAR=2, quantile, probs=qcut), nrow=1)
    y <- 1
    size <- 1
  } else {
    y <- tree.info$prediction
    size <- tree.info$size.node
  }
  
  # evaluate predictions over subsample of active rules
  nrule <- min(nrule, nrow(nf))
  ss <- sample(nrow(nf), nrule, prob=size)
  preds <- matrix(0, nrow=nrow(x), ncol=nrule)
  for (i in 1:nrule) {
    s <- ss[i]
    id.active <- nf[s, ] != 0
    tlow <- t(x[,!id.pos & id.active]) <= nf[s, !id.pos & id.active]
    if (is.null(ncol(tlow))) tlow <- matrix(0, nrow=0, ncol=nrow(x))
    
    thigh <-  t(x[,id.pos & id.active]) > nf[s, id.pos & id.active]
    if (is.null(ncol(thigh))) thigh <- matrix(0, nrow=0, ncol=nrow(x))
    
    int.active <- (colSums(tlow) + colSums(thigh)) == sum(id.active)
    preds[,i] <- int.active * size[s] * y[s]
  }
  return(rowSums(preds) / sum(size[ss]))
}

interactPredictPermute <- function(x, int, read.forest, varnames.grp=1:ncol(x), 
                                   hard.region=FALSE, qcut=0.5, nrule=1000, 
                                   min.node=1, is.split=FALSE, nperm=100
                                   ) {

  if (!is.split) int <- strsplit(int, '_')[[1]]
  int <- intUnsign(int)
  id.perm <- which(varnames.grp %in% int)
  out <- list()
  for (i in 1:length(id.perm)) {
    permPred  <- function(ii) {
      xperm <- x
      xperm[,ii] <- sample(xperm[,ii])
      out <- interactPredict(int, x=xperm, varnames.grp=varnames.grp, 
                             read.forest=read.forest, hard.region=hard.region, 
                             qcut=qcut, nrule=nrule, min.node=min.node, 
                             is.split=is.split)
    }
    
    out[[i]] <- replicate(nperm, permPred(id.perm[[i]]))
  }
  return(out)
}

isPositive <- function(x) {
  out <- rep(FALSE, length(x))
  out[grep('+', x, fixed=TRUE)] <- TRUE
  return(out)
}

intUnsign <- function(x) gsub('[-\\+]', '', x)

intId <- function(int, varnames.grp, adj) { 
  # Evaluate 1:2p index of interaction term
  int.adj <- which(varnames.grp == int) + adj * length(varnames.grp)
  return(int.adj)
}
