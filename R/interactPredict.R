interactPredict <- function(x, int, read.forest, varnames.grp=1:ncol(x), 
                            hard.region=FALSE, qcut=0.5, nrule=1000, 
                            min.node=1, is.split=FALSE) {
  # Generate RF predictions for given interactions, using only information
  # from interacting features.
  
  p <- ncol(x)
  stopifnot(p == length(varnames.grp))
  stopifnot(p == ncol(read.forest$node.feature) / 2)
  
  tree.info <- read.forest$tree.info
  nf <- read.forest$node.feature
  
  # Get interaction term indices
  if (!is.split) int <- strsplit(int, '_')[[1]]
  
  int.adj <- isPositive(int)
  int.unsgn <- intUnsign(int) 
  id <- mapply(function(i, a) {
    int2Id(int=i, varnames=varnames.grp, adj=a)
  }, int.unsgn, int.adj, SIMPLIFY=TRUE)
  
  nf <- nf[tree.info$size.node >= min.node,]
  tree.info <- tree.info[tree.info$size.node >= min.node,]
  
  # if classification, subset to class 1 leaf nodes
  if (all(tree.info$prediction %in% 1:2)) {
    nf <- nf[tree.info$prediction == 2,]
    tree.info <- tree.info[tree.info$prediction == 2,]
    tree.info$prediction <- tree.info$prediction - 1
  }
  
  # subset node feature matrix to interacting features
  int.nds <- Matrix::rowSums(nf[,id] != 0) > 0
  nf <- nf[int.nds, id]
  tree.info <- tree.info[int.nds,]
  id.pos <- id > p
  id.raw <- id %% p   
  id.raw <- id.raw + p * (id.raw == 0)
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
    tlow <- t(x[,id.raw[!id.pos & id.active]]) <= nf[s, !id.pos & id.active]
    if (is.null(ncol(tlow))) tlow <- matrix(0, nrow=0, ncol=nrow(x))
    
    thigh <-  t(x[,id.raw[id.pos & id.active]]) > nf[s, id.pos & id.active]
    if (is.null(ncol(thigh))) thigh <- matrix(0, nrow=0, ncol=nrow(x))
    
    int.active <- (colSums(tlow) + colSums(thigh)) == sum(id.active)
    preds[,i] <- int.active * size[s] * y[s]
  }
  return(rowSums(preds) / sum(size[ss]))
}

isPositive <- function(x) {
  out <- rep(FALSE, length(x))
  out[grep('+', x, fixed=TRUE)] <- TRUE
  return(out)
}

intUnsign <- function(x) gsub('[-\\+]', '', x)

int2Id <- function(int, varnames.grp, adj) { 
  # Evaluate 1:2p index of interaction term
  int.adj <- which(varnames.grp == int) + adj * length(varnames.grp)
  return(int.adj)
}
