interactPredict <- function(int, x, rd.forest, min.node=1,  
                            qcut=0.5, max.rule=1000, class=1,
                            hard.region=FALSE, varnames.grp=1:ncol(x)) {
  # Generate prediction from directed interaction based on regions 
  # corresponding to class-C leaf nodes.
  # args:
  #   int: directed interaction, features separated by '_'
  #   x: data matrix
  #   rd.forest: readForest output
  #   min.node: minimum leaf node size to use for prediction
  #   qcut: quantile to cut at if doing hard prediction
  #   class: class of interest (indicates what leaf nodes will be used)
  #   hard.region: T/F indicating whether regions should be aggregated before 
  #     prediction -- will give 0/1 prediction
  #   varnames.grp: grouped variable names
  tree.info <- rd.forest$tree.info
  nf <- rd.forest$node.feature

  if (all(tree.info$prediction %in% 1:2)) {
    id.cls <- tree.info$prediction == class + 1
  } else {
    id.cls <- tree.info$prediction > class
  }
  id.min <- tree.info$size.node >= min.node
  nf <- nf[id.cls & id.min,]
  tree.info <- tree.info[id.cls & id.min,]
  
  p <- length(unique(varnames.grp))
  stopifnot(ncol(x) == p)
  id <- int2Id(int, varnames.grp, directed=TRUE)
  id.raw <- id - (id > p) * p
  sgn <- intSign(int)

  if (length(id.raw) == 1) {
    id.int <- nf[,id] != 0
  } else {
    id.int <- apply(nf[,id], MAR=1, function(z) all(z != 0))
  }

  if (sum(id.int) == 0) {
    return(0, nrow(x))
  }

  if (hard.region) {
    qq <- function(x) quantile(x, probs=qcut)
    nf.sub <- matrix(nf[id.int, id], nrow=sum(id.int))
    nf.sub <- apply(nf.sub, MAR=2, qq)
    nf.sub <- matrix(nf.sub, ncol=1)
  } else {
    id.sel <- sample(which(id.int), min(max.rule, sum(id.int)))
    nf.sub <- t(matrix(nf[id.sel , id], nrow=length(id.sel)))
  }
  
  nf.sub[sgn == -1,] <- nf.sub[sgn == -1,] * -1
  xraw <- as.matrix(x[,id.raw])
  xraw[,sgn == -1] <- xraw[,sgn == -1] * -1

  node.pred <- apply(xraw, MAR=1, regionPredObs, thresh=nf.sub)
  return(node.pred)
}

regionPredObs <- function(x, thresh) {
  # Generates prediction for single observation based on whether it falls in 
  # regions indicated by thresholds.
  # args:
  #   x: numeric vector
  #   thresh: matrix of thresholds, rows corresponding to interacting features
  
  pred <- mean(colSums(x >= thresh) == length(x))
  return(pred)  
}

int2Id <- function(int, varnames.grp, directed=FALSE, split=FALSE) {
  # Determine integer index of named variable in nf matrix (directed or not)
  if (!split) int <- strsplit(int, '_')[[1]]
  if (directed) {
    dir <- grep('\\+$', int)  
    varnames.grp <- gsub('[-\\+]', '', varnames.grp)
    int <- gsub('[-\\+]', '', int)
  }
  
  varnames.grp <- unique(varnames.grp)
  id <- sapply(int, function(i) which(varnames.grp == i))
  if (directed) {
    adjust <- rep(0, length(int))
    adjust[dir] <- length(varnames.grp)
    id <- id + adjust
  }
  
  return(id)
}

intSign <- function(int) {
  # Evaluates the direction of each feature in an interaction
  int <- unlist(strsplit(int, '_'))
  dir <- grep('\\+$', int)
  out <- rep(-1, length(int))
  out[dir] <- 1
  return(out)
}
