interactPredict <- function(int, x, rd.forest, min.node=1, qcut=0.5, class=1,
                            hard.region=TRUE, varnames.grp=1:ncol(x)) {
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

  id.cls <- rd.forest$tree.info$prediction == class + 1
  id.min <- rd.forest$tree.info$size.node >= min.node
  nf <- rd.forest$node.feature[id.min & id.cls,]
  tree.info <- rd.forest$tree.info[id.min & id.cls,]
  
  p <- length(unique(varnames.grp))
  stopifnot(ncol(x) == p)
  id <- int2Id(int, varnames.grp, directed=TRUE)
  id.raw <- id - (id > p) * p
  sgn <- intSign(int)
  
  if (length(id) == 1) {
    id.int <- which(nf[,id] != 0)
  } else {
    id.int <- apply(nf[,id], MAR=1, function(z) all(z != 0))
  }
    
  if (sum(id.int) == 0) {
    warning('0 paths contain interaction')
    return(rep(0, nrow(x)))
  }
  
  nf.sub <- t(matrix(nf[id.int, id], nrow=sum(id.int)))
  nf.sub[sgn < 0,] <- nf.sub[sgn < 0,] * -1
  xraw <- as.matrix(x[,id.raw])
  xraw[,sgn < 0] <- xraw[,sgn < 0] * -1
  
  if (hard.region) {
    qq <- function(x) quantile(x, probs=qcut)
    nf.sub <- as.matrix(apply(nf.sub, MAR=1, qq))
  } 
  node.pred <- apply(xraw, MAR=1, regionPredObs, thresh=nf.sub)
  return(node.pred)
}

regionPredObs <- function(x, thresh) {
  # Generates prediction for single observation based on whether it falls in 
  # regions indicated by thresholds.
  # args:
  #   x: numeric vector
  #   thresh: matrix of thresholds, rows corresponding to interacting features
  stopifnot(length(x) == nrow(thresh))
  pred <- mean(colSums(x > thresh) == nrow(thresh))
  return(pred)  
}

int2Id <- function(int, varnames.grp, directed=FALSE, split=FALSE) {
  # Determine integer index of named variable in nf matrix (directed or not)
  if (!split) int <- strsplit(int, '_')[[1]]
  if (directed) {
    dir <- grep('\\+$', int)  
    varnames.grp <- gsub('[\\+\\-]', '', varnames.grp)
    int <- gsub('[\\+\\-]', '', int)
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
