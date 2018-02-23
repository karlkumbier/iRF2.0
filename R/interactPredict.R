interactPredict <- function(int, x, rd.forest, min.node=1, 
                            max.rule=1000, class=1,
                            varnames.grp=1:ncol(x)) {
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
  nf <- rd.forest$node.feature[id.cls & id.min,]
  tree.info <- rd.forest$tree.info[id.cls & id.min,]
  
  p <- length(unique(varnames.grp))
  stopifnot(ncol(x) == p)
  id <- int2Id(int, varnames.grp, directed=TRUE)
  id.raw <- id - (id > p) * p

  if (length(id.raw) == 1) {
    id.int <- which(nf[,id] != 0)
  } else {
    id.int <- apply(nf[,id], MAR=1, function(z) all(z != 0))
  }

  if (sum(id.int) == 0) {
    warning('0 paths contain interaction')
    return(rep(0, nrow(x)))
  }

  id.sel <- sample(which(id.int), min(max.rule, sum(id.int)))
  nf.sub <- t(nf[id.sel , c(id.raw, id.raw + p)])
  xraw <- as.matrix(x[,id.raw])
  nf.array <- array(nf.sub, dim=c(nrow(nf.sub)/2, 2, ncol(nf.sub)))
  print(dim(nf.array))
  # Bound regions by min/max values of each feature
  xmin <- apply(xraw, MAR=2, min)
  xmax <- apply(xraw, MAR=2, max)
  
  nf.array[,2,] <- apply(nf.array[,2,], MAR=2, function(z) {
    z[z == 0] <- xmin[z == 0]
    return(z)
  })
  
  nf.array[,1,] <- apply(nf.array[,1,], MAR=2, function(z) {
    z[z == 0] <- xmax[z == 0]
    return(z)
  })

  node.pred <- apply(xraw, MAR=1, regionPredObs, thresh=nf.array)
  return(node.pred)
}

regionPredObs <- function(x, thresh) {
  # Generates prediction for single observation based on whether it falls in 
  # regions indicated by thresholds.
  # args:
  #   x: numeric vector
  #   thresh: matrix of thresholds, rows corresponding to interacting features
  
  pred <- mean(x >= thresh[,2,] & x <= thresh[,1,])
  #pred <- mean(colSums(x >= thresh[,2,] & x <= thresh[,1,]) == nrow(thresh))
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
