intImportance <- function(int, nf, prec.nd, select.id, weight) {
  # Calculate the prevalence of an interaction across selected leaf nodes of a
  # random forest
  
  # Remove all 0-weighted leaf nodes from further analysis
  id.rm <- weight == 0
  select.id <- select.id[!id.rm]
  weight <- weight[!id.rm]
  nf <- nf[!id.rm,]
  prec.nd <- prec.nd[!id.rm]

  intord <- length(int)
  if (intord == 1)
    int.id <- nf[, int] != 0
  else
    int.id <- Matrix::rowSums(nf[, int] != 0) == intord

  if (sum(int.id) == 0)
    return(data.table(prev1=0, prev0=0, prec=0))

  # Compute prevalence and precision for given interaction
  prev <- prevalence(weight, int.id, select.id)
  prec <- mean(prec.nd[int.id & select.id])
  return(data.table(prev1=prev[1], prev0=prev[2], prec=prec))
}

prevalence <- function(weight, idint, idcl) {
  # Computes the prevalence of an interaction among class-1 and class-0 leaf 
  # nodes in the fitted RF.
  sint1 <- sum(weight[idint & idcl])
  sint0 <- sum(weight[idint & !idcl])
  s1 <- sum(weight[idcl])
  s0 <- sum(weight[!idcl])

  prev1 <- sint1 / s1
  prev0 <- sint0 / s0
  return(c(prev1, prev0))
}

nodeAttr <- function(read.forest, y, weight=rep(1, length(y))) {
  # Evaluate class proportion of class-1 observations in each leaf node of  a
  # fitted RF.
  if (is.factor(y)) y <- as.numeric(y) - 1

  ndcnt <- t(read.forest$node.obs)
  ndcntY <- Matrix::colSums(ndcnt * y * weight)
  ndcnt <- Matrix::colSums(ndcnt * weight)
  prec.nd <- ndcntY / ndcnt
  prec.nd[ndcnt == 0] <- 0
  return(list(precision=prec.nd, ndcnt=ndcnt))
}

subsetTest <- function(int, ints, importance) {
  # Compare interaction importance metrics to simple baselines.
  # Prevalence: expectation under independent selection.
  # Precision: precision of lower-order interactions.

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  
  # Determine all s-1 order interactions to evaluate
  ss <- combn(int, length(int) - 1, simplify=FALSE)
  setEq <- function(x, y) all(x %in% y) & all(y %in% x)
  setsEq <- function(x, y) sapply(y, setEq, x=x)
  getPair <- function(z) setsEq(z, ints) | setsEq(setdiff(int, z), ints)
  pairs <- lapply(ss, getPair)

  id <- setsEq(int, ints)
  prev <- importance$prev1[id]
  prev.null <- sapply(pairs, function(z) prod(importance$prev1[z]))
  prec <- importance$prec[id]
  prec.null <- importance$prec[unlist(lapply(pairs, which))]
  return(data.table(prev.test=min(prev - prev.null),
                    prec.test=min(prec - prec.null)))
}
