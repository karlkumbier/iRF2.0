#' @importFrom data.table data.table
#' @importFrom Matrix colSums rowSums t
#' @importFrom utils combn
#' @importFrom fastmatch fmatch
intImportance <- function(int, nf, precision, select.id, weight) {
  # Calculate importance metrics for an interaction across selected elaf nodes
  # of a fitted random forest.
  #   prev: the prevalence of an interaction across all selected leaf nodes,
  #   weighted by observations in each leaf node.
  #   prec: the proportion of class-1 observations in leaf nodes containing the
  #   interaction.
  
  # Determine which leaf nodes contain the given interaction
  int.id <- Reduce(intersect, nf[as.character(int)]) + 1
  if (length(int.id) == 0)
    return(data.table(prev1=0, prev0=0, prec=0))

  # Compute prevalence and precision for given interaction
  prev1 <- prevalence(weight, int.id, select.id)
  prev0 <- prevalence(weight, int.id, !select.id)
  prec <- sum(precision[int.id] * weight[int.id]) / sum(weight[int.id])
  return(data.table(prev1=prev1, prev0=prev0, prec=prec))
}

prevalence <- function(weight, idint, idsel) {
  # Computes the prevalence in selected nodes 
  sint <- sum(weight[idint][idsel[idint]])
  s <- sum(weight[idsel])
  return(sint / s)
}

nodePrecision <- function(read.forest, y, count, weights=rep(1, length(y))) {
  # Evaluate class proportion of class-1 observations in each leaf node.
  if (is.factor(y)) y <- as.numeric(y) - 1
  count.y <- Matrix::colSums(read.forest$node.obs * y * weights)
  precision <- count.y / count
  precision[count == 0] <- 0
  return(precision)
}

subsetTest <- function(int, ints, importance) {
  # Compare interaction importance metrics to simple baselines.
  # Prevalence: expectation under independent selection.
  # Precision: precision of lower-order interactions.

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  id.int <- fmatch(list(int), ints)

  # Evaluate prevalence relative to independent selection
  id <- sapply(int, fmatch, ints)
  prev.null <- prod(importance$prev1[id])
  prev <- importance$prev1[id.int]
  prev.test <- prev - prev.null

  # Determine all s-1 order interactions to evaluate
  ss <- combn(int, length(int) - 1, simplify=FALSE)
  id <- sapply(ss, function(ii) fmatch(list(ii), ints))
  prec.null <- importance$prec[id]
  prec <- importance$prec[id.int]
  prec.test <- min(prec - prec.null)

  return(data.table(prev.test=prev.test, prec.test=prec.test))
}
