#' @importFrom data.table data.table
#' @importFrom Matrix colSums rowSums t
#' @importFrom utils combn
intImportance <- function(int, nf, precision, select.id, weight) {
  # Calculate importance metrics for an interaction across selected elaf nodes
  # of a fitted random forest.
  #   prev: the prevalence of an interaction across all selected leaf nodes,
  #   weighted by observations in each leaf node.
  #   prec: the proportion of class-1 observations in leaf nodes containing the
  #   interaction.
  
  # Remove all 0-weighted leaf nodes from further analysis
  id.rm <- weight == 0
  select.id <- select.id[!id.rm]
  weight <- weight[!id.rm]
  nf <- nf[!id.rm,]
  precision <- precision[!id.rm]

  # Determine which leaf nodes contain the given interaction
  if (length(int) == 1)
    int.id <- nf[, int] != 0
  else
    int.id <- Matrix::rowMeans(nf[, int] != 0) == 1

  if (sum(int.id) == 0)
    return(data.table(prev1=0, prev0=0, prec=0))

  # Compute prevalence and precision for given interaction
  prev1 <- prevalence(weight, int.id, select.id)
  prev0 <- prevalence(weight, int.id, !select.id)
  prec <- mean(precision[int.id])
  return(data.table(prev1=prev1, prev0=prev0, prec=prec))
}

prevalence <- function(weight, idint, idsel) {
  # Computes the prevalence in selected nodes 
  sint <- sum(weight[idint & idsel])
  s <- sum(weight[idsel])
  return(sint / s)
}

nodeCount <- function(read.forest, weight=rep(1, nrow(read.forest$node.obs))) {
  # Evaluate the number of observations in each leaf node
  count <- Matrix::colSums(read.forest$node.obs * weight)
  return(count)
}

nodePrecision <- function(read.forest, y, count, weight=rep(1, length(y))) {
  # Evaluate class proportion of class-1 observations in each leaf node.
  if (is.factor(y)) y <- as.numeric(y) - 1
  count.y <- Matrix::colSums(read.forest$node.obs * y * weight)
  precision <- count.y / count
  precision[count == 0] <- 0
  return(precision)
}

subsetTest <- function(int, ints, importance) {
  # Compare interaction importance metrics to simple baselines.
  # Prevalence: expectation under independent selection.
  # Precision: precision of lower-order interactions.

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  id.int <- match(list(int), ints)

  # Evaluate prevalence relative to independent selection
  id <- sapply(int, match, ints)
  prev.null <- prod(importance$prev1[id])
  prev <- importance$prev1[id.int]
  prev.test <- prev - prev.null

  # Determine all s-1 order interactions to evaluate
  ss <- combn(int, length(int) - 1, simplify=FALSE)
  id <- sapply(ss, function(ii) match(list(ii), ints))
  prec.null <- importance$prec[id]
  prec <- importance$prec[id.int]
  prec.test <- mean(prec - prec.null)

  return(data.table(prev.test=prev.test, prec.test=prec.test))
}
