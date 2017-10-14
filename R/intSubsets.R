intsSubsets <- function(ints) return(unique(unlist(lapply(ints, intSubsets))))

intSubsets <- function(int) {
  # Generates all lower-order subsets of a given interaction
  vv <- strsplit(int, '_')[[1]]
  int.orders <- 1:length(vv)
  ints <- lapply(int.orders, combn, x=vv, simplify=FALSE)
  if (length(int.orders) > 1) ints <- unlist(ints, recursive=FALSE)
  ints <- sapply(ints, paste, collapse='_')
  return(ints)
}
