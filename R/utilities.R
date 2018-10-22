pasteInt <- function(x) {
  # Combine interaction into single string
  x <- paste(x, collapse='_')
  return(x)
}

nameInts <- function(ints, varnames, signed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by
  # variable names. Naming convention for an interaction is:
  #   <variable1(sign)>_ <variable2(sign)>_...

  varnames <- unique(varnames)
  p <- length(varnames)
  if (signed)
    signs <- lapply(ints, function(z) ifelse(z > p, '+', '-'))
  else
    signs <- ''

  # Adjust indexing to match varnames
  ints <- lapply(ints, function(z) z %% p + p * (z == p | z == 2 * p))
  ints.name <- mapply(function(i, s) nameInt(varnames, i, s), ints, signs)
  return(ints.name)
}


nameInt <- function(varnames, idx, sgn) {
  int <- paste0(varnames[idx], sgn)
  int <- paste(sort(int), collapse='_')
  return(int)
}
