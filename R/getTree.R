#' getTree
#'
#' Read out tree metadata from randomForest or ranger objects
#'
#' @param rfobj an object of class randomForest or ranger
#' @param k tree index to read
#' @param labelVar
#' @param terminal node membership, required for reading ranger
#'
#' @import dplyr
#' @importFrom fastmatch "%fin%"
getTree <- function(x, ...) UseMethod("getTree")

getTree.default <- function(...)
    stop(deparse(substitute(rfobj)), "is not of class ranger or randomForest")

getTree.ranger <- function(rfobj, k=1) {
  # Read metadata from forest
  tree.info <- ranger::treeInfo(rfobj, k) %>%
      transmute(`left daughter` = leftChild+1L,
                `right daughter` = rightChild+1L,
                `split var` = splitvarID+1L,
                `split point` = splitval,
                status = terminal,
                prediction = prediction)

  return(tree.info)
}

getTree.randomForest <- function(rfobj, k=1) {
  # Read metadata from forest
  tree.info <- randomForest::getTree(rfobj, k) %>%
      as.data.frame %>%
      mutate(status = status==-1L)

  return(tree.info)
}

