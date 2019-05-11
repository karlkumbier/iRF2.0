#' getTree
#'
#' Read out tree metadata from randomForest or ranger objects
#'
#' @param rfobj an object of class randomForest or ranger
#' @param k tree index to read
#' @param labelVar
#' @param terminal node membership, required for reading ranger
#'
#' @importFrom fastmatch "%fin%"
getTree <- function(rfobj, k=1, nodes=NULL) {
   
  if (class(rfobj) == 'randomForest') {
    out <- getTreeRF(rfobj, k)
  } else if (class(rfobj) == 'ranger') {
    out <- getTreeRanger(rfobj, k, nodes)
  } else {
    stop(deparse(substitute(rfobj)), "is not class ranger of randomForest")
  }

  return(out)
}

getTreeRanger <- function(rfobj, k=1, nodes=NULL) {
  # Check whether current tree can be read
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (is.null(nodes)) {
    stop("Terminal node membership missing")
  }
  if (k > rfobj$num.trees) {
    stop("There are fewer than ", k, "trees in the forest")
  }

  # Read metadata from forest
  nnode <- length(rfobj$forest$split.values[[k]])
  status <- 1:nnode %fin% (nodes[,k] + 1)
  predicted <- rep(0L, nnode)
  predicted[status] <- rfobj$forest$split.values[[k]][status]
  tree.info <- data.frame(rfobj$forest$child.nodeIDs[[k]][[1]] + 1,
                          rfobj$forest$child.nodeIDs[[k]][[2]] + 1,
                          rfobj$forest$split.varIDs[[k]],
                          rfobj$forest$split.values[[k]],
                          ifelse(status, -1, 1),
                          predicted)

  colnames(tree.info) <- c("left daughter", "right daughter", "split var",
                           "split point", "status", "prediction")

  return(tree.info)

}

getTreeRF <- function(rfobj, k=1) {
  # Check whether current tree can be read
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }

  # Read metadata from forest
  if (rfobj$type == "regression") {
      tree.info <- data.frame(rfobj$forest$leftDaughter[,k],
                              rfobj$forest$rightDaughter[,k],
                              rfobj$forest$bestvar[,k],
                              rfobj$forest$xbestsplit[,k],
                              rfobj$forest$nodestatus[,k],
                              rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  } else {
      tree.info <- data.frame(rfobj$forest$treemap[,,k],
                              rfobj$forest$bestvar[,k],
                              rfobj$forest$xbestsplit[,k],
                              rfobj$forest$nodestatus[,k],
                              rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  }

  colnames(tree.info) <- c("left daughter", "right daughter", "split var", 
                           "split point", "status", "prediction")

  return(tree.info)
}

