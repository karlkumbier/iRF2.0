#   PASS DATA THROUGH FOREST
readForest <- function(rfobj  # a randomForest object with forest component in it
                       , X   # n x p data matrix 
                       , return_node_feature=TRUE
                       , return_node_data=TRUE
                       , leaf_node_only = TRUE
                       ){

  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (!leaf_node_only)
    stop('readForest not currently implemented for non-leaf nodes')

  n <- nrow(X)
  p <- ncol(X)
  nrnodes <- rfobj$nrnodes
  ntree <- rfobj$ntree


  out <- list()
  K <- rfobj$ntree
  if (K > 1) {
    out$tree_info <- lapply(1:rfobj$ntree, function(k) getTree(rfobj, k))
    n.node.t <- sapply(out$tree_info, nrow)
    out$tree_info <- as.data.frame(do.call(rbind, out$tree_info), stringsAsFactors=FALSE)
    out$tree_info$tree <- rep(1:rfobj$ntree, times=n.node.t)
  } else {
    out$tree_info <- as.data.frame(getTree(rfobj, 1), stringsAsFactors=FALSE)
    n.node.t <- nrow(out$tree_info)
    out$tree_info$tree <- 1
  }

  select.node <- 1:nrow(out$tree_info)
  if (leaf_node_only) {
    select.node <- which(out$tree_info[,'status'] == -1)
    out$tree_info <- out$tree_info[select.node,]
    n.node.t.sub <- table(out$tree_info$tree)
  }

  if (return_node_feature) {
    tt <- matrix(0L, nrow=p, ncol=nrow(out$tree_info))
    var.nodes <- rfobj$feature.node
    nr <- nrow(var.nodes)
    out$node_feature <- nodeVars(var.nodes, p, ntree, nrnodes, 
                                 nr, n.node.t, tt, select.node)
  }

  tt <- matrix(0L, nrow=n, ncol=nrow(out$tree_info))
  obs.nodes <- rfobj$obs.node
  obs.nodes <- apply(obs.nodes, MAR=2, matchOrder) - 1
  leaf.obs <- nodeObs(obs.nodes, n, ntree, n.node.t.sub, tt)
  out$tree_info$size_node <- colSums(leaf.obs)

  if (return_node_data) out$node_data <- leaf.obs 

  return(out)
}

matchOrder <- function(x) {
  x.sorted <- sort(unique(x))
  x.matched <- match(x, x.sorted)
  return(x.matched)
}
