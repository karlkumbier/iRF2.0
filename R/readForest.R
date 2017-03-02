#   PASS DATA THROUGH FOREST
readForest <- function(rfobj  # a randomForest object with forest component in it
                       , X   # n x p data matrix 
                       , return_node_feature=TRUE
                       , return_node_data=TRUE
                       , leaf_node_only = TRUE
                       , subsetFun = function(x) rep(TRUE, nrow(x))
                       , wtFun = function(x) x$size_node
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
  if (ntree > 1) {
    out$tree_info <- lapply(1:rfobj$ntree, function(k) getTree(rfobj, k))
    n.node.t <- sapply(out$tree_info, nrow)
    out$tree_info <- as.data.frame(do.call(rbind, out$tree_info), stringsAsFactors=FALSE)
    out$tree_info$tree <- rep(1:rfobj$ntree, times=n.node.t)
  } else {
    out$tree_info <- as.data.frame(getTree(rfobj, 1), stringsAsFactors=FALSE)
    n.node.t <- nrow(out$tree_info)
    out$tree_info$tree <- 1
  }


  # How many times should each node be repeated in node_feature matrix:
  # size_node = importance sampling, 1 = uniform
  rep.node <- rep(0, nrow(out$tree_info))
  select.node <- rep(TRUE, nrow(out$tree_info))
  if (leaf_node_only) select.node <- out$tree_info$status == -1 
  
  tt <- matrix(0L, nrow=n, ncol=nrow(out$tree_info))
  obs.nodes <- rfobj$obs.node
  obs.nodes <- apply(obs.nodes, MAR=2, matchOrder) - 1
  leaf.obs <- nodeObs(obs.nodes, n, ntree, table(out$tree_info$tree[select.node]), tt)
  out$tree_info$size_node <- colSums(leaf.obs)
  
  select.node <- select.node & subsetFun(out$tree_info)
  out$tree_info <- out$tree_info[select.node,]
  rep.node[select.node] <- trunc(wtFun(out$tree_info))

  if (return_node_feature) {
    
    var.nodes <- rfobj$feature.node
    total.rows <- sum(rep.node[select.node])
    node.vars <- matrix(0L, nrow=(total.rows * p), ncol=2)
    sparse.idcs <- nodeVars(var.nodes, p, ntree, nrow(var.nodes),
                                 as.integer(select.node),
                                 as.integer(rep.node), 
                                 as.integer(n.node.t),
                                 node.vars)
    sparse.idcs <- sparse.idcs[!sparse.idcs[,1] == 0,]
    out$node_feature <- sparseMatrix(i=sparse.idcs[,1], j=sparse.idcs[,2], dims=c(total.rows, p))

  }

  if (return_node_data) out$node_data <- leaf.obs 

  return(out)
}

matchOrder <- function(x) {
  x.sorted <- sort(unique(x))
  x.matched <- match(x, x.sorted)
  return(x.matched)
}
