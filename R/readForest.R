readForest <- function(rfobj  # a randomForest object with forest component in it
                       , X   # n x p data matrix 
                       , return_node_feature=TRUE
                       , return_node_data=TRUE
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
  out$tree_info <- lapply(1:rfobj$ntree, function(k) getTree(rfobj, k))
  parent <- lapply(out$tree_info, getParent)
  n.node.t <- rfobj$forest$ndbigtree
  out$tree_info <- data.frame(do.call(rbind, out$tree_info), 
                              row.names=NULL)
  out$tree_info$tree <- rep(1:rfobj$ntree, times=n.node.t)
  out$tree_info$parent <- unlist(parent)
  parents <- out$tree_info$parent

  # Repeat each leaf node in node_feature based on specified sampling:
  # importance sampling = size_node
  # uniform = 1
  rep.node <- rep(0, nrow(out$tree_info))
  select.node <- out$tree_info$status == -1

  tt <- matrix(0L, nrow=n, ncol=sum(select.node))
  obs.nodes <- apply(rfobj$obs.nodes, MAR=2, matchOrder) - 1
  leaf.obs <- nodeObs(obs.nodes, n, ntree, table(out$tree_info$tree[leaf.node]), tt)
  out$tree_info$size_node <- 0
  out$tree_info$size_node[leaf.node] <- colSums(leaf.obs)
 
  select.node <- select.node & subsetFun(out$tree_info)
  out$tree_info <- out$tree_info[select.node,]
  rep.node[select.node] <- trunc(wtFun(out$tree_info))

  if (return_node_feature) {
    var.nodes <- rfobj$forest$bestvar
    total.rows <- sum(rep.node[select.node])
    node.vars <- matrix(0L, nrow=(total.rows * p), ncol=2)
    sparse.idcs <- nodeVars(var.nodes, as.integer(ntree), 
                            as.integer(nrow(var.nodes)), 
                            as.integer(p),
                            as.integer(parents),
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

getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent <- parent %% nrow(tree.info)
  parent[1] <- 0
  return(parent)
}
