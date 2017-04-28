#   PASS DATA THROUGH FOREST
readForest <- function(rfobj  # a randomForest object with forest object
                       , x   # n x p data matrix 
                       , return.node.feature=TRUE
                       , subsetFun = function(x) rep(TRUE, nrow(x))
                       , wtFun = function(x) x$size_node
                       ){
  
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')

  n <- nrow(x)
  p <- ncol(x)
  nrnodes <- rfobj$nrnodes
  ntree <- rfobj$ntree


  # Read tree level data from RF
  out <- list()
  out$tree.info <- lapply(1:rfobj$ntree, function(k) getTree(rfobj, k))
  parent <- lapply(out$tree.info, getParent)
  out$tree.info <- as.data.frame(do.call(rbind, out$tree.info), 
                                 stringsAsFactors=FALSE)
  
  n.node.t <- rfobj$forest$ndbigtree
  out$tree.info$tree <- rep(1:rfobj$ntree, times=n.node.t)
  out$tree.info$parent <- unlist(parent)
  parents <- out$tree.info$parent

  # Repeat each leaf node in node.feature based on specified sampling:
  # importance sampling = size_node
  # uniform = 1
  leaf.node <- out$tree.info$status == -1
  rep.node <- rep(0, nrow(out$tree.info))
  select.node <- leaf.node
  

  leaf.counts <- unname(unlist(apply(rfobj$obs.node, MAR=2, table)))
  out$tree.info$size_node[leaf.node] <- leaf.counts
 
  select.node <- select.node & subsetFun(out$tree.info)
  out$tree.info <- out$tree.info[select.node,]
  rep.node[select.node] <- trunc(wtFun(out$tree.info))

  # Extract decision paths from leaf nodes as binary sparse matrix
  if (return_node.feature) {
    var.nodes <- rfobj$forest$bestvar
    total.rows <- sum(rep.node[select.node])
    sparse.idcs <- nodeVars(var.nodes, as.integer(ntree), 
                            as.integer(nrow(var.nodes)), 
                            as.integer(p),
                            as.integer(parents),
                            as.integer(select.node),
                            as.integer(rep.node), 
                            as.integer(n.node.t),
                            matrix(0L, nrow=(total.rows * p), ncol=2))
    sparse.idcs <- sparse.idcs[!sparse.idcs[,1] == 0,]
    out$node.feature <- sparseMatrix(i=sparse.idcs[,1], j=sparse.idcs[,2], 
                                     dims=c(total.rows, p))

  }
  return(out)
}

getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent <- parent %% nrow(tree.info)
  parent[1] <- 0
  return(parent)
}
