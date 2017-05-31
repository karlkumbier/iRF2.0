readForest <- function(rfobj, x, y=NULL, 
                       return.node.feature=TRUE,
                       subsetFun = function(x) rep(TRUE, nrow(x)),
                       wtFun = function(x) x$size.node,
                       n.core=1){
  
  require(data.table)
  require(parallel)
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
 
  ntree <- rfobj$ntree
  p <- ncol(x)
  n <- nrow(x)
  out <- list()
  
  rd.forest <- mclapply(1:ntree, readTree, rfobj=rfobj, x=x, y=y,
                        return.node.feature=return.node.feature,
                        subsetFun=subsetFun, wtFun=wtFun,
                        mc.cores=n.core)
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # aggregate sparse feature matrix across forest
  temp <- lapply(rd.forest, function(tt) tt$node.feature)
  row.offset <- c(0, cumsum(sapply(temp, function(z) max(z[,1])))[-ntree])
  n.rows <- sapply(temp, nrow)
  temp <- do.call(rbind, temp)
  temp[,1] <- temp[,1] + rep(row.offset, times=n.rows)
  out$node.feature <- sparseMatrix(i=temp[,1], j=temp[,2], 
                                   dims=c(max(temp[,1]), p))
  return(out)
  
}

readTree <- function(rfobj, k, x, y, return.node.feature, subsetFun, wtFun) {
  
  n <- nrow(x)
  p <- ncol(x)
  ntree <- rfobj$ntree
  
  # Read tree level data from RF
  out <- list()
  out$tree.info <- as.data.frame(getTree(rfobj, k))
  parents <- getParent(out$tree.info)
  out$tree.info$parent <- parents
  out$tree.info$tree <- k
  
  # Repeat each leaf node in node.feature based on specified sampling:
  # importance sampling = size_node; uniform = 1
  select.node <- out$tree.info$status == -1
  rep.node <- rep(0, nrow(out$tree.info))
  
  if (is.null(rfobj$obs.nodes)) {
    fit.data <- passData(rfobj, x, out$tree.info, k)
    leaf.counts <- rowSums(fit.data[out$tree.info$status == -1,])
    which.leaf <- apply(fit.data[out$tree.info$status == -1,], MAR=2, which)
    if (!is.null(y)) leaf.sd <- c(by(y, which.leaf, sdNode))
  } else {
    leaf.counts <- unname(table(rfobj$obs.nodes[,k]))
    if (!is.null(y)) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
  }
  out$tree.info$size.node[select.node] <- leaf.counts
  
  select.node <- select.node & subsetFun(out$tree.info)
  out$tree.info <- out$tree.info[select.node,]
  if (!is.null(y)) {
    out$tree.info$purity <- leaf.sd
    out$tree.info$dec.purity <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  rep.node[select.node] <- trunc(wtFun(out$tree.info))
  
  # Extract decision paths from leaf nodes as binary sparse matrix
  if (return.node.feature) {
    row.offset <- 0
    var.nodes <- as.integer(rfobj$forest$bestvar[,k])
    total.rows <- sum(rep.node[select.node])
    sparse.idcs <- nodeVars(var.nodes, 
                            as.integer(length(select.node)), 
                            as.integer(p),
                            as.integer(parents),
                            as.integer(select.node),
                            as.integer(rep.node),
                            as.integer(row.offset),
                            matrix(0L, nrow=(total.rows * p), ncol=2))
    out$node.feature <- sparse.idcs[!sparse.idcs[,1] == 0,]
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


passData <- function(rfobj, x, tt, k) {
  # Pass data through rf object
  
  leaf.id <- tt$status == -1
  n <- nrow(x)
  n.node <- rfobj$forest$ndbigtree[k]
  node.composition <- matrix(FALSE, nrow=n.node, ncol=n)
  node.composition[1,] <- TRUE
  
  
  for (i in which(!leaf.id)){
    
    # determine children and split point for current node
    d.left <- tt$"left daughter"[i]
    d.right <- tt$"right daughter"[i]
    split.var <- tt$"split var"[i]
    split.pt <- tt$"split point"[i]
    
    parent.id <- node.composition[i,]
    d.left.id <- (x[,split.var] < split.pt) & parent.id
    d.right.id <- (x[,split.var] >= split.pt) & parent.id
    
    node.composition[d.left,] <- d.left.id
    node.composition[d.right,] <- d.right.id
  }
  
  return(node.composition)
}

sdNode <- function(x) ifelse(length(x) == 1, 0, sd(x))
