readForest <- function(rfobj, x, y=NULL,
                       varnames.grp,
                       return.node.feature=TRUE,
                       wt.pred.accuracy=FALSE,
                       n.core=1){
  
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (wt.pred.accuracy & is.null(y))
    stop('y required to evaluate prediction accuracy')
  
  ntree <- rfobj$ntree
  p <- ncol(x)
  n <- nrow(x)

  # read leaf node data from each tree in the forest 
  rd.forest <- mclapply(1:ntree, function(tt) {
    readTree(rfobj=rfobj, k=tt, x=x, y=y,
             varnames.grp=varnames.grp,
             return.node.feature=return.node.feature,
             wt.pred.accuracy=wt.pred.accuracy)
  }, mc.cores=n.core)
  
  out <- list()
  out$tree.info <- rbindlist(lapply(rd.forest, function(z) z$tree.info))
  
  # aggregate sparse feature matrix across forest
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  nf <- aggregateNodeFeature(nf)
  splits <- unlist(lapply(rd.forest, function(z) z$node.splits))
  out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], x=splits, dims=c(max(nf[,1]), 2 * p))
  return(out)
  
}

readTree <- function(rfobj, k, x, y, 
                     varnames.grp,
                     return.node.feature=TRUE, 
                     wt.pred.accuracy=FALSE) {
  
  n <- nrow(x)
  p <- ncol(x)
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  # Get node attributes across RF
  ntree <- rfobj$ntree
  tree.info <- as.data.frame(getTree(rfobj, k))
  tree.info$node.idx <- 1:nrow(tree.info)
  parents <- getParent(tree.info)
  tree.info$parent <- as.integer(parents)
  parent.splits <- parentSplit(tree.info, tree.info$parent)
  tree.info$parent.var <- parent.splits[,1]
  tree.info$parent.split <- parent.splits[,2]
  tree.info$direction <- ifelse(tree.info$parent > nrow(tree.info), 1, 0)
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  leaf.idcs <- tree.info$node.idx[tree.info$status == -1]
  ancestors <- lapply(leaf.idcs, getAncestorPath, tree.info=tree.info)
  
  # HACK TO REMOVE SPLITS OF SINGLE FEATURE
  #id.rm <- sapply(ancestors, function(z) nrow(z) < 2)
  #print(sum(id.rm))
  #ancestors <- ancestors[!id.rm]
  
  node.feature <- sapply(ancestors, path2Binary, p=p, varnames.grp=varnames.grp)
  select.node <- tree.info$status == -1
  rep.node <- rep(0L, nrow(tree.info))
  tree.info <- select(tree.info, prediction, node.idx, parent, tree, size.node)
  
  if (is.null(rfobj$obs.nodes)) {
    # if nodes not tracked, pass data through forest to get leaf counts.
    fit.data <- passData(rfobj, x, tree.info, k)
    leaf.counts <- rowSums(fit.data[select.node,])
    which.leaf <- apply(fit.data[select.node,], MAR=2, which)
    leaf.idx <- as.integer(which(select.node))
    if (wt.pred.accuracy) leaf.sd <- c(by(y, which.leaf, sdNode))
  } else {
    leaf.counts <- table(rfobj$obs.nodes[,k])
    leaf.idx <- as.integer(names(leaf.counts))
    if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
  }
  
  tree.info$size.node[leaf.idx] <- leaf.counts
  if (wt.pred.accuracy) {
    tree.info$dec.purity <- 0
    tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  out <- list()
  out$tree.info <- tree.info[select.node,]
  
  out$node.feature <- which(t(node.feature) != 0, arr.ind=TRUE)
  out$node.splits <- node.feature[node.feature != 0]
  return(out)
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
    d.left.id <- (x[,split.var] <= split.pt) & parent.id
    d.right.id <- (x[,split.var] > split.pt) & parent.id
    
    node.composition[d.left,] <- d.left.id
    node.composition[d.right,] <- d.right.id
  }
  
  return(node.composition)
}

aggregateNodeFeature <- function(nf) {
  # aggregate list of node feature data returned from each tree
  
  ntree <- length(nf)
  row.offset <- c(0, cumsum(sapply(nf, function(z) max(z[,1])))[-ntree])
  n.rows <- sapply(nf, nrow)
  nf <- do.call(rbind, nf)
  nf[,1] <- nf[,1] + rep(row.offset, times=n.rows)
  return(nf)
}

sdNode <- function(x) {
  sd.node <- ifelse(length(x) == 1, 0, sd(x))
  return(sd.node)
}


## NEW FUNCTIONS FOR HYPERRECTANGLE INTERSECTION
getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent[1] <- 0
  return(parent)
}

parentSplit <- function(tree.info, parents) {
  n <- nrow(tree.info)
  parent.idcs <- parents %% n
  parent.var <- tree.info[parent.idcs, 'split var']
  parent.split <- tree.info[parent.idcs, 'split point']
  parent.var <- c(0, parent.var)
  parent.split <- c(0, parent.split)
  return(cbind(parent.var, parent.split))
}

getAncestorPath <- function(tree.info, node.idx, 
                            path=matrix(0, nrow=node.idx, ncol=3), i=1) {
  path[i,] <- as.matrix(tree.info[node.idx, c('parent.var', 'direction', 'parent.split')])
  parent.idx <- tree.info$parent[node.idx] %% nrow(tree.info)
  if (parent.idx > 0) {
    return(getAncestorPath(tree.info, node.idx=parent.idx, path=path, i=(i+1)))
  } else {
    return(path[1:(i - 1),])
  }
}

path2Binary <- function(path, p, varnames.grp) {
  if (is.null(dim(path))) path <- matrix(path, nrow=1)
  binary <- rep(0L, 2 * p)
  path.vars <- path[,1]
  dup <- duplicated(varnames.grp[path.vars])
  path.vars <- path.vars[!dup]
  path.adj <- (path[!dup, 2] * p)
  path.vars <- path.adj + path.vars
  binary[path.vars] <- 1 * path[!dup,3]
  return(binary)
}
