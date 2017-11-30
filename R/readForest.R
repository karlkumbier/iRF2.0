readForest <- function(rfobj, x, y=NULL, 
                       return.node.feature=TRUE, 
                       return.node.obs=FALSE,
                       wt.pred.accuracy=FALSE,
                       varnames.grp=1:ncol(x),
                       obs.weights=NULL,
                       n.core=1){
  
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (wt.pred.accuracy & is.null(y))
    stop('y required to evaluate prediction accuracy')
  
  ntree <- rfobj$ntree
  n <- nrow(x)
  p <- ncol(x)
  out <- list()
  
  # Determine leaf nodes for observations in x
  prf <- predict(rfobj, newdata=x, nodes=TRUE)
  nodes <- attr(prf, 'nodes')
  
  # read leaf node data from each tree in the forest 
  rd.forest <- mclapply(1:ntree, readTree, rfobj=rfobj, x=x, y=y,
                        nodes=nodes,
                        varnames.grp=varnames.grp,
                        return.node.feature=return.node.feature,
                        return.node.obs=return.node.obs,
                        wt.pred.accuracy=wt.pred.accuracy,
                        obs.weights=obs.weights,
                        mc.cores=n.core)
  
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # aggregate sparse feature matrix across forest
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  nf <- aggregateNodeFeature(nf)
  
  # aggregate sparse observation matrix across forest
  if (return.node.obs) nobs <- lapply(rd.forest, function(tt) tt$node.obs)
  if (return.node.obs) nobs <- aggregateNodeFeature(nobs)
  
  out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], dims=c(max(nf[,1]), 2 * p))
  if (return.node.obs)
    out$node.obs <- sparseMatrix(i=nobs[,1], j=nobs[,2], dims=c(max(nf[,1]), n))
  return(out)
  
}

readTree <- function(rfobj, k, x, y, nodes,
                     varnames.grp=1:ncol(x), 
                     return.node.feature=TRUE,
                     return.node.obs=FALSE,
                     wt.pred.accuracy=FALSE, 
                     obs.weights=NULL) {
  n <- nrow(x) 
  p <- ncol(x)
  if (is.factor(y)) y <- as.numeric(y) - 1
  ntree <- rfobj$ntree

  # Read tree metadata from forest
  tree.info <- as.data.frame(getTree(rfobj, k))
  n.node <- nrow(tree.info)
  tree.info$node.idx <- 1:nrow(tree.info)
  parents <- getParent(tree.info)
  tree.info$parent <- as.integer(parents) %% n.node
  parent.splits <- parentSplit(tree.info, parents)
  tree.info$parent.var <- parent.splits[,1]
  tree.info$parent.split <- parent.splits[,2]
  tree.info$direction <- ifelse(parents > n.node, 1, 0)
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  # replicate each leaf node in node.feature based on specified sampling.
  select.node <- tree.info$status == -1
  rep.node <- rep(0, nrow(tree.info))
  which.leaf <- nodes[,k]

  if (return.node.feature) {
    ancestors <- lapply(which(select.node), getAncestorPath, tree.info=tree.info)
    node.feature <- lapply(ancestors, path2Idx, varnames.grp=varnames.grp)
    n.path <- sapply(node.feature, length)
    leaf.id <- rep(1:length(ancestors), times=n.path)
    node.feature <- cbind(leaf.id, unlist(node.feature))
  }

  # if specified, set node counts based on observation weights
  if (is.null(obs.weights)) {
    leaf.counts <- table(which.leaf)
    leaf.idx <- as.integer(names(leaf.counts))
  } else {
    leaf.counts <- c(by(obs.weights, which.leaf, sum))
    leaf.idx <- as.integer(names(leaf.counts))
  }
  
  # if specified, calculate purity of each node
  if (wt.pred.accuracy) leaf.sd <- c(by(y, which.leaf, sdNode))

  tree.info <- select(tree.info, prediction, node.idx, 
                      parent, tree, size.node, status)
  tree.info$size.node[leaf.idx] <- leaf.counts
  
  if (wt.pred.accuracy) {
    tree.info$dec.purity <- 0
    tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  tree.info <- tree.info[leaf.idx,]
  rep.node[leaf.idx] <- 1
 
  node.obs <- NULL
  if (return.node.obs) {
    id <- match(which.leaf, sort(unique(which.leaf)))
    node.obs <- cbind(id, 1:n)
    node.obs <- node.obs[order(node.obs[,1]),]
  }
  
  out <- list()
  out$tree.info <- tree.info
  out$node.feature <- node.feature
  out$node.obs <- node.obs
  return(out)
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

getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent[1] <- 0
  return(parent)
}

parentSplit <- function(tree.info, parents) {
  # Get parent split variable and point
  parent.idcs <- parents %% nrow(tree.info)
  parent.var <- c(0, tree.info[parent.idcs, 'split var'])
  parent.split <- c(0, tree.info[parent.idcs, 'split point'])
  return(cbind(parent.var, parent.split))
}

getAncestorPath <- function(tree.info, node.idx, path=NULL, i=1) {
  # Traverse tree to generate ancestor info for each leaf
  if (is.null(path)) path <- matrix(0, nrow=node.idx, ncol=2)
  cur.node <- tree.info[node.idx, c('parent.var', 'direction')]
  path[i,] <- as.matrix(cur.node)
  parent.idx <- tree.info$parent[node.idx]
  if (parent.idx > 0) {
    return(getAncestorPath(tree.info, node.idx=parent.idx, path=path, i=(i+1)))
  } else {
    return(path[1:(i - 1),])
  }
}

path2Idx <- function(path, varnames.grp) {
  # represent each path as a length 2 * p sparse binary vector indicating 
  # splitting feature and direction
  p <-  length(varnames.grp) # we are not grouping features here
  if (is.null(dim(path))) path <- matrix(path, nrow=1)
  path.vars <- path[,1]

  # Paths are ordered from leaf to root. Remove all but last 
  # occurance of duplicate
  dup <- rev(duplicated(rev(varnames.grp[path.vars])))
  path.vars <- path.vars[!dup]
  path.adj <- (path[!dup, 2] * p)
  path.vars <- path.adj + path.vars
  return(path.vars)
}
