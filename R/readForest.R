readForest <- function(rfobj, x, y=NULL, 
                       return.node.feature=TRUE, 
                       return.node.obs=FALSE,
                       wt.pred.accuracy=FALSE,
                       varnames.grp=NULL,
                       obs.weights=NULL,
                       n.core=1){
  
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (wt.pred.accuracy & is.null(y))
    stop('y required to evaluate prediction accuracy')
  if (is.null(varnames.grp)) varnames.grp <- 1:ncol(x)
  
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
  tree.info$parent <- getParent(tree.info) %% n.node
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  # replicate each leaf node in node.feature based on specified sampling.
  select.node <- tree.info$status == -1
  rep.node <- rep(0, nrow(tree.info))
  which.leaf <- nodes[,k]

  if (return.node.feature) {
    node.feature <- ancestorPath(tree.info, p, varnames.grp=varnames.grp)
    n.path <- sapply(node.feature, length)
    leaf.id <- rep(1:length(node.feature), times=n.path)
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
  tree.info$size.node[leaf.idx] <- leaf.counts
  
  if (wt.pred.accuracy) {
    tree.info$dec.purity <- 0
    tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  tree.info <- tree.info[select.node,]
 
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


ancestorPath <- function(tree.info, p, varnames.grp) {
  
  paths <- getAncestorPath(tree.info, p, varnames.grp)
  paths <- lapply(as.character(which(tree.info$status == -1)), 
                  function(z) which(paths[,z] == 1))
  return(paths)
}

getAncestorPath <- function(tree.info, p,  varnames.grp, node.idx=1, cur.path=NULL) {
  
  if (is.null(cur.path)) cur.path <- rep(0L, 2 * p)
  node.var <- tree.info$`split var`[node.idx]
  node.var.reps <- which(varnames.grp == varnames.grp[node.var])
  
  left.set <- cur.path
  if (all(cur.path[node.var + p] == 0)) left.set[node.var] <- 1L
  left.child <- tree.info$`left daughter`[node.idx]
  
  right.set <- cur.path
  if (all(cur.path[node.var] == 0)) right.set[node.var + p] <- 1L
  right.child <- tree.info$`right daughter`[node.idx]
  
  if (tree.info$status[node.idx] == -1) {
    out <- as.matrix(cur.path)
    colnames(out) <- node.idx
    return(out)
  } else {
    return(cbind(getAncestorPath(tree.info, p, varnames.grp, left.child, left.set),
                 getAncestorPath(tree.info, p, varnames.grp, right.child, right.set)))
  }
}

