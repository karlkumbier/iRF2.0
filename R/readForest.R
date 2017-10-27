readForest <- function(rfobj, x, y=NULL, 
                       return.node.feature=TRUE, 
                       return.node.obs=FALSE,
                       wt.pred.accuracy=FALSE,
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
  #if (is.null(rfobj$obs.nodes)) {
    prf <- predict(rfobj, newdata=x, nodes=TRUE)
    nodes <- attr(prf, 'nodes')
  #} else {
  #  nodes <- rfobj$obs.nodes
  #}
  
  # read leaf node data from each tree in the forest 
  rd.forest <- mclapply(1:ntree, readTree, rfobj=rfobj, x=x, y=y,
                        nodes=nodes,
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
  nobs <- lapply(rd.forest, function(tt) tt$node.obs)
  nobs <- aggregateNodeFeature(nobs)
  
  out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], dims=c(max(nf[,1]), p))
  out$node.obs <- sparseMatrix(i=nobs[,1], j=nobs[,2], dims=c(max(nf[,1]), n))
  return(out)
  
}

readTree <- function(rfobj, k, x, y, nodes, 
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
  tree.info$node.idx <- 1:nrow(tree.info)
  parents <- getParent(tree.info)
  tree.info$parent <- as.integer(parents)
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  # replicate each leaf node in node.feature based on specified sampling.
  select.node <- tree.info$status == -1
  rep.node <- rep(0, nrow(tree.info))
  which.leaf <- nodes[,k]

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
  
  
  # Extract decision paths from leaf nodes as binary sparse matrix
  node.feature <- NULL
  if (return.node.feature) {
    var.nodes <- as.integer(rfobj$forest$bestvar[,k])
    sparse.idcs <- nodeVars(var.nodes,
                            as.integer(length(select.node)),
                            as.integer(p),
                            as.integer(parents),
                            as.integer(select.node),
                            as.integer(rep.node),
                            0L,
                            matrix(0L, nrow=(n * p), ncol=2))
    node.feature <- sparse.idcs[!sparse.idcs[,1] == 0,]
  }
  
  node.obs <- NULL
  if (return.node.obs) {
    id <- match(which.leaf, sort(unique(which.leaf)))
    node.obs <- cbind(id, 1:n)
    node.obs <- node.obs[order(node.obs[,1]),] #TODO: do we need to order?
  }
  
  out <- list()
  out$tree.info <- tree.info
  out$node.feature <- node.feature
  out$node.obs <- node.obs
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
