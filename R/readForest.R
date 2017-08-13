readForest <- function(rfobj, x, y=NULL,
                       return.node.feature=TRUE,
                       wt.pred.accuracy=FALSE,
                       obs.weights=NULL,
                       n.core=1){

  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (wt.pred.accuracy & is.null(y))
    stop('y required to evaluate prediction accuracy')

  ntree <- rfobj$ntree
  p <- ncol(x)
  n <- nrow(x)
  out <- list()

  # read leaf node data from each tree in the forest 
  rd.forest <- mclapply(1:ntree, function(tt) {
                          tryCatch({
                            readTree(rfobj=rfobj, k=tt, x=x, y=y,
                              return.node.feature=return.node.feature,
                              wt.pred.accuracy=wt.pred.accuracy,
                              obs.weights=obs.weights
                              )
                          }, error=function(e) print(e),
                          warning=function(w) print(w))
                       }, mc.cores=n.core)

  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  # aggregate sparse feature matrix across forest
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  #id.rm <- sapply(nf, is.null)
  #nf <- nf[!id.rm]
  nf <- aggregateNodeFeature(nf)
  out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], dims=c(max(nf[,1]), 2 * p))
  return(out)

}

readTree <- function(rfobj, k, x, y, 
                     return.node.feature=TRUE, 
                     wt.pred.accuracy=FALSE,
                     obs.weights=NULL) {
  
  n <- nrow(x)
  p <- ncol(x)
  ntree <- rfobj$ntree
  tree.info <- as.data.frame(getTree(rfobj, k))
  tree.info$node.idx <- 1:nrow(tree.info)
  parents <- getParent(tree.info)
  tree.info$parent <- as.integer(parents)
  tree.info$parent.var <- getParentVar(tree.info, tree.info$parent)
  tree.info$direction <- ifelse(tree.info$parent > nrow(tree.info), 1, 0)
  tree.info$direction[1] <- 0
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  leaf.idcs <- tree.info$node.idx[tree.info$status == -1]
  ancestors <- lapply(leaf.idcs, getAncestorPath, tree.info=tree.info)
  # TODO: only take earliest splitting node for variables that are used more
  # than once
 # if (TRUE) {
 #   getUniqueFeature <- function(an) {
 #     unq <- unique(an[,1])
 #     unq.idcs <- sapply(unq, function(u) max(which(an[,1] == u)))
 #     print(unq.idcs)
 #     return(an[unq.idcs,])
 #   }
 #   ancestors <- lapply(ancestors, getUniqueFeature)
 # }
  node.feature <- sapply(ancestors, path2Binary, p=p)
  
  select.node <- tree.info$status == -1
  rep.node <- rep(0L, nrow(tree.info))
  tree.info <- select(tree.info, prediction, node.idx, parent, tree, size.node)
  
  if (is.null(rfobj$obs.nodes)) {
    # if nodes not tracked, pass data through forest to get leaf
    # counts
    # TODO: node weighted sampling here
    fit.data <- passData(rfobj, x, tree.info, k)
    leaf.counts <- rowSums(fit.data[select.node,])
    which.leaf <- apply(fit.data[select.node,], MAR=2, which)
    leaf.idx <- as.integer(which(select.node))
    if (wt.pred.accuracy) leaf.sd <- c(by(y, which.leaf, sdNode))
  } else {
    if (is.null(obs.weights)) {
      leaf.counts <- table(rfobj$obs.nodes[,k])
      leaf.idx <- as.integer(names(leaf.counts))
      if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
    } else {
      leaf.counts <- c(by(obs.weights, rfobj$obs.nodes[,k], sum))
      leaf.idx <- as.integer(names(leaf.counts))
      if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
    }
    
  }
  
  tree.info$size.node[leaf.idx] <- leaf.counts
  if (wt.pred.accuracy) {
    tree.info$dec.purity <- 0
    tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  tree.info <- tree.info[select.node,]
  out <- list()
  out$tree.info <- tree.info
  nf <- t(node.feature) == 1
  nf <- which(nf, arr.ind=TRUE)
  out$node.feature <- nf

  return(out)
}

readTree2 <- function(rfobj, k, x, y, return.node.feature, wt.pred.accuracy,
                     obs.weights) {
  n <- nrow(x)
  p <- ncol(x)
  ntree <- rfobj$ntree

  # Read tree level data from RF
  out <- list()
  out$tree.info <- as.data.frame(getTree(rfobj, k))
  out$tree.info$node.idx <- 1:nrow(out$tree.info)
  parents <- getParent(out$tree.info)
  out$tree.info$parent <- as.integer(parents)
  out$tree.info$tree <- as.integer(k)
  out$tree.info$size.node <- 0L

  # replicate each leaf node in node.feature based on specified
  # sampling.
  select.node <- out$tree.info$status == -1
  rep.node <- rep(0, nrow(out$tree.info))
  out$tree.info <- select(out$tree.info, prediction, node.idx, parent, tree, size.node)

  if (is.null(rfobj$obs.nodes)) {
    # if nodes not tracked, pass data through forest to get leaf
    # counts
    # TODO: node weighted sampling here
    fit.data <- passData(rfobj, x, out$tree.info, k)
    leaf.counts <- rowSums(fit.data[select.node,])
    which.leaf <- apply(fit.data[select.node,], MAR=2, which)
    leaf.idx <- as.integer(which(select.node))
    if (wt.pred.accuracy) leaf.sd <- c(by(y, which.leaf, sdNode))
  } else {
    if (is.null(obs.weights)) {
      leaf.counts <- table(rfobj$obs.nodes[,k])
      leaf.idx <- as.integer(names(leaf.counts))
      if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
    } else {
      leaf.counts <- c(by(obs.weights, rfobj$obs.nodes[,k], sum))
      leaf.idx <- as.integer(names(leaf.counts))
      if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
    }

  }

  out$tree.info$size.node[leaf.idx] <- leaf.counts
  if (wt.pred.accuracy) {
    out$tree.info$dec.purity <- 0
    out$tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }

  out$tree.info <- out$tree.info[select.node,]
  rep.node[select.node] <- 1

  # Extract decision paths from leaf nodes as binary sparse
  # matrix
  if (return.node.feature) {
    row.offset <- 0
    var.nodes <- as.integer(rfobj$forest$bestvar[,k])
    total.rows <- n #sum(rep.node[select.node])
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

getParent2 <- function(tree.info) {
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

getParentVar <- function(tree.info, parents) {
  n <- nrow(tree.info)
  parent.var <- tree.info[parents %% n, 'split var']
  parent.var <- c(0, parent.var)
  return(parent.var)
}

getAncestorPath <- function(tree.info, node.idx, path=matrix(0, nrow=node.idx, ncol=2), i=1) {
  
  path[i,] <- as.matrix(tree.info[node.idx, c('parent.var', 'direction')])
  parent.idx <- tree.info$parent[node.idx] %% nrow(tree.info)
  if (parent.idx > 0) {
    return(getAncestorPath(tree.info, node.idx=parent.idx, path=path, i=(i+1)))
  } else {
    return(path[1:(i - 1),])
  }
}

path2Binary <- function(path, p) {
  binary <- rep(0L, 2 * p)
  path.vars <- path[,1]
  path.adj <- (path[,2] * p)
  path.vars <- path.adj + path.vars
  binary[path.vars] <- 1L
  return(binary)
}
