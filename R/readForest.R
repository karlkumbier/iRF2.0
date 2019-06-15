#' Read forest
#'
#' Read out metadata from random forest decision paths
#'
#' @param rand.forest an object of class randomForest.
#' @param x numeric feature matrix
#' @param return.node.feature if True, will return sparse matrix indicating
#'  features used on each decision path of the rand.forest
#' @param return.node.obs if True, will return sparse matrix indicating
#'  observations in x that fall in each leaf node of rand.forest.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param first.split if True, splitting threshold will only be evaluated for
#'  the first time a feature is selected.
#' @param n.core number of cores to use. If -1, all available cores are used.
#'
#' @return a list containing the follosing entries
#' \itemize{
#'    \item{tree.info}{data frame of metadata for each leaf node in rand.forest}
#'    \item{node.feature}{optional sparse matrix indicating feature usage on
#'      each decision path}
#'    \item{node.obs}{optional sparse matrix indicating observations appearing
#'      in each leaf node}
#'  }
#'
#' @export
#'
#' @importFrom Matrix Matrix t sparseMatrix rowSums colSums
#' @importFrom data.table data.table rbindlist
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
#' @importFrom dplyr one_of select filter
#' @importFrom fastmatch fmatch
readForest <- function(rand.forest, x, 
                       return.node.feature=TRUE, 
                       return.node.obs=TRUE,
                       varnames.grp=NULL,
                       first.split=TRUE,
                       weights=rep(1, nrow(x)),
                       n.core=1){
  
  # Check for valid input RF
  if (is.null(rand.forest$forest))
    stop('No Forest component in the random forest object')

  # Set feature names if not supplied
  varnames.grp <- groupVars(varnames.grp, x)
  if (is.null(colnames(x)) & class(rand.forest) == 'ranger') {
    colnames(x) <- names(rand.forest$variable.importance)
    varnames.grp <- colnames(x)
  } else if (is.null(colnames(x)) & class(rand.forest) == 'randomForest') {
    colnames(x) <- rownames(rand.forest$importance)
    varnames.grp <- colnames(x)
  }

  # Register cores for parallelization
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)
  
  ntree <- ifelse(class(rand.forest) == 'randomForest',
                  rand.forest$ntree, rand.forest$num.trees)
  
  n <- nrow(x)
  p <- length(unique(varnames.grp))
  out <- list()
  
  # Pass observations through RF to determine leaf node membership
  nodes <- NULL
  if (return.node.obs & class(rand.forest) == 'randomForest') {
    pred <- predict(rand.forest, newdata=x, nodes=TRUE)
    nodes <- attr(pred, 'nodes')
  } else if (class(rand.forest) == 'ranger') {
    pred <- predict(rand.forest, data=x, type='terminalNodes')
    nodes <- pred$predictions
  }
  
  # Split trees across cores to read forest in parallel
  a <- floor(ntree / n.core)
  b <- ntree %% n.core
  ntree.core <- c(rep(a + 1, b), rep(a, n.core - b))
  core.id <- rep(1:n.core, times=ntree.core)

  # Read decision paths across each tree in the forest
  suppressWarnings(
  rd.forest <- foreach(id=1:n.core) %dorng% {
    tree.id <- which(core.id == id)
    readTrees(k=tree.id, rand.forest=rand.forest, x=x, 
              nodes=nodes, varnames.grp=varnames.grp,
              return.node.feature=return.node.feature, 
              return.node.obs=return.node.obs, 
              first.split=first.split, weights=weights)
  })
  rd.forest <- unlist(rd.forest, recursive=FALSE)

  # Aggregate node level metadata
  offset <- cumsum(sapply(rd.forest, function(tt) nrow(tt$tree.info)))
  offset <- c(0, offset[-length(offset)])
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # Aggregate sparse node level feature matrix
  if (return.node.feature) {
    nf.i <- mapply(function(rf, oo) rf$node.feature@i + 1 + oo, rd.forest, offset)
    nf.j <- sapply(rd.forest, function(rf) rep(1:(2 * p), times=diff(rf$node.feature@p)))
    nf.x <- sapply(rd.forest, function(rf) rf$node.feature@x)
    out$node.feature <- sparseMatrix(i=unlist(nf.i), j=unlist(nf.j), x=unlist(nf.x), 
                                     dims=c(nrow(out$tree.info), 2 * p))
  }
  
  # Aggregate sparse node level observation matrix
  if (return.node.obs) {
    col.id <- mapply(function(rf, oo) {
      id <- as.numeric(names(rf$node.obs)) + oo
      nrep <- out$tree.info$size.node[id]
      rep(id, times=nrep)
    }, rd.forest, offset)
    
    nobs <- lapply(rd.forest, function(tt) tt$node.obs)
    nobs <- unlist(nobs, recursive=FALSE)
    out$node.obs <- sparseMatrix(i=unlist(nobs), j=c(col.id),
                                 dims=c(n, nrow(out$tree.info)))
  } 
  
  # Adjust predicted value for classification
  if (class(rand.forest) == 'randomForest')
    if (rand.forest$type == 'classification') 
      out$tree.info$predicted <- out$tree.info$predicted - 1
  
  stopImplicitCluster()
  return(out)
  
}

readTrees <- function(rand.forest, k, x, nodes,
                      varnames.grp=1:ncol(x),
                      return.node.feature=TRUE,
                      return.node.obs=FALSE,
                      first.split=TRUE,
                      weights=rep(1, nrow(x))) {

  out <- lapply(k, readTree, rand.forest=rand.forest, x=x, 
                nodes=nodes,varnames.grp=varnames.grp,
                return.node.feature=return.node.feature, 
                return.node.obs=return.node.obs, 
                first.split=first.split, weights=weights)
  return(out)
}

readTree <- function(rand.forest, k, x, nodes,
                     varnames.grp=1:ncol(x), 
                     return.node.feature=TRUE,
                     return.node.obs=FALSE,
                     first.split=TRUE,
                     weights=rep(1, nrow(x))) {

  n <- nrow(x) 
  ntree <- ifelse(class(rand.forest) == 'randomForest',
                  rand.forest$ntree, rand.forest$num.trees)

  # Read tree metadata from forest
  if (class(rand.forest) == 'randomForest') {
    tree.info <- getTree(rand.forest, k)
  } else if (class(rand.forest) == 'ranger') {
    tree.info <- getTree(rand.forest, k, nodes=nodes)
  } else {
      stop(deparse(substitute(rand.forest)), 
           "is not class ranger of randomForest")
  }

  tree.info$node.idx <- 1:nrow(tree.info)
  tree.info$parent <- getParent(tree.info) %% nrow(tree.info)
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  select.node <- tree.info$status == -1
  
  # Read active features for each decision path
  if (return.node.feature) {
    node.feature <- readFeatures(tree.info, varnames.grp=varnames.grp, 
                                 first.split=first.split)
  }

  tree.info <- filter(tree.info, select.node)
  
  # Read leaf node membership for each observation
  node.obs <- NULL
  if (return.node.obs) {
    which.leaf <- nodes[,k]
    unq.leaf <- sort(unique(which.leaf))
    id <- fmatch(unq.leaf + 1, tree.info$node.idx)
    node.obs <- c(by(1:n, which.leaf, list))
    names(node.obs) <- id
    tree.info$size.node[id] <- sapply(node.obs, length)
  }
    
  out <- list()
  col.remove <- c('left daughter', 'right daughter', 'split var',
                  'split point', 'status')
  out$tree.info <- select(tree.info, -one_of(col.remove))
  out$node.feature <- node.feature
  out$node.obs <- node.obs
  return(out)
}

getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  children <- c(tree.info$`left daughter`, tree.info$`right daughter`)
  parent <- fmatch(1:nrow(tree.info), children)
  parent[1] <- 0
  return(parent)
}

readFeatures <- function(tree.info, varnames.grp, 
                         split.pt=FALSE, first.split=TRUE) {

  # Pre-allocate variables for path ancestry 
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)
  nlf <- sum(tree.info$status == -1)
  cur.path <- rep(0L, 2 * p + 1)

  # Recursively extract path info for all leaf nodes 
  paths <- ancestorPath(tree.info, varnames.grp, varnames.unq, p=p,
                        cur.path=cur.path, first.split=first.split)

  # Generate sparse matrix of decision path feature selection
  paths <- Matrix(unlist(paths), nrow=nlf, byrow=TRUE, sparse=TRUE)
  rownames(paths) <- paths[,ncol(paths)]
  paths <- paths[,1:(2 * p)]

  # Reorder leaf nodes according to tree.info
  idlf <- tree.info$node.idx[tree.info$status == -1]
  paths <- paths[fmatch(idlf, rownames(paths)),]
  return(paths)
}



ancestorPath <- function(tree.info, varnames.grp, varnames.unq, p, 
                         cur.path, node.idx=1, first.split=TRUE) {

  # Return path vector if current node is leaf
  if (tree.info$status[node.idx] == -1) {
    cur.path[length(cur.path)] <- node.idx
    return(cur.path)
  }

  # Extract split feature and threshold for current node
  id <- tree.info$`split var`[node.idx]
  sp <- tree.info$`split point`[node.idx]
  node.var <- which(varnames.grp[id] == varnames.unq)
  
  # Generate vector indicating threshold for first variable split 
  left.child <- tree.info$`left daughter`[node.idx]
  left.set <- cur.path
  if (first.split) {
    if (cur.path[node.var + p] == 0) left.set[node.var] <- sp
  } else {
    left.set[node.var] <- sp
  }
  
  right.child <- tree.info$`right daughter`[node.idx]
  right.set <- cur.path
  if (first.split) {
    if (cur.path[node.var] == 0) right.set[node.var + p] <- sp
  } else {
    right.set[node.var + p] <- sp
  }
  
  # Get acncestor path info for child nodes and combine
  lpath <- ancestorPath(tree.info, varnames.grp, varnames.unq, 
                        p, left.set, left.child, first.split)
    
  rpath <- ancestorPath(tree.info, varnames.grp, varnames.unq,
                        p, right.set, right.child, first.split)
  return(list(lpath, rpath))
}

