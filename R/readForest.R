#' Read forestt
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
                       return.node.obs=FALSE,
                       varnames.grp=NULL,
                       first.split=TRUE,
                       weights=rep(1, nrow(x)),
                       n.core=1){
  
  if (is.null(rand.forest$forest))
    stop('No Forest component in the randomForest object')
  varnames.grp <- groupVars(varnames.grp, x)

  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)
  
  ntree <- rand.forest$ntree
  n <- nrow(x)
  p <- length(unique(varnames.grp))
  out <- list()
  
  # Pass observations through RF to determine leaf node membership
  nodes <- NULL
  if (return.node.obs) {
    pred <- predict(rand.forest, newdata=x, nodes=TRUE)
    nodes <- attr(pred, 'nodes')
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
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # Aggregate sparse node level feature matrix
  if (return.node.feature) {
    nf <- lapply(rd.forest, function(tt) tt$node.feature)
    out$node.feature <- do.call(rbind, nf)
  }
  
  # Aggregate sparse node level observation matrix
  if (return.node.obs) {
    nobs <- lapply(rd.forest, function(tt) tt$node.obs)
    nobs <- unlist(nobs, recursive=FALSE)
    col.id <- rep(1:length(nobs), times=out$tree.info$size.node)
    out$node.obs <- sparseMatrix(i=rep(1:n, times=ntree), j=col.id,
                                 dims=c(n, nrow(out$tree.info)))
  } 

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
  
  ntree <- rand.forest$ntree
  n <- nrow(x) 

  # Read tree metadata from forest
  tree.info <- as.data.frame(getTree(rand.forest, k))
  n.node <- nrow(tree.info)
  tree.info$node.idx <- 1:nrow(tree.info)
  tree.info$parent <- getParent(tree.info) %% n.node
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
    unq.leaf <- unique(which.leaf)
    id <- fmatch(which.leaf, sort(unq.leaf))
    node.obs <- c(by(1:n, id, list))
    tree.info$size.node <- sapply(node.obs, length)
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
  children <- c(tree.info[,'left daughter'], tree.info[,'right daughter'])
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

