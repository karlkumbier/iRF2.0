conditionalPred <- function(rfobj, rd.forest, x, y, int, varnames.group=NULL) {
  # Evaluate interaction based on prediction accuracy where predictions are 
  # made using only leaf nodes for which the given interaction falls on the 
  # decision path.
  require(AUC)
  
  if (is.null(varnames.group) & !is.null(colnames(x)))
    varnames.group <- colnames(x)
  
  y.hat <- predIntForest(rfobj, rd.forest, x, y, int, varnames.group)
  if (is.factor(y)) {
    r <- roc(y.hat, y)
    accuracy <- auc(r)
  } else {
    accuracy <- 1 - mean((y - y.hat) ^ 2) / var(y)
  }
  return(accuracy)
}

predIntForest <- function(rfobj, rd.forest, x, y, int, varnames.group) {
  # Predict responses from a RF using only leaf nodes for which a given 
  # interaction falls on the decision path.
  avg.response <- ifelse(is.factor(y), mean(as.numeric(y) - 1), mean(y))
  rd.forest$tree.info$forest.idx <- 1:nrow(rd.forest$tree.info)

  preds <- predict(rfobj, newdata=x, predict.all=TRUE, nodes=TRUE)
  node.mat <- attr(preds, 'nodes')
  interact.nodes <- getInteractNodes(nf=rd.forest$node.feature, 
                                     x.names=varnames.group, int=int)
  
  tree.preds <- sapply(1:rfobj$ntree, predIntTree, 
                          pred.tree=preds$individual, 
                          node.mat=node.mat, 
                          interact.nodes=interact.nodes, 
                          avg.response=avg.response,
                          rd.forest=rd.forest)
  interact.pred <- rowMeans(tree.preds)
} 

predIntTree <- function(pred.tree, node.mat, interact.nodes, 
                        avg.response, rd.forest, tree.idx) {
  # Predict responses from a decision tree using only leaf 
  # nodes for which a given interaction falls on the decision 
  # path.
  require(dplyr)
  # Get node indices for paths with full interaction
  tree.info <- filter(rd.forest$tree.info, tree == tree.idx)
  tree.interact <- interact.nodes[tree.info$forest.idx]
  tree.interact.nodes <- tree.info$node.idx[tree.interact]
  
  # Get predictions of observations that fall in interaction nodes
  tree.nodes <- node.mat[,tree.idx]
  is.interact <- tree.nodes %in% tree.interact.nodes
  tree.preds <- as.numeric(pred.tree[,tree.idx])
  tree.preds[!is.interact] <- avg.response
  return(tree.preds)
}

getInteractNodes <- function(nf, x.names, int) {
  # Determine which leaf nodes contain a given interactions along their 
  # decision paths
  int.split <- strsplit(int, '_')[[1]]
  if (!is.null(x.names)) {
    # group feature matrix by replicated variables
    grp.names <- unique(x.names)
    makeGroup <- function(x, g) apply(as.matrix(x[,x.names == g]), MAR=1, max)
    nf <- sapply(grp.names, makeGroup, x=nf)
    is.interact <- apply(nf[,int.split], MAR=1, sum) == length(int.split)
  } else {
    int.split <- as.numeric(int.split)
    is.interact <- apply(nf[,int.split], MAR=1, sum) == length(int.split)
  }
  return(is.interact)
}


permImportance <- function(rfobj, x, y, int, n.perms=3, 
                           weight=rep(1, length(y)), varnames.group=NULL) {
  # Evaluate interaction based on prediction accuracy of RF when all other 
  # variables are permuted.
  require(iRF)
  
  if (is.null(varnames.group) & !is.null(colnames(x)))
    varnames.group <- colnames(x)
  
  predInt <- function(int) {
    pred.type <- ifelse(rfobj$type == 'regression', 'response', 'prob')
    x.perm <- replicate(n.perms, permuteVars(x, 1:ncol(x)), simplify=FALSE)
    x.fixed <- lapply(x.perm, fixInteract, x=x, int=int, x.names=varnames.group)
    pred.fixed <- lapply(x.fixed, predict, object=rfobj, type=pred.type)
    if (pred.type == 'prob') pred.fixed <- lapply(pred.fixed, function(p) p[,2])  
    interact.score <- mean(sapply(pred.fixed, predAccuracy, y=y, weight=weight))
    return(interact.score)
  }

  interact.scores <- sapply(int, predInt)
  return(interact.scores)
}

predictInteract <- function(rfobj, x, int, n.perms=3, varnames.group=NULL, n.cores=1) {
  # Evaluates predictions from RF when all other variables have been permuted
  require(iRF)

  if (is.null(varnames.group) & !is.null(colnames(x)))
    varnames.group <- colnames(x)
  
  predInt <- function(int) {
    pred.type <- ifelse(rfobj$type == 'regression', 'response', 'prob')
    x.perm <- replicate(n.perms, permuteVars(x, 1:ncol(x)), simplify=FALSE)
    x.fixed <- lapply(x.perm, fixInteract, x=x, int=int, x.names=varnames.group)
    pred.fixed <- lapply(x.fixed, predict, object=rfobj, type=pred.type)
    if (pred.type == 'prob') pred.fixed <- lapply(pred.fixed, function(p) p[,2])
    pred.fixed <- rowMeans(do.call(cbind, pred.fixed))
    return(pred.fixed)
  }

  pred.ints <- mclapply(int, predInt, mc.cores=n.cores)
  pred.ints <- do.call(cbind, pred.ints)
  return(pred.ints)
}

predAccuracy <- function(y.hat, y, weight=weight) {
  # Evaluate prediction accuracy, scaled to interavl [0,1]
  require(AUC)
  if (is.factor(y)) {
    accuracy <- auc(roc(pred, y))
    accuracy <- 2 * (accuracy - 0.5)
  } else {
    if (any(weight < 0)) stop('negative weights')
    weight <- weight / sum(weight)
    accuracy <- 1 - mean((y * weight - y.hat * weight) ^ 2) / var(y * weight)
    if (accuracy < 0) accuracy <- 0
  }
  return(accuracy)
}


fixInteract <- function(x, x.perm, int, x.names) {
  # Returns feature matrix where all variables in int have been fixed to their 
  # raw values and all other variables have been permuted
  int.split <- strsplit(int, '_')[[1]]
  if (!is.null(x.names)) {
    int.split <- which(x.names %in% int.split)
  } else {
    int.split <- as.numeric(int.split)
  }
  
  x.perm[,int.split] <- x[,int.split]
  return(x.perm)
}

permuteVars <- function(x, vars) {
  # Permutes the columns of x corresponding to vars
  x[,vars] <- apply(x[,vars], MAR=2, permute)
  return(x)
}

permute <- function(x) {
  x <- x[sample(length(x))]
  return(x)
}
