# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, 
                ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1, ncol(x)),
                iter.return=NULL, 
                select.iter=ifelse(is.null(iter.return), TRUE, FALSE),
                rit.param=list(depth=5, ntree=500, 
                               nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL), 
                varnames.grp=colnames(x), 
                n.bootstrap=1,
                weights=rep(1, nrow(x)),
                signed=TRUE,
                bs.sample=NULL,
                verbose=TRUE,
                interactions.return=NULL,
                wt.pred.accuracy=NULL,
                ...) {
 
  # Check for depricated arguments
  if (!is.null(interactions.return)) {
    warning('interactions.return is depricated, use iter.return instead')
    iter.return <- interactions.return
    select.iter <- FALSE
  }
  if (!is.null(wt.pred.accuracy))
    warning('wt.pred.accuracy is depricated')

  # Check input attributes for correct format
  require(doRNG, quiet=TRUE)
  if (!class(x) %in% c('data.frame', 'matrix'))
    stop('x must be matrix or data frame')
  if (nrow(x) != length(y))
    stop('x and y must contain the same number of observations')
  if (ncol(x) < 2 & (!is.null(iter.return) | select.iter))
    stop('cannot find interaction - x has less than two columns!')
  if (any(iter.return > n.iter))
    stop('selected iteration greater than n.iter')
  if (length(mtry.select.prob) != ncol(x))
    stop('length mtry.select.prob must equal number of features')
  if (length(weights) != nrow(x))
    stop('length weights differs from # training observations')
  if (!is.null(xtest)) {
    if (ncol(xtest) != ncol(x)) 
      stop('training/test data must have same number of features')
    if (is.null(ytest))
      stop('test set responses not indicated')
    if (nrow(xtest) != length(ytest))
      stop('xtest and ytest must contain the same number of observations')
  }

  # Check all RIT and set to defaul if missing
  if (is.null(rit.param$depth)) rit.param$depth <- 5
  if (is.null(rit.param$ntree)) rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) rit.param$nchild <- 2
  if (is.null(rit.param$class.id)) rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & is.numeric(y)) 
    rit.param$class.cut <- median(y)
  

  class.irf <- is.factor(y)
  importance <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
  
  # Fit a series of iteratively re-weighted RFs 
  rf.list <- list()  
  for (iter in 1:n.iter) {
    
    # Grow Random Forest on full data
    if (verbose) print(paste('iteration = ', iter))
    rf.list[[iter]] <- parRF(x, y, xtest, ytest, ntree=ntree, n.core=n.core, 
                             mtry.select.prob=mtry.select.prob, ...)
    
    # Update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]]$importance

    # Evaluate test set error if supplied
    if (!is.null(xtest) & verbose) 
      print(printAcc(rf.list[[iter]], ytest, class.irf))
  
  }
  

  # Select iteration to return interactions based on OOB error
  selected.iter <- selectIter(rf.list, y=y)
  if (select.iter) iter.return <- selected.iter
  

  # Run gRIT across full data RF and outer layer boostrap stability analysis of
  # recovered interactions 
  importance <- list()
  if (is.null(bs.sample)) bs.sample <- lreplicate(n.bootstrap, bsSample(y))
  for (iter in iter.return) {
    
    # Evaluate interactions in full data random forest 
    if (verbose) cat('finding interactions...\n')
    rit.param$ntree <- rit.param$ntree * n.bootstrap
    ints.eval <- gRIT(rf.list[[iter]], x=x, y=y,
                      weights=weights,
                      varnames.grp=varnames.grp,
                      rit.param=rit.param,
                      signed=signed,
                      n.core=n.core)
    ints.eval <- ints.eval$int
    rit.param$ntree <- rit.param$ntree / n.bootstrap

    if (verbose) cat('evaluating interactions...\n')
    # Find interactions across bootstrap replicates

    if (iter == 1) rf.weight <- rep(1, ncol(x))
    if (iter > 1) rf.weight <- rf.list[[iter - 1]]$importance
    importance[[iter]] <- stabilityScore(x, y, rf.weight, bs.sample,
                                         ints.eval=ints.eval, ntree=ntree,
                                         weights=weights, signed=signed,
                                         varnames.grp=varnames.grp,
                                         rit.param=rit.param, n.core=n.core, 
                                         ...)
  
  }
  
  # Combine reults for return
  out <- list()
  out$rf.list <- rf.list
  out$selected.iter <- selected.iter
  if (!is.null(iter.return)) {
    out$interaction <- importance
  }

  if (length(iter.return) == 1) {
    iter.wt <- iter.return - 1
    if (iter.return > 1) out$weights <- out$rf.list[[iter.wt]]$importance 
    out$rf.list <- out$rf.list[[iter.return]]
    out$interaction <- importance[[iter.return]]
    out$selected.iter <- iter.return
  }

  return(out)
}


selectIter <- function(rf.list, y) {
  # Evaluate optimal iteration based on prediction error in OOB samples.
  #   For classification: accuracy.
  #   For regression: MSE.
  predicted <- lapply(rf.list, function(z) as.numeric(z$predicted))
  
  if (is.factor(y)) {
    predicted <- lapply(predicted, '-', 1)
    y <- as.numeric(y) - 1
    eFun <- function(y, yhat) sum(y == 1 & yhat == 0) + sum(y == 0 & yhat == 1)
  } else {
    eFun <- function(y, yhat) mean((yhat - y) ^ 2, na.rm=TRUE)
  }
  
  error <- sapply(predicted, eFun, y=y)
  min.err <- min(error)
  id.select <- max(which(error == min.err))
  return(id.select)
}

sampleClass <- function(y, cl, n) {
  # Take a bootstrap sample from a particular class of observations
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

bsSample <- function(y) {
  # Generate outer layer bootstrap samples
  n <- length(y)
  if (is.factor(y)) {
    # Take bootstrap sample that maintains class balance of full data
    ncl <- table(y)
    class <- as.factor(names(ncl))
    sample.id <- mapply(function(cc, n) sampleClass(y, cc, n), class, ncl)
    sample.id <- unlist(sample.id)
  } else {
    sample.id <- sample(n, replace=TRUE)
  }
  return(sample.id)
}

