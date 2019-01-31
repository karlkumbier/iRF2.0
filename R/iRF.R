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

  # Check for depricated input arguments and warn if used
  if (!is.null(interactions.return)) {
    warning('interactions.return is depricated, use iter.return instead')
    iter.return <- interactions.return
    select.iter <- FALSE
  }

  if (!is.null(wt.pred.accuracy))
    warning('wt.pred.accuracy is depricated')

  # Check input attributes for correct format
  if (!class(x) %in% c('data.frame', 'matrix'))
    stop('x must be matrix or data frame')
  if (nrow(x) != length(y) | nrow(xtest) != length(ytest))
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
  }

  # Check all RIT and set to defaul if missing
  if (is.null(rit.param$depth)) rit.param$depth <- 5
  if (is.null(rit.param$ntree)) rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) rit.param$nchild <- 2
  if (is.null(rit.param$class.id)) rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & is.numeric(y)) 
    rit.param$class.cut <- median(y)

  # Register cores for parallel computation 
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)
  

  class.irf <- is.factor(y)
  importance <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
  

  # Iteratively fit random forests, reweighting towards important features.
  rf.list <- list()
  for (iter in 1:n.iter) {

    # Grow Random Forest on full data
    if (verbose) print(paste('iteration = ', iter))
    rf.list[[iter]] <- parRF(x=x, y=y, xtest=xtest, ytest=ytest,
                             mtry.select.prob=mtry.select.prob,
                             ntree=ntree, n.core=n.core)
   
    # Update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]]$importance

    # Report test set accuracy if test set supplied and verbose 
    if (!is.null(xtest) & verbose) 
      print(printAcc(fit=rf.list[[iter]], y=ytest, class.irf=class.irf))
  
  }


  # Select iteration based on prediction error for OOB samples
  selected.iter <- selectIter(rf.list, y=y)
  if (select.iter) {
    iter.return <- selected.iter
    if (verbose) print(paste('selected iteration:', iter.return))
  }

  stability.score <- list()
  importance.score <- list()

  if (is.null(bs.sample)) 
    bs.sample <- replicate(n.bootstrap, bootstrapSample(y), simplify=FALSE)

  for (iter in iter.return) {
    
    # Search for interactions to evaluate across full-data RF 
    rit.param$ntree <- rit.param$ntree * n.bootstrap
    ints.eval <- gRIT(rand.forest=rf.list[[iter]], x=x, y=y,
                      weights=weights, varnames.grp=varnames.grp,
                      rit.param=rit.param, signed=signed,
                      n.core=n.core)
    ints.eval <- ints.eval$int
    rit.param$ntree <- rit.param$ntree / n.bootstrap
       
    # Evaluate stability/importance of recovered interactions
    int.stab <- stabilityScore(fit=rf.list, iter=iter, x=x, y=y,
                               bs.sample=bs.sample, ints.eval=ints.eval,
                               weights=weights, varnames.grp=varnames.grp,
                               rit.param=rit.param, signed=signed, ntree=ntree,
                               n.core=n.core, ...)
    
    stability.score[[iter]] <- int.stab$stab
    importance.score[[iter]] <- int.stab$imp
  
  }
  

  # Generate list for output, selecting iteraction if specified.
  out <- list()
  out$rf.list <- rf.list
  out$selected.iter <- selected.iter
  if (!is.null(iter.return)) {
    out$interaction <- stability.score
    out$importance <- importance.score
  }

  if (length(iter.return) == 1) {
    out$rf.list <- out$rf.list[[iter.return]]
    out$interaction <- out$interaction[[iter.return]]
    out$importance <- out$importance[[iter.return]]
    out$selected.iter <- iter.return
  }

  return(out)
}

sampleClass <- function(y, cl, n) {
  # Take a bootstrap sample from a particular class of observations
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

selectIter <- function(rf.list, y) {
  # Evaluate optimal iteration based on prediction error in OOB samples.
  # For classification, error is given by missclassification rate.
  # For regression, error is given by MSE
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

bootstrapSample <- function(y) {
  # Generate outer layer bootstrap samples
  n <- length(y)
  if (is.factor(y)) {
    # Take bootstrap sample that maintains class balance in full data
    ncl <- table(y)
    class <- as.factor(names(ncl))
    sample.id <- mapply(function(cc, n) sampleClass(y, cc, n), class, ncl)
    sample.id <- unlist(sample.id)
  } else {
    sample.id <- sample(n, replace=TRUE)
  }  
  return(sample.id)
}

