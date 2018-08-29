# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1, ncol(x)),
                interactions.return=NULL, 
                rit.param=list(depth=5, ntree=500, 
                               nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL), 
                varnames.grp=NULL, 
                n.bootstrap=1,
                select.iter=FALSE,
                get.prevalence=TRUE,
                int.sign=TRUE,
                int.subs=TRUE,
                verbose=TRUE,
                bootstrap.path=NULL,
                ...) {
  
  
  if (ncol(x) < 2 & !is.null(interactions.return))
    stop('cannot find interaction - X has less than two columns!')
  
  if (any(interactions.return > n.iter))
    stop('interaction iteration to return greater than n.iter')
  
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  if (!class.irf & is.null(rit.param$class.cut)) 
    rit.param$class.cut <- median(y)
  importance <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
   
  rf.list <- list()
  if (!is.null(interactions.return) | select.iter) {
    stability.score <- list()
    if (get.prevalence) prevalence.score <- list()
  }


  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  registerDoMC(n.core)
  for (iter in 1:n.iter) {
    
    # Grow Random Forest on full data
    print(paste('iteration = ', iter))
    rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dorng% {
                                 randomForest(x, y, 
                                              xtest, ytest, 
                                              ntree=ntree.id[i], 
                                              mtry.select.prob=mtry.select.prob, 
                                              keep.forest=TRUE,
                                              ...)
                               }
    
    # Update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]]$importance

    if (!is.null(xtest) & class.irf & verbose) {
      auroc <- auc(roc(rf.list[[iter]]$test$votes[,2], ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    } else if (!is.null(xtest) & verbose) {
      pct.var <- 1 - mean((rf.list[[iter]]$test$predicted - ytest) ^ 2) / var(ytest)
      pct.var <- max(pct.var, 0)
      print(paste('% var explained:', pct.var * 100))
    }
  }
  
  # Select iteration for which to return interactions based on minimizing 
  # prediction error on OOB samples
  if (select.iter) {
    interactions.return <- selectIter(rf.list, y=y)
    if (verbose) print(paste('selected iter:', interactions.return))
  }
  
  for (iter in interactions.return) {
    
    # Find interactions across bootstrap replicates
    if (verbose) cat('finding interactions ... ')
    
    interact.list <- list()
    if (get.prevalence) prev.list <- list()

    for (i.b in 1:n.bootstrap) { 
      if (class.irf) {
        # Take bootstrap sample that maintains class balance in full data
        ncl <- table(y)
        class <- as.factor(names(ncl))
        sample.id <- mapply(function(cc, n) sampleClass(y, cc, n), class, ncl)
        sample.id <- unlist(sample.id)
      } else {
        sample.id <- sample(n, replace=TRUE)
      }
      
      # Use feature weights from current iteraction of full data RF
      if (iter == 1) 
        mtry.select.prob <- rep(1, ncol(x))
      else
        mtry.select.prob <- rf.list[[iter - 1]]$importance

      # Fit random forest on bootstrap sample
      rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                      .multicombine=TRUE, .packages='iRF') %dorng% {
                        randomForest(x[sample.id,], y[sample.id], 
                                     xtest, ytest, 
                                     ntree=ntree.id[i], 
                                     mtry.select.prob=mtry.select.prob, 
                                     keep.forest=TRUE, 
                                     ...)
                      }

      # Write out bootstrap RFs for later processing
      if (!is.null(bootstrap.path)) {
        dir.create(bootstrap.path, showWarnings=FALSE)
        out.file <- paste0(bootstrap.path, 'bs_iter', 
                           iter, '_b', i.b, '.Rdata')
      } else {
        out.file <- NULL
      }
      
     
      # Run generalized RIT on rf.b to learn interactions
      ints <- generalizedRIT(rf=rf.b, x=x, y=y,
                             varnames.grp=varnames.grp,
                             rit.param=rit.param,
                             get.prevalence=get.prevalence,
                             int.sign=int.sign,
                             int.subs=int.subs,
                             out.file=out.file,
                             n.core=n.core)
      
      interact.list[[i.b]] <- ints$int
      if (get.prevalence) prev.list[[i.b]] <- ints$prev

      rm(rf.b)       
    }
    
    # Calculate stability scores of interactions
    stability.score[[iter]] <- summarizeInteract(interact.list)
    if (get.prevalence)  prevalence.score[[iter]] <- summarizePrev(prev.list)
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) out$interaction <- stability.score
  if (get.prevalence) out$prevalence <- prevalence.score

  if (select.iter) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    if (get.prevalence) out$prevalence <- out$prevalence[[interactions.return]]
  }

  return(out)
}


summarizeInteract <- function(store.out, local=FALSE){
  # Aggregate interactions across bootstrap samples
  if (local) store.out <- lapply(store.out, names)

  n.bootstrap <- length(store.out)
  store <- unlist(store.out)
  
  if (length(store) >= 1){
    int.tbl <- sort(c(table(store)), decreasing=TRUE)
    int.tbl <- int.tbl / n.bootstrap
  } else {
    return(list(interaction=numeric(0), prevalence=numeric(0)))
  }
  
  out <- int.tbl
  return(out)
}

sampleClass <- function(y, cl, n) {
  # Sample indices specific to a given class
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

selectIter <- function(rf.list, y) {
  # Evaluate optimal iteration based on prediction error in OOB samples
  predicted <- lapply(rf.list, function(z) 
                      as.numeric(z$predicted) - is.factor(y))
  
  if (is.factor(y)){
    y <- as.numeric(y) - 1
    errFun <- function(y, yhat) {
     out <- mean(y == 1) * sum(y == 0 & yhat == 1) + 
       mean(y == 0) * sum(y == 1 & yhat == 0)
     return(out)
    }
  } else {
    errFun <- function(y, yhat) mean((yhat - y) ^ 2, na.rm=TRUE)
  }
  
  error <- sapply(predicted, errFun, y=y)
  min.err <- min(error)
  id.select <- max(which(error == min.err))
  return(id.select)
} 
