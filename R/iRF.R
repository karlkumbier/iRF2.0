# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1/ncol(x), ncol(x)), 
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                rit.param=list(depth=5, ntree=500, nchild=2, class.id=1, 
                               class.cut=NULL, class.qt=0.5), 
                varnames.grp=NULL, 
                n.bootstrap=20,
                bootstrap.forest=TRUE,
                select.iter=FALSE,
                verbose=TRUE,
                ...) {
  
  
  if (!is.matrix(x) | (!is.null(xtest) & !is.matrix(xtest)))
    stop('either x or xtest is not a matrix !')
  
  if (!is.numeric(x) | (!is.null(xtest) & !is.numeric(xtest)))
    stop('either x or xtest is not a numeric matrix!')
  
  if (ncol(x) < 2 & !is.null(interactions.return))
    stop('cannot find interaction - X has less than two columns!')
  
  if (any(interactions.return > n.iter))
    stop('interaction iteration to return greater than n.iter')
  
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  imp <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
  registerDoMC(n.core)  
  
  rf.list <- list()
  if (!is.null(interactions.return) | select.iter) {
    interact.list <- list()
    stability.score <- list()
  }
  
  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  for (iter in 1:n.iter) {
    
    # Grow Random Forest on full data
    print(paste('iteration = ', iter))
    rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dopar% {
                                 randomForest(x, y, 
                                              xtest, ytest, 
                                              ntree=ntree.id[i], 
                                              mtry.select.prob=mtry.select.prob, 
                                              keep.forest=TRUE,
                                              track.nodes=TRUE, 
                                              ...)
                               }
    
    ## update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]]$importance[,imp]
    
    if (!is.null(xtest) & class.irf & verbose){
      ypred <- predict(rf.list[[iter]], newdata=xtest, type='prob')[,2]
      auroc <- auc(roc(ypred, ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    } else if (!is.null(xtest) & verbose) {
      ypred <- predict(rf.list[[iter]], newdata=xtest)
      pct.var <- 1 - mean((ypred - ytest) ^ 2) / var(ytest)
      pct.var <- max(pct.var, 0)
      print(paste('% var explained:', pct.var * 100))
    }
  }
  
  # Select iteration for which to return interactions based on minimizing 
  # prediction error on OOB samples
  if (select.iter) interactions.return <- selectIter(rf.list, y=y)
  
  for (iter in interactions.return) {
    
    # Find interactions across bootstrap replicates
    if (verbose) cat('finding interactions ... ')
    interact.list.b <- list()      
    
    for (i.b in 1:n.bootstrap) { 
      
      if (class.irf) {
        # Take bootstrap sample that maintains class balance in full data
        n.class <- table(y)
        class <- as.factor(names(n.class))
        sample.id <- mapply(function(cc, nn) sampleClass(y, cc, nn), 
                            class, n.class)
        sample.id <- unlist(sample.id)
      } else {
        sample.id <- sample(n, replace=TRUE)
      }
      
      if (bootstrap.forest) { 
        
        # use feature weights from current iteration of full data RF
        if (iter == 1)
          mtry.select.prob <- rep(1, p)
        else
          mtry.select.prob <- rf.list[[iter - 1]]$importance[,imp]
        
        # fit random forest on bootstrap sample
        rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                        .multicombine=TRUE, .packages='iRF') %dopar% {
                          randomForest(x[sample.id,], y[sample.id], 
                                       xtest, ytest, 
                                       ntree=ntree.id[i], 
                                       mtry.select.prob=mtry.select.prob, 
                                       keep.forest=TRUE, 
                                       track.nodes=TRUE)#, 
                          #...)
                        }
      } else {
        rf.b <- rf.list[[iter]]
      }
      
      # run generalized RIT on rf.b to learn interactions
      ints <- generalizedRIT(rf=rf.b, 
                             x=x[sample.id,], y=y[sample.id],
                             wt.pred.accuracy=wt.pred.accuracy,
                             class.irf=class.irf,
                             imp=imp,
                             varnames.grp=varnames.grp,
                             rit.param=rit.param,
                             n.core=n.core)
      
      interact.list.b[[i.b]] <- ints
      rm(rf.b)       
    }
    
    # calculate stability scores of interactions
    if (!is.null(varnames.grp))
      varnames.new <- unique(varnames.grp)
    else if (!is.null(colnames(x)))
      varnames.new <- colnames(x)
    else
      varnames.new <- 1:ncol(x)
    
    interact.list[[iter]] <- interact.list.b 
    summary.interact <- summarizeInteract(interact.list[[iter]], 
                                          varnames=varnames.new)
    stability.score[[iter]] <- summary.interact
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) out$interaction <- stability.score
  
  if (select.iter) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    out$opt.k <- interactions.return
    out$weights <- rf.list[[interactions.return]]$importance[,imp]
  }
  
  return(out)
}

generalizedRIT <- function(rf, x, y, wt.pred.accuracy, class.irf, 
                           imp, varnames.grp, rit.param, n.core) {
  
  # Extract decision paths from rf as sparse binary matrix to be passed to RIT
  rforest <- readForest(rf, x=x, y=y, wt.pred.accuracy=wt.pred.accuracy, 
                        n.core=n.core)                                                        
  
  # Select class specific leaf nodes
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  
  if (class.irf) {
    class.id <- rit.param$class.id
    select.leaf.id <- rforest$tree.info$prediction == as.numeric(class.id) + 1
  } else if (!is.null(rit.param$class.cut)) {
    select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
    select.leaf.id <- rforest$tree.info$prediction > rit.param$class.cut
    
    if (sum(select.leaf.id) < 2) {
      warning('fewer than 2 leaf nodes with prediction > class.cut, using median')
      select.leaf.id <- rforest$tree.info$prediction > median(y)
    }
    
  } else if (!is.null(rit.param$class.qt)) {
    select.leaf.id <- rforest$tree.info$prediction > 
      quantile(y, rit.param$class.qt)
  }       
  
  rforest <- subsetReadForest(rforest, select.leaf.id)
  nf <- rforest$node.feature
  
  if (wt.pred.accuracy) {
    wt <- rforest$tree.info$size.node * rforest$tree.info$dec.purity
  } else {
    wt <- rforest$tree.info$size.node
  }         
  
  rm(rforest)
  
  if (sum(select.leaf.id) < 2){
    return(character(0))
  } else {  
    # group features if specified
    if (!is.null(varnames.grp)) nf <- groupFeature(nf, grp=varnames.grp)
   
    interactions <- RIT(nf, weights=wt, depth=rit.param$depth,
                        n_trees=rit.param$ntree, branch=rit.param$nchild,
                        n_cores=n.core)                                                 
    interactions$Interaction <- gsub(' ', '_', interactions$Interaction)
    return(interactions)
  }                   
}

subsetReadForest <- function(rforest, subset.idcs) {
  
  # Subset nodes from readforest output 
  if (!is.null(rforest$node.feature)) 
    rforest$node.feature <- rforest$node.feature[subset.idcs,]
  
  if(!is.null(rforest$tree.info))
    rforest$tree.info <- rforest$tree.info[subset.idcs,]
  
  return(rforest)
}

groupFeature <- function(node.feature, grp){
  # Group feature level data in node.feature 
  sparse.mat <- is(node.feature, 'Matrix')
  grp.names <- unique(grp)
  makeGroup <- function(x, g) apply(as.matrix(x[,grp == g]), MAR=1, max) 
  node.feature.new <- sapply(grp.names, makeGroup, x=node.feature)
  
  if (sparse.mat) node.feature.new <- Matrix(node.feature.new, sparse=TRUE)
  colnames(node.feature.new) <- grp.names
  
  return(node.feature.new)
}



summarizeInteract <- function(store.out, varnames=NULL){
  # Aggregate interactions across bootstrap samples
  n.bootstrap <- length(store.out)
  store <- do.call(rbind, store.out)
  
  if (length(store) >= 1){
    int.tbl <- sort(table(store$Interaction), decreasing = TRUE)
    int.tbl <- int.tbl / n.bootstrap
  } else {
    return(list(interaction=numeric(0), prevalence=numeric(0)))
  }
  
  if (!is.null(varnames)) {
    names.int <- strsplit(names(int.tbl), split='_')
    names.int <- sapply(names.int, function(n) 
      paste(varnames[as.numeric(n)], collapse='_'))
    names(int.tbl) <- names.int
    
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
  # Evaluate optimal iteration based on prediction accuracy on OOB samples
  predicted <- lapply(rf.list, function(z) z$predicted)
  predicted <- lapply(predicted, function(z) as.numeric(z) - is.factor(z))
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  mse <- function(y, py) mean((py - y) ^ 2, na.rm=TRUE)
  error <- sapply(predicted, mse, y=y)
  return(which.min(error))
} 
