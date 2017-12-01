# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1, ncol(x)),
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                rit.param=list(depth=5, ntree=500, nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL, class.qt=0.5), 
                varnames.grp=NULL, 
                n.bootstrap=20,
                select.iter=FALSE,
                verbose=TRUE,
                keep.subset.var=NULL,
                get.prevalence=FALSE,
                take.diff=FALSE,
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
  importance <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
   
  
  rf.list <- list()
  if (!is.null(interactions.return) | select.iter) stability.score <- list()
  if (get.prevalence) prev.list <- list()

  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  registerDoMC(n.core)
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
  if (select.iter) interactions.return <- selectIter(rf.list, y=y)
  
  for (iter in interactions.return) {
    
    # Find interactions across bootstrap replicates
    if (verbose) cat('finding interactions ... ')
    
    interact.list.b <- list()      
    if (get.prevalence) prev.list.b <- list()

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
      
      # Use feature weights from current iteraction of full data RF
      mtry.select.prob <- rf.list[[iter]]$importance
      
      # Fit random forest on bootstrap sample
      rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                      .multicombine=TRUE, .packages='iRF') %dopar% {
                        randomForest(x[sample.id,], y[sample.id], 
                                     xtest, ytest, 
                                     ntree=ntree.id[i], 
                                     mtry.select.prob=mtry.select.prob, 
                                     keep.forest=TRUE, 
                                     ...)
                      }
      
      # Run generalized RIT on rf.b to learn interactions
      ints <- generalizedRIT(rf=rf.b, x=x[sample.id,], y=y[sample.id],
                             wt.pred.accuracy=wt.pred.accuracy,
                             varnames.grp=varnames.grp,
                             rit.param=rit.param,
                             get.prevalence=get.prevalence,
                             take.diff=take.diff,
                             n.core=n.core)
      
      interact.list.b[[i.b]] <- ints$interactions
      if (get.prevalence) prev.list.b[[i.b]] <- ints$prevalence

      rm(rf.b)       
    }
    
    # Calculate stability scores of interactions
    stability.score[[iter]] <- summarizeInteract(interact.list.b)
    if (get.prevalence) 
      prev.list[[iter]] <- summarizePrevalence(prev.list.b, interact.list.b, n.bootstrap)

  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) out$interaction <- stability.score
  if (get.prevalence) out$prevalence <- prev.list

  if (select.iter) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    out$opt.k <- interactions.return
    out$weights <- rf.list[[interactions.return]]$importance[,importance]
    if (get.prevalence) out$prevalence <- out$prevalence[[interactions.return]]
  }
  return(out)
}

generalizedRIT <- function(rf, x, y, 
                           wt.pred.accuracy=FALSE, 
                           varnames.grp=NULL,
                           rit.param=list(depth=5, ntree=500, nchild=2, 
                                          class.id=1, min.nd=1,
                                          class.cut=NULL, class.qt=0.5), 
                           get.prevalence=FALSE,
                           take.diff=take.diff,
                           n.core=1) {
  
  class.irf <- is.factor(y)

  # Extract decision paths and tree metadata from random forest
  rforest <- readForest(rf, x=x, y=y, 
                        return.node.feature=TRUE,
                        wt.pred.accuracy=wt.pred.accuracy, 
                        varnames.grp=varnames.grp,
                        n.core=n.core)                                                        
  
  # Select class specific leaf nodes
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  
  if (class.irf) {
    # classification: subset leaf nodes by class
    class.id <- rit.param$class.id
    select.leaf.id <- rforest$tree.info$prediction == as.numeric(class.id) + 1
  } else if (!is.null(rit.param$class.cut)) {
    # regression with cutoff: sample all leaf nodes that exceed cutoff
    select.leaf.id <- rforest$tree.info$prediction > rit.param$class.cut
  } else if (!is.null(rit.param$class.qt)) {
    # regression with quantile cutoff: determine cutoff from quantile and 
    # sample leaf nodes
    qt.cut <- quantile(y, rit.param$class.qt)
    select.leaf.id <- rforest$tree.info$prediction > qt.cut
  }       
 
  rforest <- subsetReadForest(rforest, select.leaf.id)
  if (sum(select.leaf.id) < 2) {
    return(character(0))
  } else {
    out <- runRIT(subsetReadForest(rforest, select.leaf.id),
                  varnames.grp, wt.pred.accuracy, rit.param,
                  get.prevalence)
    if (take.diff) { 
      ints.diff <- runRIT(subsetReadForest(rforest, !select.leaf.id),
                            varnames.grp, wt.pred.accuracy, rit.param,
                            get.prevalence)
      id.rm <- out$interactions %in% ints.diff$interactions
      out$interactions <- out$interactions[!id.rm]
      out$prevalence <- out$prevalence[!id.rm]
    }
  }
  return(out)
}

runRIT <- function(rforest, varnames.grp, wt.pred.accuracy, 
                   rit.param, get.prevalence) {
  
  # Set weights for leaf node sampling using either size or size and accuracy
  if (wt.pred.accuracy) {
    wt <- rforest$tree.info$size.node * rforest$tree.info$dec.purity
  } else {
    wt <- rforest$tree.info$size.node
  }         
  
  # group features if specified
  if (!is.null(varnames.grp)) 
    rforest$node.feature <- groupFeature(rforest$node.feature, grp=varnames.grp)
  
  # remove nodes below specified size threshold
  id.rm <- rforest$tree.info$size.node >= rit.param$min.nd
  rforest$node.feature <- rforest$node.feature[!id.rm,]
  rforest$tree.info <- rforest$tree.info[!id.rm,]
  wt <- wt[!id.rm]
  
  interactions <- RIT(rforest$node.feature, weights=wt, depth=rit.param$depth,
                      n_trees=rit.param$ntree, branch=rit.param$nchild,
                      n_cores=n.core) 
  
  # Group interactions and rename by variable name
  if (!is.null(varnames.grp))
    varnames.new <- unique(varnames.grp)
  else if (!is.null(colnames(x)))
    varnames.new <- colnames(x)
  else
    varnames.new <- 1:ncol(x)
  
  interactions <- gsub(' ', '_', interactions$Interaction)
  if (get.prevalence) 
    prev <- sapply(interactions, prevalence, nf=rforest$node.feature, wt=wt)
  
  interactions <- nameInts(interactions, varnames.new)
  if (get.prevalence) 
    names(prev) <- interactions
  
  out <- list()
  out$interactions <- interactions
  if (get.prevalence) out$prevalence <- prev
  return(out)
}

nameInts <- function(int, varnames) {
  # Convert interactions indicated by indices to interactions indicated by name
  ints.split <- strsplit(int, '_')
  varnames.unq <- unique(varnames)
  p <- length(varnames.unq)

  ints.split <- lapply(ints.split, as.numeric)
  ints.signs <- lapply(ints.split, function(z) ifelse(z > p, '+', '-'))
  ints.split <- lapply(ints.split, '%%', p)

  ints.name <- mapply(function(i, s) paste(varnames.unq[i], s, sep=''),
                      ints.split, ints.signs, SIMPLIFY=FALSE)
  ints.name <- sapply(ints.name, paste, collapse='_')
  return(ints.name)
}

summarizePrevalence <- function(prev, int, n.bootstrap) {
  # summarize interaction prevalence across bootstrap samples
  prev <- unlist(prev)
  nn <- names(prev)
  prev.summary <- by(prev, nn, prevalenceSummary)
  prev.summary <- do.call(rbind, c(prev.summary))
  prev.summary <- data.frame(prev.summary)
  prev.summary$int <- rownames(prev.summary)
  prev.summary <- arrange(prev.summary, desc(prop), desc(mean))
  prev.summary$prop <- prev.summary$prop / n.bootstrap
  return(prev.summary)
}

prevalenceSummary <- function(x) {
  # determine the min, mean, and maximum prevalence across 
  # bootstrap samples
  x.min <- min(x)
  x.mean <- mean(x)
  x.max <- max(x)
  x.n <- length(x)
  return(c(min=x.min, mean=x.mean, max=x.max, prop=x.n))
}


prevalence <- function(int, nf, wt=rep(1, ncol(nf))) {
  # calculate the decision path prevalence of a single interaction
  int <- as.numeric(strsplit(int, '_')[[1]])
  int.nd <- apply(nf[,int], MAR=1, function(z) all(z == 1))
  prev <- sum(wt[int.nd]) / sum(wt)
  return(prev)
}


subsetReadForest <- function(rforest, subset.idcs) { 
  # Subset nodes from readforest output 
  if (!is.null(rforest$node.feature)) 
    rforest$node.feature <- rforest$node.feature[subset.idcs,]
  
  if (!is.null(rforest$tree.info))
    rforest$tree.info <- rforest$tree.info[subset.idcs,]
  
  if (!is.null(rforest$node.obs))
    rforest$node.obs <- rforest$node.obs[subset.idcs,]
  
  return(rforest)
}

groupFeature <- function(node.feature, grp){
  # Group replicated features in node.feature 
  sparse.mat <- is(node.feature, 'Matrix')
  
  grp.names <- unique(grp)
  makeGroup <- function(x, g) apply(as.matrix(x[,grp == g]), MAR=1, max) 
  node.feature.new <- sapply(grp.names, makeGroup, x=node.feature)
  
  if (sparse.mat) node.feature.new <- Matrix(node.feature.new, sparse=TRUE)
  colnames(node.feature.new) <- grp.names
  
  return(node.feature.new)
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
  # Evaluate optimal iteration based on ESCV critereon 
  predicted <- lapply(rf.list, function(z) as.numeric(z$predicted) - is.factor(y))
  
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  mse <- function(y, py) mean((py - y) ^ 2, na.rm=TRUE)
  error <- sapply(predicted, mse, y=y)
  return(which.min(error))
} 
