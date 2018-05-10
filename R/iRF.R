# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1, ncol(x)),
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                rit.param=list(depth=5, ntree=500, 
                               nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL), 
                varnames.grp=NULL, 
                n.bootstrap=1,
                select.iter=FALSE,
                get.prevalence=TRUE,
                int.direction=TRUE,
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
        #save(file=paste0(bootstrap.path, out.file), rf.b)
      } else {
        out.file <- NULL
      }
      
     
      # Run generalized RIT on rf.b to learn interactions
      ints <- generalizedRIT(rf=rf.b, x=x, y=y,
                             wt.pred.accuracy=wt.pred.accuracy,
                             varnames.grp=varnames.grp,
                             rit.param=rit.param,
                             get.prevalence=get.prevalence,
                             int.direction=int.direction,
                             out.file=out.file,
                             n.core=n.core)
      
      interact.list[[i.b]] <- ints$int
      if (get.prevalence) prev.list[[i.b]] <- ints$prev

      rm(rf.b)       
    }
    
    # Calculate stability scores of interactions
    stability.score[[iter]] <- summarizeInteract(interact.list)
    
    
    if (get.prevalence) 
      prev.list[[iter]] <- summarizePrev(prev.list)
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) out$interaction <- stability.score
  if (get.prevalence) out$prevalence <- prev.list

  if (select.iter) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    if (get.prevalence) out$prevalence <- out$prevalence[[interactions.return]]
  }
  return(out)
}

generalizedRIT <- function(rf=NULL, x=NULL, y=NULL, rforest=NULL, 
                           wt.pred.accuracy=FALSE, 
                           varnames.grp=NULL,
                           rit.param=list(depth=5, ntree=500, 
                                          nchild=2, class.id=1, 
                                          min.nd=1, class.cut=NULL), 
                           get.prevalence=FALSE,
                           int.direction=FALSE,
                           out.file=NULL,
                           n.core=1) {
  
  out <- list()
  stopifnot(!is.null(rf) & !is.null(x)| !is.null(rforest))
  if (is.null(varnames.grp) & is.null(colnames(x))) {
    varnames.grp <- as.character(1:ncol(x))
    p <- ncol(x)
  } else if (is.null(varnames.grp)) {
    varnames.grp <- colnames(x)
    p <- ncol(x)
  } else {
    p <- length(unique(varnames.grp))
  }
  
  varnames.unq <- unique(varnames.grp)

  # Extract decision paths and tree metadata from random forest
  class.irf <- is.factor(y)
  if (is.null(rforest)) {
    rforest <- readForest(rf, x=x, y=y, 
                          return.node.feature=TRUE,
                          wt.pred.accuracy=wt.pred.accuracy, 
                          varnames.grp=varnames.grp,
                          get.split=TRUE,
                          n.core=n.core)
  }

  # Collapse node feature matrix if not tracking split directions
  if (!int.direction) {
    rforest$node.feature <- rforest$node.feature[,1:p] + 
      rforest$node.feature[,(p + 1):(2 * p)]
  }  
  
  if (!is.null(out.file)) save(file=out.file, rforest)

  # Select class specific leaf nodes
  if (class.irf) 
    select.id <- rforest$tree.info$prediction == rit.param$class.id + 1
  else
    select.id <- rforest$tree.info$prediction > rit.param$class.cut

  # Run RIT on selected and non-selected nodes to determine interactions for
  # each group
  if (sum(select.id) < 2) {
    return(character(0))
  } else {
    out$int <- runRIT(subsetReadForest(rforest, select.id),
                      wt.pred.accuracy, rit.param, n.core)
    int.names <- nameInts(out$int, varnames.unq)
    
    if (get.prevalence) { 
      wt <- rforest$tree.info$size.node
      out$prev$i1 <- unlist(mclapply(out$int, prevalence, 
                            nf=rforest$node.feature[select.id,], 
                            wt=wt[select.id], mc.cores=n.core))
      out$prev$i0 <- unlist(mclapply(out$int, prevalence,
                            nf=rforest$node.feature[!select.id,],
                            wt=wt[!select.id], mc.cores=n.core))
      
      names(out$prev$i1) <- int.names
      names(out$prev$i0) <- int.names 
    }
    
    out$int <- int.names
  }
  return(out)
}

runRIT <- function(rforest, wt.pred.accuracy, rit.param, n.core=1) {
 
  # Set weights for leaf node sampling using either size or size and accuracy
  if (wt.pred.accuracy) 
    wt <- rforest$tree.info$size.node * rforest$tree.info$dec.purity
  else 
    wt <- rforest$tree.info$size.node
           
  # remove nodes below specified size threshold
  id.rm <- rforest$tree.info$size.node < rit.param$min.nd
  if (mean(id.rm) == 1) {
    warning(paste('No nodes with greater than ', rit.param$min.nd,
                  'observations. Using all nodes'))
    id.rm <- rep(FALSE, length(id.rm))
  }
  rforest <- subsetReadForest(rforest, !id.rm)
  wt <- wt[!id.rm]
  

  interactions <- RIT(rforest$node.feature, weights=wt, depth=rit.param$depth,
                      n_trees=rit.param$ntree, branch=rit.param$nchild,
                      n_cores=n.core) 
  
  # Rename interactions using variable names and evaluate prevalence
  interactions <- strsplit(interactions$Interaction, ' ')
  interactions <- lapply(interactions, intSubsets)
  interactions <- unique(unlist(interactions))
  return(interactions)
}

intSubsets <- function(int) {
  # Generate all lower order subsets of specified interaction
  int.order <- length(int)
  int.subs <- lapply(1:int.order, combn, x=int, simplify=FALSE)
  int.subs <- unlist(int.subs, recursive=FALSE)
  int.subs <- sapply(int.subs, function(z) paste(sort(z), collapse='_'))
  return(int.subs)
}


nameInts <- function(int, varnames, directed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by name
  ints.split <- strsplit(int, '_')
  varnames.unq <- unique(varnames)
  p <- length(varnames.unq)

  ints.split <- lapply(ints.split, as.numeric)
  if (directed) {
    ints.signs <- lapply(ints.split, function(z) ifelse(z > p, '+', '-'))
  } else {
    ints.signs <- ''
  }

  ints.split <- lapply(ints.split, function(z) z %% p + p * (z == p))
  ints.name <- mapply(function(i, s) paste(varnames.unq[i], s, sep=''),
                      ints.split, ints.signs, SIMPLIFY=FALSE)
  ints.name <- sapply(ints.name, function(z) paste(sort(z), collapse='_'))
  return(ints.name)
}

summarizePrev <- function(prev) {
  # summarize interaction prevalence across bootstrap samples
  n.bs <- length(prev)
  prev1 <- unlist(lapply(prev, function(z) z$i1))
  prev0 <- unlist(lapply(prev, function(z) z$i0))
  
  prev <- data.frame(int=names(prev1), prev1=prev1, prev0=prev0)
  prev <- group_by(prev, int) %>%
    summarize(prev1=mean(prev1), prev0=mean(prev0), n=n()/n.bs) %>%
    mutate(diff=(prev1-prev0)) %>%
    arrange(desc(diff))
  return(prev)
}

prevalence <- function(int, nf, wt=rep(1, ncol(nf))) {
  # calculate the decision path prevalence of a single interaction
  int <- as.numeric(strsplit(int, '_')[[1]])
  id.rm <- wt == 0
  wt <- wt[!id.rm]

  if (length(int) == 1)
    int.id <- nf[!id.rm, int] != 0
  else
    int.id <- apply(nf[!id.rm, int], MAR=1, function(z) all(z != 0))
  
  prev <- sum(wt[int.id]) / sum(wt)
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
  print(error)
  return(which.min(error))
} 
