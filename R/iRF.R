# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=matrix(1, nrow=nrow(x), ncol=ncol(x)),
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                rit.param=list(depth=5, ntree=500, nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL, class.qt=0.5), 
                varnames.grp=NULL, 
                n.bootstrap=20,
                bootstrap.forest=TRUE,
                select.iter=FALSE,
                verbose=TRUE,
                keep.subset.var=NULL,
                local=FALSE,
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
  importance.feature <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
  registerDoMC(n.core)  
  
  rf.list <- list()
  if (!is.null(interactions.return) | select.iter) {
    stability.score <- list()
  }
  
  if (local) local.list <- list()

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
                                              ...)
                               }
    
    ## update feature selection probabilities
    if (local){
      mtry.select.prob <- rf.list[[iter]]$obsgini
    } else {
      mtry.select.prob <- matrix(colSums(rf.list[[iter]]$obsgini),
                                 nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
    }

    if (!is.null(xtest) & class.irf & verbose){
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
        
        # use feature weights from current iteraction of full data RF
        if (local) {
          mtry.select.prob <- rf.list[[iter]]$obsgini[sample.id,]
        } else {
          mtry.select.prob <- matrix(colSums(rf.list[[iter]]$obsgini[sample.id,]), 
                                   nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
        }
       

        # fit random forest on bootstrap sample
        rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                        .multicombine=TRUE, .packages='iRF') %dopar% {
                          randomForest(x[sample.id,], y[sample.id], 
                                       xtest, ytest, 
                                       ntree=ntree.id[i], 
                                       mtry.select.prob=mtry.select.prob, 
                                       keep.forest=TRUE, 
                          ...)
                        }
      } else {
        rf.b <- rf.list[[iter]]
      }
      
      # run generalized RIT on rf.b to learn interactions
      ints <- generalizedRIT(rf=rf.b, x=x[sample.id,], y=y[sample.id],
                             wt.pred.accuracy=wt.pred.accuracy,
                             varnames.grp=varnames.grp,
                             rit.param=rit.param,
                             local=FALSE, # NOT USING LOCAL RIT
                             n.core=n.core)
      print('backsampling')
      if (local) ints <- lapply(ints, function(z) unique(sample.id[z]))
      interact.list.b[[i.b]] <- ints
      rm(rf.b)       
    }
    
    # calculate stability scores of interactions
    stability.score[[iter]] <- summarizeInteract(interact.list.b, local=local)
    if (FALSE) { # NOT USING LOCAL RIT
      print('grouping')
      xx <- unlist(interact.list.b, recursive=FALSE)
      local.list[[iter]] <- lapply(unique(names(xx)), function(i)
        unique(unlist(xx[names(xx) == i])))      
      names(local.list[[iter]]) <- unique(names(xx))
    }
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) out$interaction <- stability.score
  #if (local) out$local <- local.list NOT USING LCOAL RIT

  if (select.iter) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    out$opt.k <- interactions.return
    out$weights <- rf.list[[interactions.return]]$importance[,importance.feature]
    #if (local) out$local <- local.list[[interactions.return]]
  }
  
  return(out)
}

generalizedRIT <- function(rf, x, y, 
                           wt.pred.accuracy=FALSE, 
                           varnames.grp=NULL,
                           rit.param=list(depth=5, ntree=500, nchild=2, 
                                          class.id=1, min.nd=10,
                                          class.cut=NULL, class.qt=0.5), 
                           local=FALSE, 
                           n.core=1) {
  
  class.irf <- is.factor(y)

  # Extract decision paths and tree metadata from random forest
  rforest <- readForest(rf, x=x, y=y, 
                        return.node.feature=TRUE,
                        return.node.obs=local,
                        wt.pred.accuracy=wt.pred.accuracy, 
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
  nf <- rforest$node.feature

  # Set weights for leaf node sampling using either size or size and accuracy
  if (wt.pred.accuracy) {
    wt <- rforest$tree.info$size.node * rforest$tree.info$dec.purity
  } else {
    wt <- rforest$tree.info$size.node
  }         
  
  if (sum(select.leaf.id) < 2) {
    return(character(0))
  } else {  
    
    # group features if specified
    if (!is.null(varnames.grp)) nf <- groupFeature(nf, grp=varnames.grp)
    p <- ncol(nf)
    
    # add observation data for each node
    if (local) nf <- cbind(nf, rforest$node.obs)
      
    id <- rforest$tree.info$size.node >= rit.param$min.nd
    nf <- nf[id,]
    rforest$tree.info <- rforest$tree.info[id,]
    wt <- wt[id]
    
    interactions <- RIT(nf, weights=wt, depth=rit.param$depth,
                        n_trees=rit.param$ntree, branch=rit.param$nchild,
                        n_cores=n.core) 
    
    # Group interactions and rename by variable name
    if (!is.null(varnames.grp))
      varnames.new <- unique(varnames.grp)
    else if (!is.null(colnames(x)))
      varnames.new <- colnames(x)
    else
      varnames.new <- 1:ncol(x)
    
    if (local) {
      interactions <- groupLocalInteracts(interactions$Interaction, n, p)
      names(interactions) <- nameInts(names(interactions), varnames.new)
    } else {
      interactions <- gsub(' ', '_', interactions$Interaction)
      interactions <- nameInts(interactions, varnames.new)
    }
    
    return(interactions)
  }                   
}

nameInts <- function(int, varnames) {
  
  ints.split <- strsplit(int, '_')
  varnames.unq <- unique(varnames)
  ints.name <- lapply(ints.split, function(i) varnames.unq[as.numeric(i)])
  ints.name <- sapply(ints.name, paste, collapse='_')
  return(ints.name)
}


groupLocalInteracts <- function(interactions, n, p) {
  # Aggregate local interactions
  
  ints.split <- strsplit(interactions, ' ')
  ints.split <- sapply(ints.split, as.numeric)
  
  # split interactions into feature and observation
  ints.feat <- sapply(ints.split, function(z) paste(z[z <= p], collapse='_'))
  ints.obs <- lapply(ints.split, function(z) z[z > p])
  
  # remove interactions between only observations
  id.rm <- ints.feat == ''
  ints.obs <- ints.obs[!id.rm]
  ints.feat <- ints.feat[!id.rm]
  
  # aggreate observations by interaction
  agg <- lapply(unique(ints.feat), function(i) 
    unique(unlist(ints.obs[ints.feat == i])))
  id.rm <- sapply(agg, function(z) length(z) == 0)
  agg <- agg[!id.rm]
  names(agg) <- unique(ints.feat)[!id.rm]
  agg <- lapply(agg, '-', p)
  
  return(agg)
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
  # Group feature level data in node.feature 
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
  print('summarizing')
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
