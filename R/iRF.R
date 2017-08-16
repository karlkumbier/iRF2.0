# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1/ncol(x), ncol(x)), 
                keep.impvar.quantile=NULL, 
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                cutoff.unimp.feature=0,  
                rit.param=list(depth=5, ntree=100, nchild=2, class.id=1, class.cut=NULL), 
                varnames.grp=NULL, 
                n.bootstrap=30,
                bootstrap.forest=TRUE,
                escv.select=FALSE,
                verbose=TRUE,
                keep.subset.var=NULL,
                return.node.feature=FALSE,
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
  if (n.core > 1) registerDoMC(n.core)  
  
  rf.list <- list()
  if (!is.null(interactions.return) | escv.select) {
    interact.list <- list()
    stability.score <- list()
    prevalence <- list()
  }
  
  weight.mat <- matrix(0, nrow=p, ncol=(n.iter + 1))
  weight.mat[,1] <- mtry.select.prob
  
  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  for (iter in 1:n.iter) {
    
    ## 1: Grow Random Forest on full data
    print(paste('iteration = ', iter))
    rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dopar% {
                                 randomForest(x, y, xtest, ytest, 
                                              ntree=ntree.id[i], 
                                              mtry.select.prob=weight.mat[,iter], 
                                              keep.forest=TRUE,
                                              track.nodes=TRUE, 
                                              ...)
                               }
    
    ## 3: update mtry.select.prob 
    if (!class.irf) 
      weight.mat[,iter + 1] <- rf.list[[iter]]$importance[,'IncNodePurity']
    else
      weight.mat[, iter + 1] <- rf.list[[iter]]$importance[,'MeanDecreaseGini']
    
    
    if (!is.null(xtest) & class.irf){
      auroc <- auc(roc(rf.list[[iter]]$test$votes[,2], ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    } else if (!is.null(xtest)) {
      pct.var <- 1 - mean((rf.list[[iter]]$test$predicted - ytest) ^ 2) / var(ytest)
      print(paste('% var explained:', pct.var * 100))
    }
  }
  
  # If escv.select = TRUE, determine optimal iteration #
  if (escv.select) {
    opt.k <- escv(rf.list, x=x, y=y)
    print(opt.k)
    interactions.return <- opt.k
  }
  
  
  for (iter in 1:n.iter) {
    ## 2.1: Find interactions across bootstrap replicates
    if (iter %in% interactions.return){
      if (verbose){cat('finding interactions ... ')}
      
      interact.list.b <- list()      
      for (i.b in 1:n.bootstrap) { 
        
        if (class.irf) {
          n.class <- table(y)
          sample.id <- mapply(function(cc, nn) sampleClass(y, cc, nn),
                              as.factor(names(n.class)), n.class)
          sample.id <- unlist(sample.id)
        } else {
          sample.id <- sample(n, n, replace=TRUE)
        }
        
        if (bootstrap.forest) { 
          #2.1.1: fit random forest on bootstrap sample
          rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                          .multicombine=TRUE, .packages='iRF') %dopar% {
                            randomForest(x[sample.id,], y[sample.id], xtest, ytest, 
                                         ntree=ntree.id[i], 
                                         mtry.select.prob=weight.mat[,iter], 
                                         keep.forest=TRUE, 
                                         track.nodes=TRUE, 
                                         ...)
                          }
        } else {
          rf.b <- rf.list[[iter]]
        }
        
        #2.1.2: run generalized RIT on rf.b to learn interactions
        ints <- generalizedRIT(rf=rf.b, 
                               x=x[sample.id,], y=y[sample.id],  
                               wt.pred.accuracy=wt.pred.accuracy,
                               class.irf=class.irf, 
                               varnames.grp=varnames.grp,
                               cutoff.unimp.feature=cutoff.unimp.feature,
                               rit.param=rit.param,
                               n.core=n.core)

        interact.list.b[[i.b]] <- ints
        rm(rf.b)       
      }
      
      # 2.2: calculate stability scores of interactions
      if (!is.null(varnames.grp))
        varnames.new <- unique(varnames.grp)
      else if (!is.null(colnames(x)))
        varnames.new <- colnames(x)
      else
        varnames.new <- 1:ncol(x)
      
      interact.list[[iter]] <- interact.list.b
      pp <- length(varnames.new)
      summary.interact <- summarizeInteract(interact.list[[iter]], 
                                             varnames=varnames.new,
                                             p=pp)
      stability.score[[iter]] <- summary.interact
      
    } # end if (find_interaction)
    
    
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) {
    out$interaction <- stability.score
  }
  
  if (escv.select) {
    out$rf.list <- out$rf.list[[opt.k]]
    out$interaction <- out$interaction[[opt.k]]
    out$opt.k <- opt.k
    out$weights <- weight.mat[,opt.k]
  }
  
  if (return.node.feature) {
    rdf <- readForest(rf.list[[n.iter]], x=x, y=y,
                     return.node.feature=TRUE,
                     wt.pred.accuracy=wt.pred.accuracy,
                     varnames.grp=varnames.grp,
                     n.core=n.core)
    out$nf <- rdf$node.feature
  }
  
  return(out)
}

generalizedRIT <- function(rf, x, y, wt.pred.accuracy, class.irf, varnames.grp,
                           cutoff.unimp.feature, rit.param, n.core) {
  
  # Extract decision paths from rf as sparse binary matrix to be passed to RIT
  rforest <- readForest(rf, x=x, y=y,
                        return.node.feature=TRUE,
                        wt.pred.accuracy=wt.pred.accuracy,             
                        n.core=n.core,
                        varnames.grp=varnames.grp)                                                        
  class.id <- rit.param$class.id

  # Select class specific leaf nodes
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  if (class.irf) { 
    select.leaf.id1 <- rforest$tree.info$prediction == 2
    select.leaf.id0 <- rforest$tree.info$prediction == 1
  } else if (is.null(rit.param$class.cut)) {
    select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  } else {
    select.leaf.id1 <- rforest$tree.info$prediction > rit.param$class.cut
    select.leaf.id0 <- rforest$tree.info$prediction <= rit.param$class.cut
  }       
  
  
  rforest1 <- subsetReadForest(rforest, select.leaf.id1)
  rforest0 <- subsetReadForest(rforest, select.leaf.id0)
  if (wt.pred.accuracy) {
    wt1 <- rforest1$tree.info$size.node * rforest1$tree.info$dec.purity
    wt0 <- rforest0$tree.info$size.node * rforest0$tree.info$dec.purity
  } else {
    wt1 <- rforest1$tree.info$size.node
    wt0 <- rforest0$tree.info$size.node
  }         
  #rm(rforest)
  
  # sample w/ replacement...
  nf1 <- rforest1$node.feature
  nf0 <- rforest0$node.feature
  nf1 <- nf1[sample(nrow(nf1), prob=wt1, replace=TRUE),]
  nf0 <- nf0[sample(nrow(nf0), prob=wt0, replace=TRUE),]
  
  if (sum(select.leaf.id) < 2){
    return(character(0))
  } else {  
    # group features if specified
    if (!is.null(varnames.grp)) nf1 <- groupFeature(nf1, grp=varnames.grp)
    if (!is.null(varnames.grp)) nf0 <- groupFeature(nf0, grp=varnames.grp)
    
    # drop feature if cutoff.unimp.feature is specified
    if (cutoff.unimp.feature > 0){
      if (!class.irf)   
        rfimp <- rf$importance[,'IncNodePurity']
      else            
        rfimp <- rf$importance[,'MeanDecreaseGini']
      drop.id <- which(rfimp < quantile(rfimp, prob=cutoff.unimp.feature))
      nf1[,drop.id] <- FALSE 
      nf0[,drop.id] <- FALSE
    }                     
    
    if (rit.param$class.id == 1) {
      z <- nf1 != 0
      z0 <- nf0 != 0
    } else {
      z <- nf0 != 0
      z0 <- nf1 != 0
    }
    
    interactions <- RIT(z=z, z0=z0, depth=rit.param$depth, n_trees=rit.param$ntree, 
                         branch=rit.param$nchild, n_cores=n.core) 
    
    interactions <- interactions$Class1
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
  p <- ncol(node.feature)
  pp <- p / 2
  sparse.mat <- is(node.feature, 'Matrix')
  grp.names <- unique(grp)
  # TODO: is max correct for hyperrectangle? Should depend on direction
  makeGroup <- function(x, g) apply(as.matrix(x[,grp == g]), MAR=1, max) 
  node.feature.new0 <- sapply(grp.names, makeGroup, x=node.feature[,1:pp])
  node.feature.new1 <- sapply(grp.names, makeGroup, x=node.feature[,(pp + 1):p])
  node.feature.new <- cbind(node.feature.new0, node.feature.new1)
  
  if (sparse.mat) node.feature.new <- Matrix(node.feature.new, sparse=TRUE)
  
  colnames(node.feature.new) <- c(grp.names, grp.names)
  return(node.feature.new)
}



summarizeInteract <- function(store.out, varnames=NULL, p){
  # Aggregate interactions across bootstrap samples
  n.bootstrap <- length(store.out)
  store <- do.call(rbind, store.out)
  
  if (length(store) >= 1){
    int.tbl <- sort(table(store$Interaction), decreasing = TRUE)
    int.tbl <- int.tbl / n.bootstrap
    
    prev.tbl <- c(by(store$Prevalence, store$Interaction, sum))
    prev.tbl <- prev.tbl / n.bootstrap
    prev.tbl <- prev.tbl[names(int.tbl)]
  } else {
    return(list(interaction=numeric(0), prevalence=numeric(0)))
  }
  
  if (!is.null(varnames)) {
    stopifnot (names(int.tbl) == names(prev.tbl))
    names.int <- lapply(names(int.tbl), strsplit, split='_')
    names.int <- lapply(names.int, unlist)
    names.int <- sapply(names.int, function(n) {
      nn <- as.numeric(n)
      direction <- as.numeric(nn > p)
      direction <- c('-', '+')[direction + 1]
      nn <- nn %% p
      nn[nn == 0] <- p
      nn.direction <- paste0(varnames[nn], direction)
      return(paste(nn.direction, collapse='_'))
    })
    names(int.tbl) <- names.int
    names(prev.tbl) <- names.int
    
  }
  out <- list(interaction=int.tbl)
  return(out)
}

sampleClass <- function(y, cl, n) {
  # Sample indices specific to a given class
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

estStabilityIter <- function(rfobj, x, y) {
  # Evaluate estimation stability for an RF object
  y.hat <- predict(rfobj, newdata=x, predict.all=TRUE)
  y.bar <- y.hat$aggregate
  y.agg <- apply(y.hat$individual, MAR=2, function(z) sum((z - y.bar) ^ 2))
  es <- mean(y.agg) 
  
  error <- mean((rfobj$predicted - y) ^ 2, na.rm=TRUE)
  return(c(es=es, cv=error))
}

escv <- function(rf.list, x, y) {
  # Evaluate optimal iteration based on ESCV critereon 
  escv.list <- sapply(rf.list, estStabilityIter, x=x, y=y)
  min.cv <- which.min(escv.list[2,])
  min.es <- which.min(escv.list[1,min.cv:ncol(escv.list)])
  escv.opt <- min.cv + min.es - 1
  return(escv.opt)
} 
