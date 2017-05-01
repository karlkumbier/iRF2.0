# iteratively grows random forests, finds case specific feature interactions
# in binary classification problems
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1/ncol(x), ncol(x)), 
                keep.impvar.quantile=NULL, 
                interactions.return=NULL, 
                node.sample=list(subset=function(x) rep(TRUE, nrow(x)),
                                wt=function(x) x$size_node), 
                cutoff.unimp.feature=0, 
                class.id=1, 
                rit.param=c(5, 100, 2), 
                varnames.grp=NULL, 
                n.bootstrap=30, 
                verbose=TRUE, 
                ...) {

  require(data.table)  
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  
  if (!is.matrix(x) | (!is.null(xtest) & !is.matrix(xtest)))
    stop('either x or xtest is not a matrix !')
  if (!is.numeric(x) | (!is.null(xtest) & !is.numeric(xtest)))
    stop('either x or xtest is not a numeric matrix!')
  if (ncol(x) < 2 & !is.null(interactions.return))
    stop('cannot find interaction - X has less than two columns!')
  if (any(interactions.return > n.iter))
    stop('interaction iteration to return greater than n.iter')
  
  # initialize outputs
  rf.list <- list()
  if (!is.null(interactions.return)) {
    interact.list <- list()
    stability.score <- list()
  }
  
  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  for (iter in 1:n.iter){
    
    ## 1: Grow Random Forest on full data
    print(paste('iteration = ', iter))
    rf.list[[iter]] <- foreach(nt=ntree.id, .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dopar% {
                                 randomForest(x, y, 
                                              xtest, ytest, 
                                              ntree=nt, 
                                              mtry_select_prob=mtry.select.prob, 
                                              keep.forest=TRUE, 
                                              ...)
                               }
    
    ## 2.1: Find interactions across bootstrap replicates
    if (iter %in% interactions.return){
      if (verbose){cat('finding interactions ... ')}
      
      stability.score[[iter]] <- list()      
      interact.list[[iter]] <- lapply(1:n.bootstrap, function(i.b) {
        
        if (class.irf) {
          n.class <- table(y)
          sampleCl <- function(y, cc, n) sample(which(y == cc), n, replace=TRUE)
          sample.id <- mapply(function(cc, nn) sampleCl(y, cc, nn),
            as.factor(names(n.class)), n.class)
          sample.id <- unlist(sample.id)
        } else {
          sample.id <- sample(n, n, replace=TRUE)
        }
        
        #2.1.1: fit random forest on bootstrap sample
        rf.b <- foreach(nt=ntree.id, .combine=combine, 
                        .multicombine=TRUE, .packages='iRF') %dopar% {
                          randomForest(x[sample.id,], y[sample.id], 
                                       xtest, ytest, 
                                       ntree=nt, 
                                       mtry_select_prob=mtry.select.prob, 
                                       keep.forest=TRUE, 
                                       track.nodes=TRUE,
                                       ...)
                        }
          
        
        #2.1.2: run generalized RIT on rf.b to learn interactions
        ints <- generalizedRIT(rf=rf.b, x=x[sample.id,], y=y[sample.id], 
                               subsetFun=node.sample$subset, 
                               wtFun=node.sample$wt,
                               class.irf=class.irf, 
                               class.id=class.id, 
                               varnames.grp=varnames.grp,
                               cutoff.unimp.feature=cutoff.unimp.feature,
                               rit.param=rit.param,
                               n.core=n.core) 
        return(ints)       
        
      })
      
      # 2.2: calculate stability scores of interactions
      if (!is.null(varnames.grp))
        varnames.new <- unique(varnames.grp)
      else if (!is.null(colnames(x)))
        varnames.new <- colnames(x)
      else
        varnames.new <- 1:ncol(x)
      
      stability.score[[iter]] <- summarizeInteract(interact.list[[iter]], 
                                                    varnames=varnames.new)      
      
    } # end if (find_interaction)
    
    ## 3: update mtry.select.prob 
    if (!class.irf) 
      mtry.select.prob <- rf.list[[iter]]$importance[,'IncNodePurity']
    else
      mtry.select.prob <- rf.list[[iter]]$importance[,'MeanDecreaseGini']
    
    
    if (!is.null(xtest) & class.irf){
      auroc <- auc(roc(rf.list[[iter]]$test$votes[,2], ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    }
    
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)){
    out$interaction <- stability.score
  }
  return(out)
}


generalizedRIT <- function(rf, x, y,
                           subsetFun, 
                           wtFun, 
                           class.irf, 
                           class.id, 
                           varnames.grp,
                           cutoff.unimp.feature, 
                           rit.param,
                           n.core) {
  
  # Extract decision paths from rf as binary matrix to be passed to RIT
  rforest <- readForest(rf, x=x, 
                        return.node.feature=TRUE,
                        subsetFun=subsetFun, 
                        wtFun=wtFun,
                        n.core=n.core)
  
  # Select class specific leaf nodes if classification
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  if (class.irf & all(y %in% c(0, 1))) {
    select.leaf.id <- rforest$tree.info$prediction == as.numeric(class.id) + 1
    rforest <- subsetReadForest(rforest, select.leaf.id)
  }
  nf <- rforest$node.feature
  
  
  if (sum(select.leaf.id) < 2){
    return(character(0))
  } else {
    # group features if specified
    if (!is.null(varnames.grp)) nf <- groupFeature(nf, grp = varnames.grp)
    
    # drop feature if cutoff.unimp.feature is specified
    if (cutoff.unimp.feature > 0){
      if (!class.irf)
        rfimp <- rf$importance[,'IncNodePurity']
      else
        rfimp <- rf$importance[,'MeanDecreaseGini']
      
      drop.id <- which(rfimp < quantile(rfimp, prob=cutoff.unimp.feature))
      nf[,drop.id] <- FALSE
    }
    
    interactions <- RIT(nf, depth=rit.param[1], n_trees=rit.param[2], 
                        branch=rit.param[3], n_cores=n.core)
    interactions <- interactions$Interaction
    interactions <- gsub(' ', '_', interactions)
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
  store <- unlist(store.out)
  
  if (length(store) >= 1){
    store.tbl <- sort(table(store), decreasing = TRUE)
    store.tbl <- store.tbl / n.bootstrap
  } else {
    return(character(0))
  }
  
  if (!is.null(varnames)) {
    names.int <- lapply(names(store.tbl), strsplit, split='_')
    names.int <- lapply(names.int, unlist)
    names.int <- sapply(names.int, function(n) {
      nn <- as.numeric(n)
      return(paste(varnames[nn], collapse='_'))
    })
    names(store.tbl) <- names.int
  }
  
  return(store.tbl)
}
