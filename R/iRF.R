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
  
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  
  # Check whether iRF can be run for given inputs
  if (!is.matrix(x) | (!is.null(xtest) & !is.matrix(xtest)))
    stop('either x or xtest is not a matrix !')
  if (!is.numeric(x) | (!is.null(xtest) & !is.numeric(xtest)))
    stop('either x or xtest is not a numeric matrix!')
  if (ncol(x) < 2 & !is.null(interactions.return))
    stop('cannot find interaction - X has less than two columns!')
  if (any(interactions.return > n.iter))
    stop('interactions to return greater than niter')
  
  # initialize outputs
  rf.list <- list()
  if (!is.null(interactions.return)) {
    interact.list <- list()
    stability.score <- list()
  }
  
  # set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  for (iter in 1:n.iter){
    ## 1: Grow Random Forest on full data using all available cores
    
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
    
    ## 2: Find interactions that are consistently recovered over bootstrap 
    ## replicates
    if (iter %in% interactions.return){
      if (verbose){cat('finding interactions ... ')}
      
      stability.score[[iter]] <- list()      
      interact.list[[iter]] <- mclapply(1:n.bootstrap, function(i.b) {
        
        if (class.irf) {
          n.class <- table(y)
          sample.id <- mapply(function(cc, nn) 
            sample(which(y == cc), nn, replace=TRUE),
            as.factor(names(n.class)), n.class)
          sample.id <- unlist(sample.id)
        } else {
          sample.id <- sample(n, n, replace=TRUE)
        }
        
        #2.1.1: fit random forest on bootstrap sample
        rf.b <- randomForest(x[sample.id,], y[sample.id], 
                             xtest, ytest, 
                             ntree=ntree, 
                             mtry_select_prob=mtry.select.prob, 
                             keep.forest=TRUE, 
                             track.nodes=TRUE, 
                             ...)
        
        ints <- generalizedRIT(rf=rf.b, x=x[sample.id,], y=y[sample.id], 
                               subsetFun=node.sample$subset, 
                               wtFun=node.sample$wt,
                               class.irf=class.irf, 
                               class.id=class.id, 
                               varnames.grp=varnames.grp, 
                               rit.param=rit.param,
                               cutoff.unimp.feature=cutoff.unimp.feature) 
        return(ints)       
        
      }, mc.cores=n.core)
      
      # 2.2: calculate stability scores of interactions
      if (!is.null(varnames.grp))
        varnames.new <- unique(varnames.grp)
      else if (!is.null(colnames(x)))
        varnames.new <- colnames(x)
      else
        varnames.new <- 1:ncol(x)
      
      stability.score[[iter]] <- summarize_interact(interact.list[[iter]], 
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
                           rit.param) {
  
  # Extract decision paths from rf as binary matrix to be passed to RIT
  rforest <- readForest(rf, x=x, 
                        return_node_feature=TRUE,
                        subsetFun=subsetFun, 
                        wtFun=wtFun)
  
  # Select class specific leaf nodes if classification
  select.leaf.id <- rep(TRUE, nrow(rforest$tree_info))
  if (class.irf & all(y %in% c(0, 1))) {
    select.leaf.id <- rforest$tree_info$prediction == as.numeric(class.id) + 1
    rforest <- subsetReadForest(rforest, select.leaf.id)
  }
  nf <- rforest$node_feature
  
  
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
      
      drop_id <- which(rfimp < quantile(rfimp, prob=cutoff.unimp.feature))
      nf[,drop_id] <- FALSE
    }
    
    interactions <- RIT(nf, depth=rit.param[1], n_trees=rit.param[2], 
                        branch=rit.param[3])
    interactions <- interactions$Interaction
    interactions <- gsub(' ', '_', interactions)
    return(interactions)
  }
}

subsetReadForest <- function(rforest, subset.idcs) {
  # subset nodes from readforest output 
  if (!is.null(rforest$node_feature)) 
    rforest$node_feature <- rforest$node_feature[subset.idcs,]
  if (!is.null(rforest$node_data)) 
    rforest$node_data <- t(rforest$node_data[,subset.idcs])
  if(!is.null(rforest$tree_info))
    rforest$tree_info <- rforest$tree_info[subset.idcs,]
  return(rforest)
}
