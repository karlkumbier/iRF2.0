# iteratively grows random forests, finds case specific feature interactions
# in binary classification problems
iRF <- function(x, y, xtest = NULL, ytest = NULL, 
                n_iter = 5, 
                ntree = 500, 
                n_core = 1, 
                mtry_select_prob = rep(1/ncol(x), ncol(x)), 
                keep_impvar_quantile = NULL, 
                interactions_return = NULL, 
                node_sample = list(subset=function(x) rep(TRUE, nrow(x)),
                                   wt=function(x) x$size_node), 
                cutoff_unimp_feature = 0, 
                class_id = 1, 
                rit_param = c(5, 100, 2), 
                varnames_grp = NULL, 
                n_bootstrap = 30, 
                verbose = TRUE, 
                ...){

n <- nrow(x)
p <- ncol(x)
keep_subset_var <- NULL
class.irf <- is.factor(y)

# Check whether iRF can be run for given inputs
if (!is.matrix(x) | ((!is.null(xtest)) & (!is.matrix(xtest))))
    stop('either x or xtest is not a matrix !')
if ((!is.numeric(x)) | ((!is.null(xtest) & (!is.numeric(xtest)))))
    stop('either x or xtest is not a numeric matrix!')
if ((ncol(x) < 2) & !is.null(interactions_return))
    stop('cannot find interaction - X has less than two columns!')
if (any(interactions_return > n_iter))
  stop('interactions to return greater than niter')
    
# initialize outputs
rf_list <- list()
if (!is.null(interactions_return)){
    interact_list <- list()
    Stability_Score <- list()
}

# set number of trees to grow in each core
a <- floor(ntree/n_core) 
b <- ntree %% n_core
ntree_id <- c(rep(a+1, b), rep(a, n_core-b))


for (iter in 1:n_iter){
## 1: run random forest
   print(paste('iteration = ', iter))

   rf_list[[iter]] <- foreach(nt=ntree_id, .combine=combine, 
                              .multicombine=TRUE, .packages='iRF') %dopar% {
     randomForest(x, y, xtest, ytest, 
                  ntree=nt, 
                  mtry_select_prob=mtry_select_prob, 
                  keep_subset_var=keep_subset_var, 
                  keep.forest=TRUE, 
                  ...)
   }

   ## 2: find interactions stable over bootstrap replicates
   if (iter %in% interactions_return){
      if (verbose){cat('finding interactions ... ')}

      Stability_Score[[iter]] <- list()      
      # 2.1: find interactions in  multiple bootstrap samples to assess stability
      interact_list[[iter]] <- mclapply(1:n_bootstrap, function(i_b) {
         if (verbose) cat(paste('b = ', i_b, ';  ', '\n', sep=''))

         if (class.irf) {
           n.class <- table(y)
           sample_id <- mapply(function(cc, nn) sample(which(y == cc), nn, replace=TRUE),
                               as.factor(names(n.class)), n.class)
           sample_id <- unlist(sample_id)
         } else {
           sample_id = sample(n, n, replace=TRUE)
         }

         #2.1.1: fit random forest
         rf_b <- randomForest(x[sample_id,], y[sample_id], xtest, ytest, 
                              ntree=ntree, 
                              mtry_select_prob=mtry_select_prob, 
                              keep_subset_var = keep_subset_var, 
                              keep.forest=TRUE, 
                              track.nodes=TRUE, 
                              ...)

         ints <- generalizedRIT(rf=rf_b, x=x[sample_id,], y=y[sample_id], 
                                subsetFun=node_sample$subset, 
                                wtFun=node_sample$wt,
                                class.irf=class.irf, 
                                class_id=class_id, 
                                varnames_grp=varnames_grp, 
                                rit_param=rit_param,
                                cutoff_unimp_feature=cutoff_unimp_feature) 

         
        return(ints)       

      }, mc.cores=n_core) # end for (i_b in ... )

   # 2.2: calculate stability scores of interactions
   if (!is.null(varnames_grp))
       vnames_new = unique(varnames_grp)
   else if (!is.null(colnames(x)))
       vnames_new = colnames(x)
   else
       vnames_new = seq(ncol(x))
   
   Stability_Score[[iter]] <- summarize_interact(interact_list[[iter]], varnames=vnames_new)      
   
   } # end if (find_interaction)
   
   ## 3: update mtry_select_prob 
   if (!class.irf) 
     mtry_select_prob <- rf$importance[,'IncNodePurity']
   else
     mtry_select_prob <- rf$importance[,'MeanDecreaseGini']

   
   if (!is.null(xtest) & class.irf){
       auroc <- auc(roc(rf$test$votes[,2], ytest))
       print(paste('AUROC: ', round(auroc, 2)))
   }

} # end for (iter in ... )


out <- list()
out$rf_list <- rf_list
if (!is.null(interactions_return)){
    out$interaction <- Stability_Score
}
    
return(out)
}


generalizedRIT <- function(rf, x, y,
                           subsetFun, 
                           wtFun, 
                           class.irf, 
                           class_id, 
                           varnames_grp,
                           cutoff_unimp_feature, 
                           rit_param) {
  
  rforest <- readForest(rf, X=x, 
                        return_node_feature=TRUE, 
                        return_node_data=FALSE, 
                        subsetFun=subsetFun, 
                        wtFun=wtFun)

  # Select class specific leaf nodes if classification
  select_leaf_id <- rep(TRUE, nrow(rforest$tree_info))
  if (class.irf & all(y %in% c(0, 1)))
    select_leaf_id <- rforest$tree_info$prediction == as.numeric(class_id) + 1

  # Apply random intersection trees
  rforest <- subsetReadForest(rforest, select_leaf_id)
  nf <- rforest$node_feature

  if (sum(select_leaf_id) < 2){
    return(character(0))
  } else {

    # group features if specified
    if (!is.null(varnames_grp)) nf <- groupFeature(nf, grp = varnames_grp)

    # drop feature, if cutoff_unimp_feature is specified
    if (cutoff_unimp_feature > 0){
      if (!class.irf)
        rfimp <- rf$importance[,'IncNodePurity']
      else
        rfimp <- rf$importance[,'MeanDecreaseGini']

      drop_id <- which(rfimp < quantile(rfimp, prob=cutoff_unimp_feature))
      nf[,drop_id] <- FALSE
    }

    interactions <- RIT(nf, depth=rit_param[1], n_trees=rit_param[2], branch=rit_param[3])
    interactions <- interactions$Interaction
    interactions <- gsub(' ', '_', interactions)
    return(interactions)
  }
}

subsetReadForest <- function(rforest, subset_idcs) {
  # subset nodes from readforest output 
  if (!is.null(rforest$node_feature)) 
    rforest$node_feature <- rforest$node_feature[subset_idcs,]
  if (!is.null(rforest$node_data)) 
    rforest$node_data <- t(rforest$node_data[,subset_idcs])
  if(!is.null(rforest$tree_info))
  rforest$tree_info <- rforest$tree_info[subset_idcs,]
  return(rforest)
}
