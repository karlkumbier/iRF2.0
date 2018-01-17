## Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y,
                xtest=NULL, ytest=NULL,
                n.iter=5,
                ntree=500,
                n.core=1,
                n.cpupercore=1,
                mtry= if (!is.null(y) && !is.factor(y))
                  max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
                mtry.select.prob=rep(1/ncol(x), ncol(x)),
                keep.impvar.quantile=NULL,
                interactions.return=NULL,
                wt.pred.accuracy=FALSE,
                cutoff.unimp.feature=0,
                rit.param=list(depth=4, ntree=50, nchild=2, class.id=1, class.cut=NULL),
                varnames.grp=NULL,
                n.bootstrap= if (n.core > 30) n.core else 30,
                bootstrap.forest=FALSE,
                escv.select=FALSE,
                verbose=TRUE,
                keep.subset.var=NULL,
                obs.weights=NULL,
                samplesize=nrow(x),
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
  if (n.cpupercore > 1) registerDoMC(n.cpupercore)
  n.totalcores=(n.core*n.cpupercore)
  rf.list <- list()
  if (!is.null(interactions.return) | escv.select) {
    interact.list <- list()
    stability.score <- list()
    prevalence <- list()
  }

  weight.mat <- matrix(0, nrow=p, ncol=n.iter+1)
  weight.mat[,1] <- mtry.select.prob

  # Set number of trees to grow in each core
  a <- floor(ntree / n.totalcores)
  b <- ntree %% n.totalcores
  ntree.id <- a
  #ntree.id <- rep(a,n.core)
  for (iter in 1:n.iter){

    ## 1: Grow Random Forest on full data
    print(paste('iteration = ', iter))
    tree.idcs <- c(1, cumsum(ntree.id) + 1)
    tree.idcs <- lapply(2:length(tree.idcs), function(i)
      tree.idcs[i-1]:(tree.idcs[i] - 1))


    if (iter %in% interactions.return){
        trackforinteractions=TRUE
      }else{
        trackforinteractions=FALSE
      }

    #rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine,
    #                           .multicombine=TRUE, .packages='iRF') %dopar% {
    #                             randomForest(x, y,
    #                                          xtest, ytest,
    #                                          ntree=ntree.id[i],
    #                                          mtry.select.prob=mtry.select.prob,
    #                                          keep.forest=TRUE,
    #                                          track.nodes=trackforinteractions,
    #                                          keep.subset.var=keep.subset.var[tree.idcs[[i]]],
    #                                          ...)
    #                           }
    if (verbose) { print('begin pbdLapply and mclapply for randomForest')}
    if (n.totalcores >1 ){
    forestlist <- pbdLapply(1:n.core, function(i) {
                                  forests <- mclapply(1:n.cpupercore, function(j)
                                              {randomForest(x, y,
                                              xtest, ytest,
                                              ntree=ntree.id,
                                              mtry=mtry,
                                              mtry.select.prob=weight.mat[,iter],
                                              keep.forest=TRUE,
                                              track.nodes=trackforinteractions,
                                              sampsize=samplesize,
                                              ...)},
                                              mc.cores = n.cpupercore
                                            )
                                  forest <- do.call(combine,forests)
                                  return(forest)}
                            )

    if (verbose) { print('unlist and combine forestlist')}
    gatheredModels <- unlist(allgather(forestlist), recursive = FALSE)
    if (verbose) { print(paste('type of gatheredModels', typeof(gatheredModels)))}
    #if (verbose) { print(paste('type of gatheredModels[[1]]', typeof(gatheredModels[[1]])))}
    rf.list[[iter]] <- do.call(combine,gatheredModels)
    if (verbose) { print('finished unlist and combine forestlist')}
    #rf.list[[iter]] <- combine(unlist(forestlist))
    rm(forestlist)
    rm(gatheredModels)
    gc(verbose = FALSE)
    } else {
      rf.list[[iter]] <- randomForest(x, y,
                                                xtest, ytest,
                                                ntree=ntree.id[1],
                                                mtry=mtry,
                                                mtry.select.prob=weight.mat[,iter],
                                                keep.forest=TRUE,
                                                track.nodes=trackforinteractions,
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

  if (verbose) { print('iterations complete')}

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
      #interact.list[[iter]] <- lapply(1:n.bootstrap, function(i.b) {
      #for (i.b in 1:n.bootstrap) {
      gatheredRITs <- pbdLapply(1:n.bootstrap, function(i.b) {
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

            #rf.b <- foreach(i=1:length(ntree.id), .combine=combine,
            #      .multicombine=TRUE, .packages='iRF') %dopar% {
            #
            #        randomForest(x[sample.id,], y[sample.id],
            #                     xtest, ytest,
            #                     ntree=ntree.id[i],
            #                     mtry.select.prob=weight.mat[,iter],
            #                     keep.forest=TRUE,
            #                     track.nodes=TRUE,
            #                     ...)
            #    }


            if (n.totalcores >1 ){
            forestlist <- mclapply(1:length(ntree.id), function(i)
                                                      randomForest(x[sample.id], y[sample.id],
                                                      xtest, ytest,
                                                      ntree=ntree.id[i],
                                                      mtry=mtry,
                                                      mtry.select.prob=weight.mat[,iter],
                                                      keep.forest=TRUE,
                                                      track.nodes=TRUE,
                                                      ...),
                                                      mc.cores = n.cpupercore
                                  )


            rf.b <- do.call(combine,forestlist)

            #rf.list[[iter]] <- combine(unlist(forestlist))
            rm(forestlist)
            gc(verbose = FALSE)
            } else {
              rf.b <- randomForest(x[sample.id], y[sample.id],
                                                        xtest, ytest,
                                                        ntree=ntree.id[1],
                                                        mtry=mtry,
                                                        mtry.select.prob=weight.mat[,iter],
                                                        keep.forest=TRUE,
                                                        track.nodes=TRUE,
                                                        ...)
            }

        } else {
          rf.b <- rf.list[[iter]]
        }
        #2.1.2: run generalized RIT on rf.b to learn interactions
        if (verbose){cat('before generalizedRIT ... ')}
        if (is.null(obs.weights)) {
          ints <- generalizedRIT(rf=rf.b,
                                 x=x[sample.id,], y=y[sample.id],
                                 wt.pred.accuracy=wt.pred.accuracy,
                                 class.irf=class.irf,
                                 varnames.grp=varnames.grp,
                                 cutoff.unimp.feature=cutoff.unimp.feature,
                                 rit.param=rit.param,
                                 obs.weight=obs.weights,
                                 n.core=n.core)
        } else {
          ints <- lapply(obs.weights, function(w) {
                           generalizedRIT(rf=rf.b,
                                          x=x[sample.id,], y=y[sample.id],
                                          wt.pred.accuracy=wt.pred.accuracy,
                                          class.irf=class.irf,
                                          varnames.grp=varnames.grp,
                                          cutoff.unimp.feature=cutoff.unimp.feature,
                                          rit.param=rit.param,
                                          obs.weight=w,
                                          n.core=n.core)
                                 })
        }

        interact.list.b[[i.b]] <- ints
        rm(rf.b)
        if (verbose){cat('after generalizedRIT ... ')}
        return(ints)
      })
      if (verbose){cat('gathering interactions ... ')}
      interact.list[[iter]] <- unlist(allgather(gatheredRITs),recursive = FALSE)
      interact.list.b <- interact.list[[iter]]
      # 2.2: calculate stability scores of interactions
      if (!is.null(varnames.grp))
        varnames.new <- unique(varnames.grp)
      else if (!is.null(colnames(x)))
        varnames.new <- colnames(x)
      else
        varnames.new <- 1:ncol(x)

      if (verbose){cat('summarizing interactions ... ')}
      if (is.null(obs.weights)) {
        interact.list[[iter]] <- interact.list.b
        summary.interact <- summarizeInteract(interact.list[[iter]],
                                              varnames=varnames.new)
        stability.score[[iter]] <- summary.interact$interaction
      } else {
        interact.list[[iter]] <- interact.list.b
        interact.list[[iter]] <- lapply(1:length(interact.list[[iter]][[1]]), function(i) {
                                           lapply(interact.list[[iter]], function(ll) {
                                                    return(ll[[i]])})
                                          })
        summary.interact <- lapply(interact.list[[iter]], summarizeInteract,
                                   varnames=varnames.new)
        stability.score[[iter]] <- lapply(summary.interact, function(i) i$interaction)
        #prevalence[[iter]] <- lapply(summary.interact, function(i) i$prevalence)

      }

    } # end if (find_interaction)


  } # end for (iter in ... )


  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)){
    out$interaction <- stability.score
    #out$prevalence <- prevalence
  }

  if (escv.select) {
    out$rf.list <- out$rf.list[[opt.k]]
    out$interaction <- out$interaction[[opt.k]]
    out$opt.k <- opt.k
    #out$prevalence <- out$prevalence[[opt.k]]
    out$weights <- weight.mat[,opt.k]
  }
  return(out)
}

generalizedRIT <- function(rf, x, y, wt.pred.accuracy, class.irf, varnames.grp,
                           cutoff.unimp.feature, rit.param, obs.weight, n.core) {

  # Extract decision paths from rf as sparse binary matrix to be passed to RIT
  #print('before readForest')
  rforest <- readForest(rf, x=x, y=y,
                        return.node.feature=TRUE,
                        wt.pred.accuracy=wt.pred.accuracy,
                        n.core=n.core)
  class.id <- rit.param$class.id
  #print('after readForest')
  # Select class specific leaf nodes
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  if (class.irf) {
    select.leaf.id <- rforest$tree.info$prediction == as.numeric(class.id) + 1
  } else if (is.null(rit.param$class.cut)) {
    select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  } else {
    select.leaf.id <- rforest$tree.info$prediction > rit.param$class.cut
    if (sum(select.leaf.id) < 2)
      select.leaf.id <- rforest$tree.info$prediction > rit.param$class.cut / 2
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
    #print('before groupFeature')
    if (!is.null(varnames.grp)) nf <- groupFeature(nf, grp=varnames.grp)

    # drop feature if cutoff.unimp.feature is specified
    if (cutoff.unimp.feature > 0){
      if (!class.irf)
        rfimp <- rf$importance[,'IncNodePurity']
      else
        rfimp <- rf$importance[,'MeanDecreaseGini']
      drop.id <- which(rfimp < quantile(rfimp, prob=cutoff.unimp.feature))
      nf[,drop.id] <- FALSE
    }
    #print(paste('before registerDoMC with cores: ',n.core))
    if (n.core > 1) registerDoMC(n.core)
    #print('before RIT')
    interactions <- RIT(nf, weights=wt, depth=rit.param$depth,
                        n_trees=rit.param$ntree, branch=rit.param$nchild,
                        n_cores=n.core)
    #print('after RIT')
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
      return(paste(varnames[nn], collapse='_'))
    })
    names(int.tbl) <- names.int
    names(prev.tbl) <- names.int

  }
  out <- list(interaction=int.tbl, prevalence=prev.tbl)
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
