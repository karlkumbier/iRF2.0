gRIT <- function(x, y, 
                 rand.forest=NULL, 
                 read.forest=NULL,
                 weights=rep(1, nrow(x)),
                 varnames.grp=colnames(x),
                 rit.param=list(depth=5,
                         ntree=500,
                         nchild=2,
                         class.id=1,
                         min.nd=1,
                         class.cut=NULL),
                 signed=TRUE,
                 ints.full=NULL,
                 n.core=1) {

  class.irf <- is.factor(y)
  if (n.core == -1) n.core <- detectCores()  
  if (n.core > 1) registerDoParallel(n.core)

  # Check rit parameters and set default values if not specified
  if (is.null(rand.forest) & is.null(read.forest))
    stop('Supply random forest or read forest output')
  if (is.null(rit.param$depth)) 
    rit.param$depth <- 5
  if (is.null(rit.param$ntree)) 
    rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) 
    rit.param$nchild <- 2
  if (is.null(rit.param$class.id) & class.irf) 
    rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) 
    rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & !class.irf) 
    rit.param$class.cut <- median(y)
  if (!class.irf) 
    y <- as.numeric(y >= rit.param$class.cut)


  # Set feature names for grouping interactions
  if (is.null(varnames.grp) & !is.null(colnames(x)))
    varnames.grp <- colnames(x)
  else if (is.null(varnames.grp))
    varnames.grp <- as.character(1:ncol(x))
  
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)

  print('read forest')
  print(detectCores())

  # Read RF object to extract decision path metadata
  if (is.null(read.forest)) {
    read.forest <- readForest(rand.forest, x=x,
                              return.node.feature=TRUE,
                              return.node.obs=TRUE,
                              varnames.grp=varnames.grp,
                              get.split=TRUE,
                              n.core=n.core)
  }

  print('summarize nodes')
  print(detectCores())
  # Collapse node feature matrix for unsigned iRF
  if (!signed) read.forest$node.feature <- collapseNF(read.forest$node.feature)

  # Evaluate leaf node attributes
  nd.attr <- nodeAttr(read.forest, y, weights)
  prec.nd <- nd.attr$precision
  ndcnt <- nd.attr$ndcnt
  
  # Subset leaf nodes based on mimimum size
  idcnt <- ndcnt >= rit.param$min.nd
  read.forest <- subsetReadForest(read.forest, idcnt)
  ndcnt <- ndcnt[idcnt]
  prec.nd <- prec.nd[idcnt]

  # Select class specific leaf nodes
  if (class.irf)
    idcl <- read.forest$tree.info$prediction == rit.param$class.id + 1
  else
    idcl <- read.forest$tree.info$prediction > rit.param$class.cut

  print('RIT')
  print(detectCores())
  if (sum(idcl) < 2) {
    return(nullReturn())
  } else {
  
    # Run RIT on leaf nodes of selected class  
    ints <- runRIT(subsetReadForest(read.forest, idcl), weights=ndcnt[idcl],
                  rit.param=rit.param, n.core=n.core)
    
    if (is.null(ints)) return(nullReturn())
    
    # Set recovered interactions or convert to indices if supplied
    if (is.null(ints.full)) {
      ints.full <- ints
    } else {
      ints.full <- lapply(ints.full, int2Id, signed=signed,
                          varnames.grp=varnames.grp)
    }
   
    # Evaluate importance metrics for supplied/recovered interactions and lower
    # order subsets.
    ints.sub <- lapply(ints.full, intSubsets)
    ints.sub <- unique(unlist(ints.sub, recursive=FALSE))
    
    print('Importance')
    print(detectCores())
    suppressWarnings(
    ximp <- foreach(int=ints.sub) %dorng% {
      intImportance(int, nf=read.forest$node.feature, weight=ndcnt,
                    prec.nd=prec.nd, select.id=idcl)
    })
    ximp <- rbindlist(ximp) 

    print('Test')
    imp.test <- lapply(ints.full, subsetTest, importance=ximp, ints=ints.sub)
    imp.test <- rbindlist(imp.test)
  }

  print('aggregate')
  # Aggregate evaluated interations for return
  ints.recovered <- nameInts(ints, varnames.unq, signed=signed)
  
  imp.test <- imp.test %>%
    mutate(int=nameInts(ints.full, varnames.unq, signed=signed))
  
  ximp <- ximp %>%
    mutate(int=nameInts(ints.sub, varnames.unq, signed=signed),
           recovered=int %in% ints.recovered) %>%
    right_join(imp.test, by='int')

  return(ximp)
}

runRIT <- function(read.forest, weights, rit.param, n.core=1) {
  # Run a weighted version of RIT across RF decision paths
  xrit <- cbind(read.forest$node.feature[weights > 0,])
  interactions <- RIT(xrit,
                      weights=weights[weights > 0],
                      depth=rit.param$depth,
                      n_trees=rit.param$ntree,
                      branch=rit.param$nchild,
                      output_list=TRUE,
                      n_cores=n.core)$Interaction
  
  return(interactions)
}

subsetReadForest <- function(read.forest, subset.idcs) {
  # Subset nodes from readforest output 
  if (!is.null(read.forest$node.feature))
    read.forest$node.feature <- read.forest$node.feature[subset.idcs,]

  if (!is.null(read.forest$tree.info))
    read.forest$tree.info <- read.forest$tree.info[subset.idcs,]

  if (!is.null(read.forest$node.obs))
    read.forest$node.obs <- read.forest$node.obs[subset.idcs,]

  return(read.forest)
}

collapseNF <- function(x) {
  # Convert from signed to unsigned node feature matrix
  p <- ncol(x) / 2
  x <- x[,1:p] + x[,(p + 1):(2 * p)]
  return(x)
}

nullReturn <- function() {
  # Return empty interaction and importance
  out <- list()
  out$int <- character(0)
  out$imp <- data.table(prev1=numeric(0), prev0=numeric(0),
                        prec=numeric(0), prev.test=numeric(0),
                        prec.test=numeric(0), int=character(0))
  return(out)
}
