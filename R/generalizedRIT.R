generalizedRIT <- function(x, rand.forest=NULL, read.forest=NULL,
                           weights=rep(1, nrow(x)),
                           varnames.grp=colnames(x),
                           rit.param=list(depth=5,
                                          ntree=500,
                                          nchild=2,
                                          class.id=1,
                                          min.nd=1,
                                          class.cut=NULL),
                           get.prevalence=FALSE,
                           int.sign=TRUE,
                           obs.rit=FALSE,
                           out.file=NULL,
                           n.core=1) {

  out <- list()
  if (is.null(rand.forest) & is.null(read.forest))
    stop('Supply random forest or read forest output')

  if (!is.null(rand.forest)) 
    class.irf <- rand.forest$type == 'classification'
  else
    class.irf <- all(read.forest$tree.info$prediction %in% 1:2)

  ntrain <- nrow(x)

  # Check all RIT and set to defaul if missing
  if (is.null(rit.param$depth)) rit.param$depth <- 5
  if (is.null(rit.param$ntree)) rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) rit.param$nchild <- 2
  if (is.null(rit.param$class.id) & class.irf) rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & !class.irf) 
    stop('Specifiy class.cut for regression')

  # Set feature names for grouping interactions
  if (is.null(varnames.grp)) varnames.grp <- as.character(1:ncol(x))
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)
  
  # Read RF object to extract decision path metadata
  if (is.null(read.forest)) {
    read.forest <- readForest(rand.forest, x=x,
                              return.node.feature=TRUE,
                              return.node.obs=TRUE,
                              varnames.grp=varnames.grp,
                              get.split=TRUE,
                              n.core=n.core)
  }

  # Collapse node feature matrix for unsigned iRF
  if (!int.sign) {
    read.forest$node.feature <- read.forest$node.feature[,1:p] +
      read.forest$node.feature[,(p + 1):(2 * p)]
  }

  if (!is.null(out.file)) save(file=out.file, read.forest, ntrain)

  # Select class specific leaf nodes
  pred <- read.forest$tree.info$prediction
  if (class.irf)
    select.id <- pred == rit.param$class.id + 1
  else
    select.id <- pred > rit.param$class.cut

  # Run RIT on leaf nodes of selected class to find
  # prevalent interactions
  if (sum(select.id) < 2) {
    return(character(0))
  } else {
    
    wt <- Matrix::colSums(t(read.forest$node.obs) * weights)
    out$int <- runRIT(subsetReadForest(read.forest, select.id),
                      weights=wt[select.id],
                      rit.param=rit.param,
                      obs.rit=obs.rit,
                      n.core=n.core)
    
    if (is.null(out$int)) return(NULL)

    # If running observation level RIT, split feature and observation
    # interacitons and group by feature interactions
    if (obs.rit) {
      pp <- ncol(read.forest$node.feature)
      out$obs <- lapply(out$int, function(z) pasteInt(z[z > pp] - pp))
      out$int <- sapply(out$int, function(z) pasteInt(z[z <= pp]))

      # Group observation interactions by feature interacitons
      dgroup <- data.table(int=out$int, obs=out$obs) %>%
        group_by(int) %>%
        summarize(obs=pasteInt(obs)) %>%
        filter(int != '')

      out$int <- str_split(dgroup$int, '_')
      out$int <- lapply(out$int, as.numeric)
      out$obs <- str_split(dgroup$obs, '_')
      out$obs <- lapply(out$obs, function(z) as.numeric(unique(z)))
    }


    if (get.prevalence) {
      out$prev <- mclapply(out$int, prevalence,
                           nf=read.forest$node.feature,
                           select.id=select.id, wt=wt,
                           mc.cores=n.core)
      out$prev <- rbindlist(out$prev) 
    }
     
    int.names <- nameInts(out$int, varnames.unq, signed=int.sign)
    out$int <- int.names
    if (get.prevalence) out$prev$int <- int.names
  }
  return(out)
}

runRIT <- function(read.forest, weights, rit.param, 
                   obs.rit=FALSE, n.core=1) {

  # Remove nodes below specified size threshold
  id.rm <- read.forest$tree.info$size.node < rit.param$min.nd
  if (all(id.rm)) {
    warning(paste('No nodes with greater than ', rit.param$min.nd,
                  'observations. Using all nodes.'))
    id.rm <- rep(FALSE, length(id.rm))
  }
  
  read.forest <- subsetReadForest(read.forest, !id.rm)
  wt <- weights[!id.rm]
  
  # Input for RIT: observations and features or only features 
  if (obs.rit) {
    xrit <- cbind(read.forest$node.feature[wt > 0,],
                  read.forest$node.obs[wt > 0,])
  } else {
    xrit <- cbind(read.forest$node.feature[wt > 0,])
  }


  interactions <- RIT(xrit,
                      weights=wt[wt > 0],
                      depth=rit.param$depth,
                      n_trees=rit.param$ntree,
                      branch=rit.param$nchild,
                      output_list=TRUE,
                      n_cores=n.core)$Interaction
  
  return(interactions)
}

nameInts <- function(ints, varnames, signed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by
  # variable names. Naming convention for an interaction is:
  #   <variable1(sign)>_ <variable2(sign)>_...
  
  varnames <- unique(varnames)
  p <- length(varnames)
  if (signed) 
    ints.signs <- lapply(ints, function(z) ifelse(z > p, '+', '-'))
  else 
    ints.signs <- ''
  
  # Adjust indexing to match varnames
  ints <- lapply(ints, function(z) z %% p + p * (z == p | z == 2 * p))
  ints.name <- mapply(function(i, s) paste0(varnames[i], s),
                      ints, ints.signs, SIMPLIFY=FALSE)
  ints.name <- sapply(ints.name, function(z) paste(sort(z), collapse='_'))
  return(ints.name)
}

prevalence <- function(int, nf, select.id, wt=rep(1, ncol(nf))) {
  # Calculate the decision path prevalence of an interaction
  
  id.rm <- wt == 0
  select.id <- select.id[!id.rm]
  wt <- wt[!id.rm]
  nf <- nf[!id.rm,]

  int <- as.numeric(int)
  intord <- length(int)
  if (intord == 1)
    int.id <- nf[, int] != 0
  else
    int.id <- Matrix::rowSums(nf[, int] != 0) == intord

  # Compute conditional probabilities of interaction given class and class given
  # interaction
  if (sum(int.id) == 0) 
    return(data.table(prev1=0, prev0=0, prop1=0))
  
  sint1 <- sum(wt[int.id & select.id])
  sint0 <- sum(wt[int.id & !select.id])
  s1 <- sum(wt[select.id])
  s0 <- sum(wt[!select.id])
  sint <- sum(wt[int.id])
  
  prev1 <- sint1 / s1
  prev0 <- sint0 / s0
  prop1 <- sint1 / sint
  return(data.table(prev1, prev0, prop1))
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

pasteInt <- function(x) paste(x, collapse='_')
