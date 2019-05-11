#' random forest
#'
#' @export
#'
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
#' @importFrom ranger ranger
"randomForest" <- function(x, ...) UseMethod("randomForest")

parRF <- function(x, y, xtest=NULL, ytest=NULL, ntree=500,
                  n.core=1, mtry.select.prob=rep(1, ncol(x)),
                  type='randomForest', ...) {
  
  # Wrapper function to run RF in parallel using randomForest or ranger
  if (type == 'randomForest') {
    rf <- randomForestPar(x, y, xtest, ytest, ntree, n.core, 
                          mtry.select.prob, ...)
  } else if (type == 'ranger') {
    rf <- rangerPar(x, y, xtest, ytest, ntree, n.core, 
                    mtry.select.prob, ...)
  } else {
    stop('type must be one of "randomForest" or "ranger"')
  }

  return(rf)
}

rangerPar <- function(x, y, xtest=NULL, ytest=NULL, ntree=500,
                      n.core=1, mtry.select.prob=rep(1, ncol(x)),
                      ...) {
  
  # Run feature weighted ranger in parallel
  rf <- ranger(y ~ ., data=data.frame(x, y), num.trees=ntree, 
               num.threads=n.core, importance='impurity', 
               split.select.weights=mtry.select.prob, ...)

  return(rf)
}

randomForestPar <- function(x, y, xtest=NULL, ytest=NULL, ntree=500, 
                            n.core=1, mtry.select.prob=rep(1, ncol(x)), 
                            ...) {  
  
  # Run randomForest in parallel using foreach and dorng
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)

  # Set number of trees per RF for each core
  a <- floor(ntree / n.core)
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  suppressWarnings(
    rf <- foreach(i=1:length(ntree.id), .combine=combine, 
                  .multicombine=TRUE, .packages='iRF') %dorng% {
                    randomForest(x, y, xtest, ytest,
                                 ntree=ntree.id[i],
                                 mtry.select.prob=mtry.select.prob,
                                 keep.forest=TRUE,
                                 ...)                         
    }
  )

  stopImplicitCluster()
  return(rf)
}



