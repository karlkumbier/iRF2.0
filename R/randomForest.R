"randomForest" <-
function(x, ...)
  UseMethod("randomForest")


parRF <- function(x, y, 
                  xtest=NULL, 
                  ytest=NULL, 
                  ntree=500,
                  mtry.select.prob=rep(1, ncol(x)), 
                  n.core=1, 
                  ...) {
  
  # Wrapper function to run randomForest in parallel using foreach and dorng

  # Set number of trees per RF for each core
  a <- floor(ntree / n.core)
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  suppressWarnings(
    rf <- foreach(i=1:length(ntree.id), .combine=combine, 
                  .multicombine=TRUE, .packages='iRF') %dorng% {
                    randomForest(x, y, xtest, ytest,
                                 ntree=ntree.id,
                                 mtry.select.prob=mtry.select.prob,
                                 keep.forest=TRUE,
                                 ...)                         
    }
  )
  return(rf)
}



