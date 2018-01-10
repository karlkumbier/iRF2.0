unit_test_individual_paths_frequent_feature_interactions <- function(p = 50, n = 500){
  library(iRF)
  set.seed(47)
  thres = 20
  # Generate 50-dim normals with two active components
  # split between well separated populations
  x1.i <- matrix(rnorm(n * (p - 2)) ^ 2, nrow=n)
  x1.a <- matrix(rnorm(n * 2, 5) ^ 2, nrow=n)
  y1 <- as.numeric(x1.a[,1] > thres & x1.a[,2] > thres)
  
  x2.i <- matrix(rnorm(n * (p - 2)) ^ 2, nrow=n)
  x2.a <- matrix(rnorm(n * 2, 5) ^ 2, nrow=n)
  y2 <- as.numeric(x2.a[,1] > thres & x2.a[,2] > thres)
  
  x <- rbind(cbind(x1.a, x1.i), cbind(x2.i, x2.a))
  y <- as.factor(c(y1, y2))
  
  f <- iRF(x=x, y=y, n.iter=3, n.bootstrap=5, 
           n.core=4, interactions.return=c(1, 3), local=TRUE,
           get.prevalence=TRUE)
  out = list()
  for (i in 1:10){
    out[[i]] = individual_paths_frequent_feature_interactions(f, x, i)
  }
  return(out)
}
unit_test_individual_paths_frequent_feature_interactions()
