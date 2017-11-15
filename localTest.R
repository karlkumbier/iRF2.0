library(iRF)
set.seed(47)

# Generate 50-dim normals with two active components
# split between well separated populations
x1.i <- matrix(rnorm(500 * 48) ^ 2, nrow=500)
x1.a <- matrix(rnorm(500 * 2, 5) ^ 2, nrow=500)
y1 <- as.numeric(x1.a[,1] > 20 & x1.a[,2] > 20)

x2.i <- matrix(rnorm(500 * 48) ^ 2, nrow=500)
x2.a <- matrix(rnorm(500 * 2, 5) ^ 2, nrow=500)
y2 <- as.numeric(x2.a[,1] > 20 & x2.a[,2] > 20)

x <- rbind(cbind(x1.a, x1.i), cbind(x2.i, x2.a))
y <- as.factor(c(y1, y2))

f <- iRF(x=x, y=y, n.iter=3, n.bootstrap=5, 
  n.core=4, interactions.return=c(1, 3), local=TRUE,
  get.prevalence=TRUE)

# Get interaction prevalence among decision paths for each 
# observation 
rf <- readForest(f$rf.list[[3]], x=x, return.node.obs=TRUE)
int.name <- names(f$interaction[[3]])
ints <- lapply(strsplit(int.name, '_'), as.numeric)
prev.mat <- matrix(0, nrow=nrow(x), ncol=length(ints))
colnames(prev.mat) <- int.name
for (i in 1:nrow(x)) {
  id <- rf$node.obs[,i]
  rfsub <- subsetReadForest(rf, id)
  for (j in 1:length(ints)) {
    jj <- ints[[j]]
    prev.mat[i,j] <- mean(apply(as.matrix(rfsub$node.feature[,jj]), MAR=1, all))
  }
}

# Calculate individual feature prevalences
selected <- unique(unlist(ints))
selected.prev <- apply(rf$node.feature[,selected], MAR=2, mean)

# Scale interaction prevalence by individual feature selsection 
# probability
scaling <- sapply(ints, function(z) prod(selected.prev[selected %in% z]))
prev.mat <- t(t(prev.mat) / scaling)
heatmap(prev.mat)

