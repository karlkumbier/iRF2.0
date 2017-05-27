library(dplyr)

# subset node feature so that it contains 1 row per leaf node
rd.f <- readForest(rf, x=x)
tree.info <- rd.f$tree.info
nf <- rd.f$node.feature
nf.sub <- cumsum(tree.info$size_node) + 1
nf.sub <- c(1, nf.sub[-length(nf.sub)])
nf <- nf[nf.sub,]

# get leaf node indices for oob observations in each tree
oob <- apply(rf$inbag, MAR=2, function(z) which(z == 0))
obs.nodes.split <- split(rf$obs.nodes, rep(1:rf$ntree, each=n))
if (is.factor(y)) y <- as.numeric(y) - 1
ob.nodes <- mapply(function(ln, ob, tr) {
  data.frame(node=ln[ob], y=y[ob], tree=tr, tree.purity=sd(y[ob]),
             lf.node=match(ln[ob], sort(unique(ln))))
  }, obs.nodes.split, oob, 1:length(oob), SIMPLIFY=FALSE)
ob.nodes <- do.call(rbind, ob.nodes)

# Determine best nodes based on size and purity
n <- length(y)
purity <- function(y) ifelse(length(y) == 1, 0, sd(y))
decreasePurity <- function(tp, y) max((tp - purity(y))/tp, 0)
nd.purity <- group_by(ob.nodes, tree, lf.node) %>%
  summarise(dec.purity=decreasePurity(tree.purity, y), size.node=length(y))

# Weight nodes by decrease in variability and size
wt <- nd.purity$size.node * nd.purity$dec.purity ^ 2
nd.ord <- order(wt, decreasing = T)
best.nodes <- nd.purity[nd.ord,]

getNodeFeat <- function(nf, node, tr) {
  offset <- ifelse(tr==1, 0, sum(n.node.tr[1:(tr-1)]))
  idx <- offset + node
  feat <- nf[idx,]
}

ff <- mapply(function(nd, tt) getNodeFeat(nf, nd, tt), 
             best.nodes$lf.node, best.nodes$tree)

# Sample subset of nf nodes
n.node.tr <- table(tree.info$tree)
tree.offset <- cumsum(n.node.tr)
tree.offset <- c(0, tree.offset[-rf$ntree])
nf.keep <- sapply(1:rf$ntree, function(i) {
  return(nd.purity$lf.node[nd.purity$tree == i] + tree.offset[i])
})
nf.keep <- unlist(nf.keep)
nf.samp <- sample(1:length(wt), 50000, prob=wt, replace=TRUE)
nf.ss <- nf.keep[nf.samp]
nf.subset <- nf[nf.ss,]
rit <- RIT(nf.subset)
