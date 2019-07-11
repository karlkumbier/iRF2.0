context('test readForest')

library(randomForest)
library(ranger)


test_that('readForest works for randomForest', {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  test_readForest(rand.forest, read.forest)

  countLeaf <- function(k)
      sum(randomForest::getTree(rand.forest, k)[, 'status'] == -1)
  nleaves <- sum(sapply(1:rand.forest$ntree, countLeaf))

  expect_true('data.table' %in%
              class(read.forest$tree.info))
  expect_equal(nrow(read.forest$tree.info),
               nleaves)
  expect_equal(names(read.forest$tree.info),
               c('prediction', 'node.idx', 'parent', 'tree', 'size.node'))

  expect_true('dgCMatrix' %in%
              class(read.forest$node.feature))
  expect_equal(nrow(read.forest$node.feature),
               nleaves)
  expect_equal(ncol(read.forest$node.feature),
               2 * length(rand.forest$importance))

  expect_true('ngCMatrix' %in%
              class(read.forest$node.obs))
  expect_equal(nrow(read.forest$node.obs),
               length(rand.forest$predicted))
  expect_equal(ncol(read.forest$node.obs),
               nleaves)
  expect_equal(rowSums(read.forest$node.obs),
               rep(rand.forest$ntree, length(rand.forest$predicted)))
})

test_that('readForest works for ranger', {
  x <- iris[, -5]
  y <- iris[, 5]
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  rand.forest <- ranger(data=cbind(x, y),
               dependent.variable.name='y',
               classification=class.irf)
  read.forest <- readForest(rand.forest, x=iris[, -5])

  countLeaf <- function(k)
      sum(ranger::treeInfo(rand.forest, k)$terminal)
  nleaves <- sum(sapply(1:rand.forest$num.trees, countLeaf))

  expect_true('data.table' %in%
              class(read.forest$tree.info))
  expect_equal(nrow(read.forest$tree.info),
               nleaves)
  expect_equal(names(read.forest$tree.info),
               c('prediction', 'node.idx', 'parent', 'tree', 'size.node'))

  expect_true('dgCMatrix' %in%
              class(read.forest$node.feature))
  expect_equal(nrow(read.forest$node.feature),
               nleaves)
  expect_equal(ncol(read.forest$node.feature),
               2 * rand.forest$num.independent.variables)

  expect_true('ngCMatrix' %in%
              class(read.forest$node.obs))
  expect_equal(nrow(read.forest$node.obs),
               rand.forest$num.samples)
  expect_equal(ncol(read.forest$node.obs),
               nleaves)
  expect_equal(rowSums(read.forest$node.obs),
               rep(rand.forest$num.trees, rand.forest$num.samples))
})

