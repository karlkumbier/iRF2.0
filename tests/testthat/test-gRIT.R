context('test gRIT')

library(randomForest)
library(ranger)


test_that("gRIT works for randomForest", {
  rand.forest <- randomForest(Species ~ ., iris)
  fit <- gRIT(iris[, -5], iris[, 5], rand.forest)
  expect_true('data.table' %in% class(fit))
})

test_that("gRIT works for ranger", {
  x <- iris[, -5]
  y <- iris[, 5]
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  rand.forest <- ranger(data=cbind(x, y),
               dependent.variable.name='y',
               classification=class.irf)
  fit <- gRIT(x, y, rand.forest)
  expect_true('data.table' %in% class(fit))
})

test_that("runRIT works for randomForest", {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  rit.param <- list(depth=5, ntree=500, nchild=2,
                    class.id=1, min.nd=1, class.cut=NULL)
  weights <- rep(1, length(read.forest$tree.info$size.node))

  interactions <- runRIT(read.forest, weights, rit.param, 1)
  expect_equal(mode(interactions), 'list')
  expect_gte(length(interactions), 3)
})

test_that("runRIT works for ranger", {
  x <- iris[, -5]
  y <- iris[, 5]
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  rand.forest <- ranger(data=cbind(x, y),
               dependent.variable.name='y',
               classification=class.irf)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  rit.param <- list(depth=5, ntree=500, nchild=2,
                    class.id=1, min.nd=1, class.cut=NULL)
  weights <- rep(1, length(read.forest$tree.info$size.node))

  interactions <- runRIT(read.forest, weights, rit.param, 1)
  expect_equal(mode(interactions), 'list')
  expect_gte(length(interactions), 3)
})

