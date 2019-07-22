SUITE <- 'test-margin'


##############
#   Set up   #
##############

makeTarget('rand.forest.RF', {
  set.seed(0)
  rand.forest.RF <- randomForest::randomForest(Species ~ ., iris)
}, skip=TRUE)

makeTarget('read.forest.RF', {
  read.forest.RF <- readForest(rand.forest.RF, x=iris[, -5])
})

makeTarget('runRIT.RF', {
  set.seed(0)
  rit.param <- list(depth=5, ntree=500, nchild=2,
                    class.id=1, min.nd=1, class.cut=NULL)
  weights <- rep(1, length(read.forest.RF$tree.info$size.node))

  runRIT.RF <- runRIT(read.forest.RF, weights, rit.param, 1)
}, skip=TRUE)

makeTarget('gRIT.RF', {
  set.seed(0)
  gRIT.RF <- gRIT(iris[, -5], iris[, 5], rand.forest.RF)
}, skip=TRUE)

makeTarget('rand.forest.ranger', {
  set.seed(0)
  x <- iris[, -5]
  y <- iris[, 5]
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  rand.forest.ranger <- ranger::ranger(data=cbind(x, y),
                                       dependent.variable.name='y',
                                       classification=class.irf)
}, skip=TRUE)

makeTarget('read.forest.ranger', {
  read.forest.ranger <- readForest(rand.forest.ranger, x=iris[, -5])
}, skip=TRUE)

makeTarget('runRIT.ranger', {
  set.seed(0)
  rit.param <- list(depth=5, ntree=500, nchild=2,
                    class.id=1, min.nd=1, class.cut=NULL)
  weights <- rep(1, length(read.forest.ranger$tree.info$size.node))

  runRIT.ranger <- runRIT(read.forest.ranger, weights, rit.param, 1)
}, skip=TRUE)

makeTarget('gRIT.ranger', {
  set.seed(0)
  x <- iris[, -5]
  y <- iris[, 5]
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  gRIT.ranger <- gRIT(x, y, rand.forest.ranger)
}, skip=TRUE)


#############
#   Tests   #
#############

test_that("runRIT works for randomForest", {
  expect_equal(mode(runRIT.RF), 'list')
  expect_gte(length(runRIT.RF), 3)
})

test_that("runRIT works for ranger", {
  expect_equal(mode(runRIT.ranger), 'list')
  expect_gte(length(runRIT.ranger), 3)
})

test_that("gRIT works for randomForest", {
  expect_true('data.table' %in% class(gRIT.RF))
  expect_equal(names(gRIT.RF),
               c("prev1", "prev0", "prec",
                 "int", "recovered",
                 "prev.test", "prec.test"))
})

test_that("gRIT works for ranger", {
  expect_true('data.table' %in% class(gRIT.ranger))
  expect_equal(names(gRIT.ranger),
               c("prev1", "prev0", "prec",
                 "int", "recovered",
                 "prev.test", "prec.test"))
})

