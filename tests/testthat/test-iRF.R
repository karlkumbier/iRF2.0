context('test iRF')

test_that('iRF works', {
  n <- 100
  p <- 10
  x <- matrix(rnorm(n * p), nrow=n)
  y <- as.numeric(rowMeans(x[,1:2] > 0) == 1)

  fit <- iRF(x=x, y=as.factor(y), select.iter=TRUE, verbose=FALSE)
  expect_equal(names(fit),
               c("rf.list", "selected.iter", "interaction", "weights"))
  expect_true(class(fit$rf.list) %in%
              c('randomForest', 'ranger'))
  expect_true(is.numeric(fit$selected.iter))
  expect_equal(length(fit$selected.iter), 1)
  expect_true('data.table' %in%
              class(fit$interaction))
  expect_equal(names(fit$interaction),
               c("int", "prevalence", "precision", "cpe",
                 "sta.cpe", "fsd", "sta.fsd", "mip",
                 "sta.mip", "stability"))
  expect_equal(mode(fit$weights), 'numeric')
  expect_equal(dim(fit$weights), c(p, 1))

  fit <- iRF(x=x, y=as.factor(y), n.iter=5, int.return=5, verbose=FALSE)
  expect_equal(names(fit), c("rf.list", "interaction", "weights"))
  expect_true(class(fit$rf.list) %in%
              c('randomForest', 'ranger'))
  expect_true('data.table' %in%
              class(fit$interaction))
  expect_equal(names(fit$interaction),
               c("int", "prevalence", "precision", "cpe",
                 "sta.cpe", "fsd", "sta.fsd", "mip",
                 "sta.mip", "stability"))
  expect_equal(mode(fit$weights), 'numeric')
  expect_equal(dim(fit$weights), c(p, 1))

  fit <- iRF(x=x, y=as.factor(y), verbose=FALSE)
  expect_equal(names(fit), c("rf.list", "weights"))
  expect_true(class(fit$rf.list) %in%
              c('randomForest', 'ranger'))
  expect_equal(mode(fit$weights), 'numeric')
  expect_equal(dim(fit$weights), c(p, 1))
})

test_that('stabilityScore works', {
  ss <- stabilityScore(iris[, -5], iris[, 5])
  expect_equal(names(ss),
               c("int", "prevalence", "precision",
                 "cpe", "sta.cpe", "fsd",
                 "sta.fsd", "mip", "sta.mip",
                 "stability"))
})

test_that('sampleClass works', {
  y <- iris$Species
  cl <- iris[51, 5]
  n <- 42
  sc <- sampleClass(y, cl, n)
  expect_equal(length(sc), n)
  expect_equal(y[sc], rep(cl, n))
}) 

test_that('bsSample works for classification', {
  bs <- bsSample(iris$Species)
  expect_equal(table(iris[bs, 5]),
               structure(c(setosa = 50L,
                           versicolor = 50L,
                           virginica = 50L),
                         .Dim = 3L,
                         .Dimnames = structure(list( c("setosa",
                                                       "versicolor",
                                                       "virginica")),
                                               .Names = ""),
                         class = "table"))
})

test_that('bsSample works for regression', {
  bs <- bsSample(mtcars$mpg)
  expect_equal(length(bs), nrow(mtcars))
})
