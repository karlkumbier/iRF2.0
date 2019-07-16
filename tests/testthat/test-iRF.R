context('iRF works')

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
