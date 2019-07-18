context('test interactPredict')


test_that('sampleTree works', {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  info <- read.forest$tree.info

  leaf.id <- sampleTree(42, info$tree, info$size.node)
  expect_equal(length(leaf.id), 1)
})

test_that('interactPredict works', {
  set.seed(0)
  x <- iris[, -5]
  int <- 'Petal.Length+_Petal.Width+'
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  info <- read.forest$tree.info

  ip <- interactPredict(x, int, read.forest)
  expect_equal(length(ip), nrow(x))
  expect_equal(sum(ip), 85.920177383592)
})

