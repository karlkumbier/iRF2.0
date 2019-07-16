context('test interactImportance')

test_that("nodePrecision works", {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  count <- read.forest$tree.info$size.node
  np <- nodePrecision(read.forest, iris[, 5], count)
  expect_equal(length(np),
               nrow(read.forest$tree.info))
})

test_that("intImportance works", {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  nf.list <- by(read.forest$node.feature@i,
                rep(1:ncol(read.forest$node.feature),
                    times=diff(read.forest$node.feature@p)), list)

  int.eval <- c(5L, 8L)

  count <- read.forest$tree.info$size.node
  idcnt <- count >= 1
  read.forest <- subsetReadForest(read.forest, idcnt)
  count <- count[idcnt]

  precision <- nodePrecision(read.forest, iris[, 5], count)

  ximp <- intImportance(int.eval, nf.list, precision, count)

  expect_true('data.table' %in% class(ximp))
  expect_equal(names(ximp),
               c("prev1", "prev0", "prec"))
})

test_that("subsetTest works", {
  rand.forest <- randomForest(Species ~ ., iris)
  read.forest <- readForest(rand.forest, x=iris[, -5])
  nf.list <- by(read.forest$node.feature@i,
                rep(1:ncol(read.forest$node.feature),
                    times=diff(read.forest$node.feature@p)), list)
  int.eval <- c(5L, 8L)
  ints.sub <- list(5L, 7L, 8L,
                   c(5L, 7L), c(5L, 8L), c(7L, 8L),
                   c(5L, 7L, 8L))
  count <- read.forest$tree.info$size.node
  idcnt <- count >= 1
  read.forest <- subsetReadForest(read.forest, idcnt)
  count <- count[idcnt]
  precision <- nodePrecision(read.forest, iris[, 5], count)
  ximp <- lapply(ints.sub, intImportance, nf=nf.list, weight=count,
                 precision=precision)
  ximp <- rbindlist(ximp)

  imp.test <- subsetTest(int.eval, ints.sub, ximp)
  expect_true('data.table' %in% class(imp.test))
  expect_equal(names(imp.test),
               c("prev.test", "prec.test"))
})

