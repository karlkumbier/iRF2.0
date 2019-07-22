SUITE <- 'test-margin'


##############
#   Set up   #
##############


makeTarget('rand.forest', {
  set.seed(0)
  rand.forest <- randomForest::randomForest(Species ~ ., iris)
}, skip=TRUE)

makeTarget('read.forest', {
  read.forest <- readForest(rand.forest, x=iris[, -5])
})

makeTarget('np', {
  count <- read.forest$tree.info$size.node
  np <- nodePrecision(read.forest, iris[, 5], count)
})

makeTarget('ximp', {
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
})

makeTarget('imp.test', {
  ints.sub <- list(5L, 7L, 8L,
                   c(5L, 7L), c(5L, 8L), c(7L, 8L),
                   c(5L, 7L, 8L))
  imp.test <- subsetTest(int.eval, ints.sub, ximp)
})


#############
#   Tests   #
#############

test_that("nodePrecision works", {
  expect_equal(length(np),
               nrow(read.forest$tree.info))
})

test_that("intImportance works", {
  expect_true('data.table' %in% class(ximp))
  expect_equal(names(ximp),
               c("prev1", "prev0", "prec"))
})

test_that("subsetTest works", {
  expect_true('data.table' %in% class(imp.test))
  expect_equal(names(imp.test),
               c("prev.test", "prec.test"))
})

