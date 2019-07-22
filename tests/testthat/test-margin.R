SUITE <- 'test-margin'


##############
#   Set up   #
##############

makeTarget('randomForest_v', {
  set.seed(0)
  rand.forest <- randomForest::randomForest(iris[, -5], iris[, 5])
}, skip=TRUE)

makeTarget('margin_v', {
  margin.obs <- margin(randomForest_v)
})

makeTarget('plot_margin_v', {
  postscript(file='plot_margin.ps')
  plot(margin_v)
  dev.off()

  hash.value <- tools::md5sum('plot_margin.ps')
})


#############
#   Tests   #
#############

test_that("margin works for randomForest", {
  expect_equal(length(margin_v), 150)
})

