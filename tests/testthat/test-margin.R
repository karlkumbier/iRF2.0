SUITE <- 'test-margin'
REGENERATE <- FALSE


##############
#   Set up   #
##############

makeTarget('randomForest_v', {
  set.seed(0)
  randomForest(iris[, -5], iris[, 5])
})

makeTarget('margin_v', {
  margin(randomForest_v)
})

makeTarget('plot_margin_v', {
  postscript(file='plot_margin.ps')
  plot(margin_v)
  dev.off()

  tools::md5sum('plot_margin.ps')
})


#############
#   Tests   #
#############

test_that("margin works for randomForest", {
  expect_equal(length(margin_v), 150)
})

