suite <- 'test-gRIT'


x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)
rit.param <- list(depth=5, ntree=500, nchild=2,
                  class.id=1, min.nd=1, class.cut=NULL)

for (RF in names(RF.collection)) {
  `%<-cache%` <- `%<-meta.cache%`(suite, RF, FALSE)
  `%<-verify%` <- `%<-meta.cache%`(suite, RF, TRUE)

  rand.forest <- RF.collection[[RF]]
  read.forest %<-verify% readForest(rand.forest, x=x)
  weights <- rep(1, length(read.forest$tree.info$size.node))

  test_that(paste('runRIT works for', RF), {
    runRIT.RF %<-cache% runRIT(read.forest, weights, rit.param, 1)
    expect_equal(mode(runRIT.RF), 'list')
    expect_gte(length(runRIT.RF), 3)
  })

  test_that(paste('signed gRIT works for', RF), {
    gRIT.signed %<-cache% gRIT(x, y, rand.forest)
    expect_true('data.table' %in% class(gRIT.signed))
    expect_equal(names(gRIT.signed),
                 c('prev1', 'prev0', 'prec',
                   'int.idx', 'int', 'recovered',
                   'prev.test', 'prec.test'))
  })

  test_that(paste('unsigned gRIT works for', RF), {
    gRIT.unsigned %<-cache% gRIT(x, y, rand.forest, signed=FALSE)
    expect_true('data.table' %in% class(gRIT.unsigned))
    expect_equal(names(gRIT.unsigned),
                 c('prev1', 'prev0', 'prec',
                   'int.idx', 'int', 'recovered',
                   'prev.test', 'prec.test'))
  })
}

