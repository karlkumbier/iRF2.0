suite <- 'test-iRF'


n <- 100
p <- 10
x <- matrix(rnorm(n * p), nrow=n)
y <- as.numeric(rowMeans(x[,1:2] > 0) == 1)
RF.collection <- make.RF.collection(x, y)
cls <- structure(2L, .Label = c("setosa", "versicolor", "virginica"),
                 class = "factor")
num <- 50

for (RF in names(RF.collection)) {
  `%<-cache%` <- `%<-meta.cache%`(suite, RF, FALSE)
  `%<-verify%` <- `%<-meta.cache%`(suite, RF, TRUE)

  test_that('iRF works in the first use case', {
    fit1 %<-cache% iRF(x=x, y=as.factor(y), select.iter=TRUE, verbose=FALSE)
    expect_equal(names(fit1),
                 c("rf.list", "selected.iter", "interaction", "weights"))
    expect_true(class(fit1$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true(is.numeric(fit1$selected.iter))
    expect_equal(length(fit1$selected.iter), 1)
    expect_true('data.table' %in%
                class(fit1$interaction))
    expect_equal(names(fit1$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit1$weights), 'numeric')
    expect_equal(dim(fit1$weights), c(p, 1))
  })
  
  test_that('iRF works in the second use case', {
    fit2 %<-cache% iRF(x=x, y=as.factor(y), n.iter=5, int.return=5, verbose=FALSE)
    expect_equal(names(fit2), c("rf.list", "interaction", "weights"))
    expect_true(class(fit2$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true('data.table' %in%
                class(fit2$interaction))
    expect_equal(names(fit2$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit2$weights), 'numeric')
    expect_equal(dim(fit2$weights), c(p, 1))
  })
  
  test_that('iRF works in the third use case', {
    fit3 %<-cache% iRF(x=x, y=as.factor(y), verbose=FALSE)
    expect_equal(names(fit3), c("rf.list", "weights"))
    expect_true(class(fit3$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_equal(mode(fit3$weights), 'numeric')
    expect_equal(dim(fit3$weights), c(p, 1))
  })

  test_that('stabilityScore works', {
    ss %<-cache% stabilityScore(x, y)
    expect_equal(names(ss),
                 c("int", "prevalence", "precision",
                   "cpe", "sta.cpe", "fsd",
                   "sta.fsd", "mip", "sta.mip",
                   "stability"))
  })
  
  test_that('sampleClass works', {
    set.seed(42)
    class.labels %<-verify% iris$Species
    sampled %<-verify% sampleClass(class.labels, cls, num)
    expect_equal(length(sampled), num)
    expect_equal(class.labels[sampled], rep(cls, num))
  }) 
  
  test_that('bsSample works for classification', {
    set.seed(42)
    bsCls %<-verify% bsSample(iris$Species)
    sampledCls %<-verify% table(iris[bsCls, 5])
    expect_equal(as.vector(sampledCls), rep(50L, 3))
  })
  
  test_that('bsSample works for regression', {
    set.seed(42)
    bsReg %<-verify% bsSample(mtcars$mpg)
    expect_equal(length(bsReg), nrow(mtcars))
  })
}

