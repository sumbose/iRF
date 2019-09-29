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
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)

  test_that('signed iRF works in the first use case', {
    set.seed(42L)
    fit1.signed %<-% iRF(x=x, y=as.factor(y),
                         select.iter=TRUE, verbose=FALSE,
                         oob.importance=FALSE)
    expect_equal(names(fit1.signed),
                 c("rf.list", "selected.iter", "interaction", "weights"))
    expect_true(class(fit1.signed$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true(is.numeric(fit1.signed$selected.iter))
    expect_equal(length(fit1.signed$selected.iter), 1)
    expect_true('data.table' %in%
                class(fit1.signed$interaction))
    expect_equal(names(fit1.signed$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit1.signed$weights), 'numeric')
    expect_equal(dim(fit1.signed$weights), c(p, 1))
  })

  test_that('unsigned iRF works in the first use case', {
    set.seed(42L)
    fit1.unsigned %<-% iRF(x=x, y=as.factor(y), select.iter=TRUE,
                           verbose=FALSE, signed=FALSE,
                           oob.importance=FALSE)
    expect_equal(names(fit1.unsigned),
                 c("rf.list", "selected.iter", "interaction", "weights"))
    expect_true(class(fit1.unsigned$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true(is.numeric(fit1.unsigned$selected.iter))
    expect_equal(length(fit1.unsigned$selected.iter), 1)
    expect_true('data.table' %in%
                class(fit1.unsigned$interaction))
    expect_equal(names(fit1.unsigned$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit1.unsigned$weights), 'numeric')
    expect_equal(dim(fit1.unsigned$weights), c(p, 1))
  })
  
  test_that('signed iRF works in the second use case', {
    set.seed(42L)
    fit2.signed %<-% iRF(x=x, y=as.factor(y),
                         n.iter=5, int.return=5, verbose=FALSE,
                         oob.importance=FALSE)
    expect_equal(names(fit2.signed), c("rf.list", "interaction", "weights"))
    expect_true(class(fit2.signed$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true('data.table' %in%
                class(fit2.signed$interaction))
    expect_equal(names(fit2.signed$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit2.signed$weights), 'numeric')
    expect_equal(dim(fit2.signed$weights), c(p, 1))
  })
  
  test_that('unsigned iRF works in the second use case', {
    set.seed(42L)
    fit2.unsigned %<-% iRF(x=x, y=as.factor(y),
                          n.iter=5, int.return=5, verbose=FALSE, signed=FALSE,
                          oob.importance=FALSE)
    expect_equal(names(fit2.unsigned), c("rf.list", "interaction", "weights"))
    expect_true(class(fit2.unsigned$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_true('data.table' %in%
                class(fit2.unsigned$interaction))
    expect_equal(names(fit2.unsigned$interaction),
                 c("int", "prevalence", "precision", "cpe",
                   "sta.cpe", "fsd", "sta.fsd", "mip",
                   "sta.mip", "stability"))
    expect_equal(mode(fit2.unsigned$weights), 'numeric')
    expect_equal(dim(fit2.unsigned$weights), c(p, 1))
  })
  
  test_that('signed iRF works in the third use case', {
    set.seed(42L)
    fit3.signed %<-% iRF(x=x, y=as.factor(y), verbose=FALSE,
                         oob.importance=FALSE)
    expect_equal(names(fit3.signed), c("rf.list", "weights"))
    expect_true(class(fit3.signed$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_equal(mode(fit3.signed$weights), 'numeric')
    expect_equal(dim(fit3.signed$weights), c(p, 1))
  })
  
  test_that('unsigned iRF works in the third use case', {
    set.seed(42L)
    fit3.unsigned %<-% iRF(x=x, y=as.factor(y),
                           verbose=FALSE, signed=FALSE,
                           oob.importance=FALSE)
    expect_equal(names(fit3.unsigned), c("rf.list", "weights"))
    expect_true(class(fit3.unsigned$rf.list) %in%
                c('randomForest', 'ranger'))
    expect_equal(mode(fit3.unsigned$weights), 'numeric')
    expect_equal(dim(fit3.unsigned$weights), c(p, 1))
  })

  test_that('stabilityScore works', {
    set.seed(42L)
    ss %<-% stabilityScore(x, y, oob.importance=FALSE)
    expect_equal(names(ss),
                 c("int", "prevalence", "precision",
                   "cpe", "sta.cpe", "fsd",
                   "sta.fsd", "mip", "sta.mip",
                   "stability"))
  })
  
  test_that('sampleClass works', {
    set.seed(42L)
    class.labels %<-% iris$Species
    sampled %<-% sampleClass(class.labels, cls, num)
    expect_equal(length(sampled), num)
    expect_equal(class.labels[sampled], rep(cls, num))
  }) 
  
  test_that('bsSample works for classification', {
    set.seed(42L)
    bsCls %<-% bsSample(iris$Species)
    sampledCls %<-% table(iris[bsCls, 5])
    expect_equal(as.vector(sampledCls), rep(50L, 3))
  })
  
  test_that('bsSample works for regression', {
    set.seed(42L)
    bsReg %<-% bsSample(mtcars$mpg)
    expect_equal(length(bsReg), nrow(mtcars))
  })
}

