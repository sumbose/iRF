context("test iRF")

n <- 100
p <- 10
x <- matrix(rnorm(n * p), nrow=n)
y <- as.numeric(rowMeans(x[,1:2] > 0) == 1)

testRF <- function(fit) {
  # Test iRF output for correct returns
  expect_equal(length(fit), 4)
  expect_true(class(fit$rf.list) %in% c('randomForest', 'ranger'))
  expect_equal(class(fit$interaction)[1], 'data.table')
}

testgRIT <- function(grit) {
  expect_equal(class(grit), 'data.frame')
  cols <- c("prev1", "prev0", "prec", "int", "recovered", 
            "prev.test", "prec.test")
  expect_equal(names(grit), cols)
}

test_that("classification iRF works", {
  skip('Not today')
  fit <- iRF(x=x, y=as.factor(y), n.bootstrap=3)
  testRF(fit)
  expect_error(iRF(x=x[1:10,], y=as.factor(y)),
      regexp='x and y must contain the same number of observations')
  expect_error(iRF(x=as.matrix(x[,1]), y=as.factor(y)),
      regexp='cannot find interaction - x has less than two columns!')
  expect_error(iRF(x=x, y=as.factor(y), interactions.return=10),
      regexp='selected iteration greater than n.iter')
  expect_error(iRF(x=x, y=as.factor(y), mtry.select.prob=1:2),
      regexp='length mtry.select.prob must equal number of features')
})

test_that("ranger classification iRF works", {
  skip('Not today')
  fit <- iRF(x=x, y=as.factor(y), n.bootstrap=3, type='ranger')
  testRF(fit)
  expect_error(iRF(x=x[1:10,], y=as.factor(y)),
      regexp='x and y must contain the same number of observations')
  expect_error(iRF(x=as.matrix(x[,1]), y=as.factor(y)),
      regexp='cannot find interaction - x has less than two columns!')
  expect_error(iRF(x=x, y=as.factor(y), interactions.return=10),
      regexp='selected iteration greater than n.iter')
  expect_error(iRF(x=x, y=as.factor(y), mtry.select.prob=1:2),
      regexp='length mtry.select.prop must equal number of features')
})

test_that("regression iRF works", {
  skip('Not today')
  yreg <- y + rnorm(n, sd=0.5)
  fit1 <- iRF(x=x, y=yreg, n.bootstrap=3)

  rit.param <- list(depth=5, ntree=500, nchild=2, class.id=1,
                    min.nd=1, class.cut=quantile(y, 0.1))
  fit2 <- iRF(x=x, y=yreg, n.bootstrap=3, rit.param=rit.param)

  testRF(fit1)
  testRF(fit2)
})

test_that("readForest works", {
  skip('Not today')
  fit <- iRF(x=x, y=as.factor(y), n.bootstrap=3)$rf.list
  nleaf <- sum(fit$forest$nodestatus == -1)
  read.forest <- readForest(fit, x=x)
  expect_equal(dim(read.forest$node.feature), c(nleaf, 2 * p))
  expect_equal(dim(read.forest$node.obs), c(n, nleaf))
  expect_equal(dim(read.forest$tree.info), c(nleaf, 5))
})

test_that("gRIT works", {
  skip('Not today')
  fit <- iRF(x=x, y=as.factor(y), n.bootstrap=3)$rf.list
  read.forest <- readForest(fit, x=x, return.node.obs=TRUE)
  grit1 <- gRIT(x=x, y=as.factor(y), rand.forest=fit)
  grit2 <- gRIT(x=x, y=as.factor(y), read.forest=read.forest)
  testgRIT(grit1)
  testgRIT(grit2)
})

