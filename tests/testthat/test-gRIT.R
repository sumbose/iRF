suite <- 'test-gRIT'


x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)
rit.param <- list(depth=5, ntree=500, nchild=2,
                  class.id=1, min.nd=1, class.cut=NULL)

for (RF in names(RF.collection)) {
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)

  rand.forest <- RF.collection[[RF]]
  read.forest %<-% readForest(rand.forest, x=x, oob.importance=FALSE)
  weights <- rep(1, length(read.forest$tree.info$size.node))

  test_that(paste('runRIT works for', RF), {
    set.seed(42L)
    runRIT.RF %<-% runRIT(read.forest, weights, rit.param, 1)
    expect_equal(mode(runRIT.RF), 'list')
    expect_gte(length(runRIT.RF), 3)
  })

  test_that(paste('signed gRIT works for', RF), {
    set.seed(42L)
    gRIT.signed %<-% gRIT(x, y, rand.forest,
                          oob.importance=FALSE)
    expect_true('data.table' %in% class(gRIT.signed))
    expect_equal(names(gRIT.signed),
                 c('prev1', 'prev0', 'prec',
                   'int.idx', 'int', 'recovered',
                   'prev.test', 'prec.test'))
    expect_gt(nrow(gRIT.signed), 3)
  })

  test_that(paste('unsigned gRIT works for', RF), {
    set.seed(42L)
    gRIT.unsigned %<-% gRIT(x, y, rand.forest, signed=FALSE,
                            oob.importance=FALSE)
    expect_true('data.table' %in% class(gRIT.unsigned))
    expect_equal(names(gRIT.unsigned),
                 c('prev1', 'prev0', 'prec',
                   'int.idx', 'int', 'recovered',
                   'prev.test', 'prec.test'))
    expect_gt(nrow(gRIT.unsigned), 3)
  })
}

