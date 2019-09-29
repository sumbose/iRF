suite <- 'test-interactionImportance'


x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)

int.eval <- c(5L, 8L)
ints.sub <- list(5L, 7L, 8L,
                 c(5L, 7L), c(5L, 8L), c(7L, 8L),
                 c(5L, 7L, 8L))

for (RF in names(RF.collection)) {
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)

  rand.forest <- RF.collection[[RF]]
  read.forest %<-% readForest(rand.forest, x=x, oob.importance=FALSE)

  cnt %<-% read.forest$tree.info$size.node
  precision %<-% nodePrecision(read.forest, y, cnt)

  test_that("nodePrecision works", {
    expect_equal(length(precision),
                 nrow(read.forest$tree.info))
  })

  idcnt %<-% cnt >= 1
  selected.read.forest %<-% subsetReadForest(read.forest, idcnt)
  selected.cnt %<-% cnt[idcnt]
  
  node.feature <- as(selected.read.forest$node.feature, 'dgTMatrix')
  nf.list <- split(node.feature@i + 1L, node.feature@j + 1L)

  ximp %<-% intImportance(int.eval, nf.list, precision, selected.cnt)

  test_that("intImportance works", {
    expect_true('data.frame' %in% class(ximp))
    expect_equal(names(ximp),
                 c("prev1", "prev0", "prec"))
  })

  imp.test %<-% subsetTest(int.eval, ints.sub, ximp)

  test_that("subsetTest works", {
    expect_true('data.table' %in% class(imp.test))
    expect_equal(names(imp.test),
                 c("prev.test", "prec.test"))
  })
}

