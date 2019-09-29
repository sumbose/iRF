suite <- 'test-readForest'

x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)


for (RF in names(RF.collection)) {
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)
  rand.forest <- RF.collection[[RF]]
  read.forest %<-% readForest(rand.forest, x=x, oob.importance=FALSE)
}


test_that('readForest works for randomForest', {
  rand.forest <- RF.collection[['randomForest']]
  read.forest <- readForest(rand.forest, x=iris[, -5], oob.importance=FALSE)

  countLeaf <- function(k)
      sum(randomForest::getTree(rand.forest, k)[, 'status'] == -1)
  nleaves <- sum(sapply(1:rand.forest$ntree, countLeaf))

  expect_true('data.table' %in%
              class(read.forest$tree.info))
  expect_equal(nrow(read.forest$tree.info),
               nleaves)
  expect_equal(names(read.forest$tree.info),
               c('prediction', 'node.idx', 'parent', 'tree', 'size.node'))

  expect_true('dgCMatrix' %in%
              class(read.forest$node.feature))
  expect_equal(nrow(read.forest$node.feature),
               nleaves)
  expect_equal(ncol(read.forest$node.feature),
               2 * length(rand.forest$importance))

  expect_true('ngCMatrix' %in%
              class(read.forest$node.obs))
  expect_equal(nrow(read.forest$node.obs),
               length(rand.forest$predicted))
  expect_equal(ncol(read.forest$node.obs),
               nleaves)
  expect_equal(rowSums(read.forest$node.obs),
               rep(rand.forest$ntree, length(rand.forest$predicted)))
})


test_that('readForest works for ranger', {
  rand.forest <- RF.collection[['ranger']]
  read.forest <- readForest(rand.forest, x=iris[, -5], oob.importance=FALSE)

  countLeaf <- function(k)
      sum(ranger::treeInfo(rand.forest, k)$terminal)
  nleaves <- sum(sapply(1:rand.forest$num.trees, countLeaf))

  expect_true('data.table' %in%
              class(read.forest$tree.info))
  expect_equal(nrow(read.forest$tree.info),
               nleaves)
  expect_equal(names(read.forest$tree.info),
               c('prediction', 'node.idx', 'parent', 'tree', 'size.node'))

  expect_true('dgCMatrix' %in%
              class(read.forest$node.feature))
  expect_equal(nrow(read.forest$node.feature),
               nleaves)
  expect_equal(ncol(read.forest$node.feature),
               2 * rand.forest$num.independent.variables)

  expect_true('ngCMatrix' %in%
              class(read.forest$node.obs))
  expect_equal(nrow(read.forest$node.obs),
               rand.forest$num.samples)
  expect_equal(ncol(read.forest$node.obs),
               nleaves)
  expect_equal(rowSums(read.forest$node.obs),
               rep(rand.forest$num.trees, rand.forest$num.samples))
})

