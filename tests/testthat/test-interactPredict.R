suite <- 'test-interactPredict'


x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)
int <- 'Petal.Length+_Petal.Width+'

for (RF in names(RF.collection)) {
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)

  rand.forest <- RF.collection[[RF]]
  read.forest %<-% readForest(rand.forest, x=x, oob.importance=FALSE)
  info %<-% read.forest$tree.info

  test_that('sampleTree works', {
    set.seed(42)
    leaf.id %<-% sampleTree(42, info$tree, info$size.node)
    expect_equal(length(leaf.id), 1)
  })

  test_that('interactPredict works', {
    set.seed(42)
    ip %<-% interactPredict(x, int, read.forest)
    expect_equal(length(ip), nrow(x))
  })
}

