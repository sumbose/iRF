suite <- 'test-margin'


x <- iris[, -5]
y <- iris[, 5]
RF.collection <- make.RF.collection(x, y)

for (RF in names(RF.collection)) {
  `%<-%` <- `%<-meta.cache%`(suite, RF, TRUE)

  if (RF == 'ranger') {
    # Not implemented yet
    next
  }

  rand.forest <- RF.collection[[RF]]

  test_that("margin works for randomForest", {
    margin.obs %<-% margin(rand.forest)
    expect_equal(length(margin.obs), 150)
  })

  postscript(file='plot_margin.ps')
  plot(margin.obs)
  dev.off()
  hash.value %<-% tools::md5sum('plot_margin.ps')
}

