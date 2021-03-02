require(ranger)
require(iRF)

x <- iris[, -5]
y <- iris[, 5]
class.irf <- is.factor(y)
if (class.irf) y <- as.numeric(y) - 1
rand.forest <- ranger(data=cbind(x, y),
                      dependent.variable.name='y',
                      classification=class.irf)

require(microbenchmark)
microbenchmark(readForest(rand.forest, x=x, return.node.obs=FALSE), times=5L)

require(profvis)
p <- profvis(readForest(rand.forest, x=x, return.node.obs=FALSE))
htmlwidgets::saveWidget(p, "profile.html", selfcontained=FALSE)

