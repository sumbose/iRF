#' random forest
#'
#' @export
#'
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
#' @importFrom ranger ranger
"randomForest" <- function(x, ...) UseMethod("randomForest")

parRF <- function(x, y, xtest=NULL, ytest=NULL, ntree=500,
                  n.core=1, mtry.select.prob=rep(1, ncol(x)),
                  type='randomForest', keep.inbag=TRUE, ...) {
  
  # Wrapper function to run RF in parallel using randomForest or ranger
  if (type == 'randomForest') {
    rf <- randomForestPar(x, y, xtest, ytest, ntree, n.core, 
                          mtry.select.prob, keep.inbag, ...)
  } else if (type == 'ranger') {
    rf <- rangerPar(x, y, xtest, ytest, ntree, n.core, 
                    mtry.select.prob, keep.inbag, ...)
  } else {
    stop('type must be one of "randomForest" or "ranger"')
  }

  return(rf)
}

rangerPar <- function(x, y, xtest=NULL, ytest=NULL, ntree=500,
                      n.core=1, mtry.select.prob=rep(1, ncol(x)),
                      keep.inbag=TRUE, ...) {
  
  # Run feature weighted ranger in parallel
  mtry.select.prob <- mtry.select.prob / sum(mtry.select.prob)
  class.irf <- is.factor(y)
  if (class.irf) y <- as.numeric(y) - 1
  rf <- ranger(data=cbind(x, y), num.trees=ntree, verbose=FALSE,
               dependent.variable.name='y', classification=class.irf,
               num.threads=n.core, importance='impurity', keep.inbag=keep.inbag,
               split.select.weights=mtry.select.prob, ...)
  return(rf)
}

randomForestPar <- function(x, y, xtest=NULL, ytest=NULL, ntree=500, 
                            n.core=1, mtry.select.prob=rep(1, ncol(x)), 
                            keep.inbag=TRUE, ...) {  
  
  # Run randomForest in parallel using foreach and dorng
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)

  # Set number of trees per RF for each core
  a <- floor(ntree / n.core)
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  suppressWarnings(
    rf <- foreach(i=1:length(ntree.id), .combine=combine, 
                  .multicombine=TRUE, .packages='iRF') %dorng% {
                    randomForest(x, y, xtest, ytest,
                                 ntree=ntree.id[i],
                                 mtry.select.prob=mtry.select.prob,
                                 keep.forest=TRUE, keep.inbag=keep.inbag,
                                 ...)                         
    }
  )

  stopImplicitCluster()
  return(rf)
}

