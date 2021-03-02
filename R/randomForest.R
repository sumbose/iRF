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

parRF <- function(x, y, 
                  xtest=NULL, 
                  ytest=NULL, 
                  ntree=500,
                  n.core=1, 
                  mtry=floor(sqrt(ncol(x))), 
                  mtry.select.prob=rep(1, ncol(x)),
                  type='randomForest', 
                  keep.inbag=TRUE, ...) {

  mtry.select.prob <- pmax(mtry.select.prob, 0)  
  mtry.select.prob <- mtry.select.prob / sum(mtry.select.prob)

  # Wrapper function to run RF in parallel using randomForest or ranger
  if (type == 'randomForest') {
    rf <- randomForestPar(x, y, xtest, ytest, ntree, n.core, 
                          mtry, mtry.select.prob, keep.inbag, ...)
  } else if (type == 'ranger') {
    rf <- rangerPar(x, y, xtest, ytest, ntree, n.core, 
                    mtry, mtry.select.prob, keep.inbag, ...)
  } else {
    stop('type must be one of "randomForest" or "ranger"')
  }

  return(rf)
}

rangerPar <- function(x, y, 
                      xtest=NULL, 
                      ytest=NULL, 
                      ntree=500,
                      n.core=1, 
                      mtry=floor(sqrt(ncol(x))), 
                      mtry.select.prob=rep(1, ncol(x)),
                      keep.inbag=TRUE, 
                      ...) {
  
  # Threshold feature weights at 0
  mtry.select.prob[mtry.select.prob < 0] <- 0
  mtry.select.prob <- mtry.select.prob / sum(mtry.select.prob)
  
  # Format response for classification/regression
  class.irf <- is.factor(y)
  if (class.irf) y <- as.numeric(y) - 1

  # Check if running local importance and set feature importance mode
  dot.args <- list(...)

  if('importance' %in% names(dot.args)) {
      importance <- dot.args[['importance']]
  } else {
      importance <- 'impurity'
  }

  if ('local.importance' %in% names(dot.args)) {
      if (dot.args[['local.importance']]) {
          importance <- 'permutation'
      }
  }

  dot.args[['importance']] <- NULL
      
  # Run ranger
  rf <- ranger(data=cbind(x, y), 
               num.trees=ntree, 
               verbose=FALSE,
               dependent.variable.name='y', 
               classification=class.irf,
               num.threads=n.core, 
               keep.inbag=keep.inbag, 
               mtry=mtry, 
               split.select.weights=mtry.select.prob, 
               importance=importance,
               ...)
  return(rf)
}

randomForestPar <- function(x, y, 
                            xtest=NULL, 
                            ytest=NULL, 
                            ntree=500, 
                            n.core=1, 
                            mtry=floor(sqrt(ncol(x))), 
                            mtry.select.prob=rep(1, ncol(x)), 
                            keep.inbag=TRUE, 
                            ...) {  
  
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
                                 mtry=mtry,
                                 mtry.select.prob=mtry.select.prob,
                                 keep.forest=TRUE, keep.inbag=keep.inbag,
                                 ...)                         
    }
  )

  stopImplicitCluster()
  return(rf)
}

