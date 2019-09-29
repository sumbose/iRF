#' Iterative random forests (iRF)
#'
#' Iteratively grow feature weighted random forests and search for prevalent
#' interactions on decision paths.
#' 
#' @param x numeric feature matrix.
#' @param y response vector. If factor, classification is assumed.
#' @param xtest numeric feature matrix for test set.
#' @param ytest response vector for test set.
#' @param n.iter number of iterations to run.
#' @param ntree number of random forest trees.
#' @param mtry.select.prob feature weights for first iteration. Defaults to
#'  equal weights
#' @param iter.return which iterations should the RF be returned for.
#'  Defaults to iteration with highest OOB accuracy.
#' @param int.return which iterations should interacitons be returned for.
#' @param select.iter if TRUE, returns interactions from iteration with highest
#'  OOB accuracy.
#' @param rit.param named list specifying RIT parameters. Entries include
#'  \code{depth}: depths of RITs, \code{ntree}: number of RITs, \code{nchild}:
#'  number of child nodes for each RIT, \code{class.id}: 0-1 indicating which
#'  leaf nodes RIT should be run over, \code{min.nd}: minimum node size to run
#'  RIT over, \code{class.cut}: threshold for converting leaf nodes in
#'  regression to binary classes.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with 
#'  the same name will be treated as identical for interaction search.
#' @param n.bootstrap number of bootstrap samples to calculate stability
#'  scores.
#' @param bs.sample list of observation indices to use for bootstrap samples.
#'  If NULL, iRF will take standard bootstrap samples of observations.
#' @param weights numeric weight for each observation. Leaf nodes will be
#'  sampled for RIT with probability proprtional to the total weight of
#'  observations they contain.
#' @param signed if TRUE, signed interactions will be returned.
#' @param oob.importance if TRUE, importance measures are evaluated on OOB
#'  samples.
#' @param verbose if TRUE, display progress of iRF fit.
#' @param n.core number of cores to use. If -1, all available cores are used.
#' @param ... additional arguments passed to iRF::randomForest.
#'
#' @return A list containing the following entries:
#' \itemize{
#'    \item{rf.list}{a list of randomForest objects}
#'    \item{interaction}{a data table containing recovered interactions and
#'      importance scores}
#'    \item{selected.iter}{iterations returned by iRF}
#'    \item{weights}{feature weights used to fit each entry of rf.list}
#'  }
#'
#' @export
#'
#' @useDynLib iRF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom AUC auc roc
iRF <- function(x, y, 
                xtest=NULL, 
                ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                mtry.select.prob=rep(1, ncol(x)),
                iter.return=n.iter, 
                int.return=NULL,
                select.iter=FALSE,
                rit.param=list(depth=5, ntree=500, 
                               nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL), 
                varnames.grp=colnames(x), 
                n.bootstrap=1,
                bs.sample=NULL,
                weights=rep(1, nrow(x)),
                signed=TRUE,
                oob.importance=TRUE,
                type='randomForest',
                verbose=TRUE,
                n.core=1, 
                interactions.return=NULL,
                wt.pred.accuracy=NULL,
                ...) {
 
  # Check for depricated arguments
  if (!is.null(interactions.return)) {
    warning('interactions.return is depricated, use iter.return instead')
    iter.return <- interactions.return
    int.return <- interactions.return
    select.iter <- FALSE
  }

  if (!is.null(wt.pred.accuracy))
    warning('wt.pred.accuracy is depricated')

  # Check input attributes for correct format
  require(doRNG, quiet=TRUE)
  if (!class(x) %in% c('data.frame', 'matrix')) {
    sp.mat <- attr(class(x), 'package') == 'Matrix'
    if (is.null(sp.mat) || !sp.mat)
      stop('x must be matrix or data frame')
  }

  if (nrow(x) != length(y))
    stop('x and y must contain the same number of observations')
  if (ncol(x) < 2 && (!is.null(int.return) || select.iter))
    stop('cannot find interaction - x has less than two columns!')
  if (any(iter.return > n.iter) || any(int.return > n.iter))
    stop('selected iteration to return greater than n.iter')
  if (!is.null(varnames.grp) && length(varnames.grp) != ncol(x))
    stop('length(varnames.grp) must be equal to ncol(x)')
  if (length(mtry.select.prob) != ncol(x))
    stop('length mtry.select.prob must equal number of features')
  if (length(weights) != nrow(x))
    stop('length weights differs from # training observations')
  if (!is.null(xtest)) {
    if (ncol(xtest) != ncol(x)) 
      stop('training/test data must have same number of features')
    if (is.null(ytest))
      stop('test set responses not indicated')
    if (nrow(xtest) != length(ytest))
      stop('xtest and ytest must contain the same number of observations')
  }

  # Check all RIT and set to defaul if missing
  if (is.null(rit.param$depth)) rit.param$depth <- 5
  if (is.null(rit.param$ntree)) rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) rit.param$nchild <- 2
  if (is.null(rit.param$class.id)) rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) && is.numeric(y)) 
    rit.param$class.cut <- median(y)

  # Set variable and grouping names if not supplied
  if (is.null(colnames(x))) 
    colnames(x) <- paste0('X', 1:ncol(x))  
  if (is.null(varnames.grp)) 
    varnames.grp <- colnames(x)

  class.irf <- is.factor(y)
  imp.str <- ifelse(type == 'ranger', 'variable.importance', 'importance')
  
  # Fit a series of iteratively re-weighted RFs 
  rf.list <- list()  
  for (iter in 1:n.iter) {
    
    # Grow Random Forest on full data
    if (verbose) print(paste('iteration = ', iter))
    rf.list[[iter]] <- parRF(x, y, xtest, ytest, ntree=ntree, n.core=n.core, 
                             type=type, mtry.select.prob=mtry.select.prob, 
                             keep.inbag=oob.importance, ...)
    
    # Update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]][[imp.str]]
  }
  

  # Select iteration to return interactions based on OOB error
  if (select.iter) {
    selected.iter <- selectIter(rf.list, y=y)
    iter.return <- selected.iter 
    int.return <- selected.iter
  }

  # Generate bootstrap samples for stability analysis
  if (is.null(bs.sample) && !is.null(int.return)) 
    bs.sample <- lreplicate(n.bootstrap, bsSample(y))
    
  importance <- list()
  for (iter in int.return) {
    
    # Run gRIT across RF grown on full dataset to extract interactions.
    if (verbose) cat('finding interactions...\n')
    rit.param$ntree <- rit.param$ntree# * n.bootstrap
    ints.eval <- gRIT(rf.list[[iter]], x=x, y=y,
                      weights=weights,
                      rit.param=rit.param,
                      varnames.grp=varnames.grp,
                      signed=signed,
                      oob.importance=oob.importance,
                      n.core=n.core)

    ints.idx.eval <- ints.eval$int.idx

    # Grow RFs on BS samples to evaluate stability of recovered interactions.
    if (length(ints.eval) > 0) {
      if (verbose) cat('evaluating interactions...\n')
      if (iter == 1) rf.weight <- rep(1, ncol(x))
      if (iter > 1) rf.weight <- rf.list[[iter - 1]][[imp.str]]
      importance[[iter]] <- stabilityScore(x, y, 
                                           ntree=ntree,
                                           mtry.select.prob=rf.weight,
                                           ints.idx.eval=ints.idx.eval,
                                           rit.param=rit.param,
                                           varnames.grp=varnames.grp,
                                           bs.sample=bs.sample,
                                           weights=weights, 
                                           signed=signed,
                                           oob.importance=oob.importance,
                                           type=type,
                                           n.core=n.core, 
                                           ...)
    } else {
      importance[[iter]] <- nullReturnStab()
    }
  
  }
  
  # Combine reults for return
  out <- list()
  out$rf.list <- rf.list
  if (select.iter) out$selected.iter <- selected.iter
  if (!is.null(int.return)) out$interaction <- importance

  if (length(iter.return) == 1) {
    iter.wt <- iter.return - 1
    if (iter.return > 1) out$weights <- out$rf.list[[iter.wt]][[imp.str]] 
    out$rf.list <- out$rf.list[[iter.return]]
  }

  if (length(int.return) == 1) {
    out$interaction <- importance[[int.return]]
  }

  return(out)
}


selectIter <- function(rf.list, y) {
  # Evaluate optimal iteration based on prediction error in OOB samples.
  # For classification: accuracy. For regression: MSE.
  type <- class(rf.list[[1]])

  if (type == 'randomForest') {
    predicted <- lapply(rf.list, function(z) as.numeric(z$predicted))
  } else if (type == 'ranger') {
    predicted <- lapply(rf.list, function(z) as.numeric(z$predictions))
  } else {
    stop('rf.list must contain ranger or randomForest objects')
  }
  
  if (is.factor(y)) {
    predicted <- lapply(predicted, '-', 1)
    y <- as.numeric(y) - 1
    eFun <- function(y, yhat) sum(xor(y, yhat))
  } else {
    eFun <- function(y, yhat) mean((yhat - y) ^ 2, na.rm=TRUE)
  }
  
  error <- sapply(predicted, eFun, y=y)
  min.err <- min(error)
  id.select <- max(which(error == min.err))
  return(id.select)
}

sampleClass <- function(y, cl, n) {
  # Take a bootstrap sample from a particular class of observations
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

bsSample <- function(y) {
  # Generate outer layer bootstrap samples
  n <- length(y)
  if (is.factor(y)) {
    # Take bootstrap sample that maintains class balance of full data
    ncl <- table(y)
    class <- as.factor(names(ncl))
    sample.id <- mapply(function(cc, n) sampleClass(y, cc, n), class, ncl)
    sample.id <- c(unlist(sample.id))
  } else {
    sample.id <- sample(n, replace=TRUE)
  }
  return(sample.id)
}

