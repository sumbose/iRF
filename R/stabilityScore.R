#' Stability score
#'
#' Run outer layer bootstrap stability analysis
#'
#' @param x numeric feature matrix
#' @param y response vector. If factor, classification is assumed.
#' @param ntree number of random forest trees.
#' @param mtry.select.prob feature weights for first iteration. Defaults to
#'  equal weights
#' @param ints.eval interactions to evaluate. If specified, importance metrics
#'  will be evaluated for these interactions instead of those recovered by RIT.
#' @param ints.idx.eval like \code{ints.eval}, but specifies the indice of the
#'  interactions instead of their names. Intended for internal use only.
#' @param rit.param named list specifying RIT parameters. Entries include
#'  \code{depth}: depths of RITs, \code{ntree}: number of RITs, \code{nchild}:
#'  number of child nodes for each RIT, \code{class.id}: 0-1 indicating which
#'  leaf nodes RIT should be run over, \code{min.nd}: minimum node size to run
#'  RIT over, \code{class.cut}: threshold for converting leaf nodes in
#'  regression to binary classes.
#' @param n.bootstrap number of bootstrap samples to calculate stability scores.
#' @param bs.sample list of observation indices to use for bootstrap samples. If
#'  NULL, iRF will take standard bootstrap samples of observations.
#' @param weights numeric weight for each observation. Leaf nodes will be
#'  sampled for RIT with probability proprtional to the total weight of
#'  observations they contain.
#' @param signed if TRUE, signed interactions will be returned
#' @param oob.importance if TRUE, importance measures are evaluated on OOB
#'  samples.
#' @param n.core number of cores to use. If -1, all available cores are used.
#' @param ... additional arguments passed to iRF::randomForest.
#'
#' @return a data table containing the recovered interactions and importance
#'  scores.
#'
#' @export
#'
#' @importFrom dplyr group_by summarize arrange desc '%>%'
#' @importFrom data.table data.table
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
stabilityScore <- function(x, y,
                           ntree=500,
                           mtry.select.prob=rep(1, ncol(x)), 
                           ints.idx.eval=NULL,
                           ints.eval=NULL,
                           rit.param=list(depth=5, ntree=500,
                                          nchild=2, class.id=1,
                                          min.nd=1, class.cut=NULL), 
                           varnames.grp=NULL,
                           n.bootstrap=1,
                           bs.sample=NULL,
                           weights=rep(1, nrow(x)),
                           signed=TRUE,
                           oob.importance=TRUE,
                           type='randomForest',
                           n.core=1,
                           ...) {
  
  # Check for valid input parameters
  if (!class(x) %in% c('data.frame', 'matrix')) {
    sp.mat <- attr(class(x), 'package') == 'Matrix'
    if (!is.null(sp.mat)) {
      if (!sp.mat) stop('x must be matrix or data frame')
    } else {
      stop('x must be matrix or data frame')
    }
  }
  if (nrow(x) != length(y))
    stop('x and y must contain the same number of observations')
  if (length(mtry.select.prob) != ncol(x))
    stop('length mtry.select.prob must equal number of features')
  if (length(weights) != nrow(x))
    stop('length weights differs from # training observations')


  # Set feature names for grouping interactions
  varnames.grp <- groupVars(varnames.grp, x)

  # Register cores for parallelization
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)

  # Generate bootstrap samples for stability analysis
  if (is.null(bs.sample))
      bs.sample <- lreplicate(n.bootstrap, bsSample(y))

  suppressWarnings(
  out <- foreach(sample.id=bs.sample) %dorng% {
    # Use only 1 core for each bsgRIT, as the loop is already parallelized
    # Note that reproducibility is guaranteed even with ``n.core=n.core'',
    # so feel free to use more cores if you benchmark says otherwise.
    bsgRIT(x, y, mtry.select.prob, sample.id,
           ints.idx.eval=ints.idx.eval,
           ints.eval=ints.eval, ntree=ntree, weights=weights,
           rit.param=rit.param, varnames.grp=varnames.grp,
           signed=signed, oob.importance=oob.importance,
           n.core=1L, type=type, ...)
  })

  stopImplicitCluster()

  # Summarize stability and importance metrics across bootstrap replicates
  out <- summarizeInteract(out)
  return(out)
}


bsgRIT <- function(x, y, mtry.select.prob, sample.id, ints.idx.eval,
                   ints.eval, weights, ntree, varnames.grp, rit.param,
                   signed, oob.importance, type, n.core, ...) {

  # Remove replicates in bs sample for OOB importance
  if (oob.importance) sample.id <- unique(sample.id)
 
  # Generate bootstrap sample for stability analysis
  x <- x[sample.id,]
  y <- y[sample.id]

  # Fit random forest on bootstrap sample
  rf <- parRF(x, y, ntree=ntree, n.core=n.core, 
              mtry.select.prob=mtry.select.prob, type=type, 
              keep.inbag=oob.importance, ...)
  
  # Run generalized RIT on rf.b to learn interactions
  ints <- gRIT(rand.forest=rf, x=x, y=y,
               weights=weights[sample.id],
               varnames.grp=varnames.grp,
               rit.param=rit.param,
               signed=signed,
               oob.importance=oob.importance,
               ints.idx.eval=ints.idx.eval,
               ints.eval=ints.eval,
               n.core=n.core)

  return(ints)
}

summarizeInteract <- function(x) {
  # Summarize interaction importance metrics across bootstrap samples 
  n.bootstrap <- length(x)
  x <- rbindlist(x)

  if (nrow(x) > 0) {
    imp <- mutate(x, diff=(prev1-prev0)) %>%
      group_by(int) %>%
      summarize(prevalence=mean(prev1),
                precision=mean(prec),
                cpe=mean(diff),
                sta.cpe=mean(diff > 0),
                fsd=mean(prev.test),
                sta.fsd=mean(prev.test > 0),
                mip=mean(prec.test),
                sta.mip=mean(prec.test > 0),
                stability=mean(recovered)) %>%
      arrange(desc(cpe))
  } else {
    nullReturnStab()
  }

  return(data.table(imp))
}

nullReturnStab <- function() {
  # Returns empty data table of scored interactions
  out <- data.table(prevalence=numeric(0),
                    precision=numeric(0),
                    cpe=numeric(0),
                    sta.cpe=numeric(0),
                    fsd=numeric(0),
                    sta.fsd=numeric(0),
                    mip=numeric(0),
                    sta.mip=numeric(0),
                    stability=numeric(0))
  return(out)
}
