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
stabilityScore <- function(x, y, 
                           ntree=500,
                           mtry.select.prob=rep(1, ncol(x)), 
                           ints.eval=NULL,
                           rit.param=list(depth=5, ntree=500,
                                          nchild=2, class.id=1,
                                          min.nd=1, class.cut=NULL), 
                           varnames.grp=NULL,
                           n.bootstrap=1,
                           bs.sample=NULL,
                           weights=rep(1, nrow(x)),
                           signed=TRUE, 
                           n.core=1,
                           ...) {
  
  # Set feature names for grouping interactions
  varnames.grp <- groupVars(varnames.grp, x)

  # Generate bootstrap samples for stability analysis
  if (is.null(bs.sample)) bs.sample <- lreplicate(n.bootstrap, bsSample(y))

  out <- list()
  for (i in 1:length(bs.sample)) {

    sample.id <- bs.sample[[i]]
    out[[i]] <- bsgRIT(x, y, mtry.select.prob, sample.id, ints.eval=ints.eval, 
                       ntree=ntree, weights=weights, rit.param=rit.param,
                       varnames.grp=varnames.grp, signed=signed, n.core=n.core,
                       ...)

  }

  # Summarize stability and importance metrics across bootstrap replicates
  out <- summarizeInteract(out)
  return(out)
}


bsgRIT <- function(x, y, mtry.select.prob, sample.id, ints.eval, weights, ntree,
                   varnames.grp, rit.param, signed, n.core, ...) {
  
  # Fit random forest on bootstrap sample
  rf <- parRF(x[sample.id,], y[sample.id], ntree=ntree, n.core=n.core, 
              mtry.select.prob=mtry.select.prob, ...)

  # Run generalized RIT on rf.b to learn interactions
  ints <- gRIT(rand.forest=rf, x=x, y=y,
               weights=weights,
               varnames.grp=varnames.grp,
               rit.param=rit.param,
               signed=signed,
               ints.eval=ints.eval,
               n.core=n.core)

  return(ints)
}

summarizeInteract <- function(x) {
  # Summarize interaction importance metrics across bootstrap samples 
  save(file='~/interact_temp.Rdata', x)
  n.bootstrap <- length(x)
  x <- rbindlist(x)

  if (nrow(x) > 0) {
    imp <- mutate(x, diff=(prev1-prev0)) %>%
      group_by(int) %>%
      summarize(prevalence.diff=mean(diff),
                sta.diff=mean(diff >= 0),
                independence=mean(prev.test),
                sta.independence=mean(prev.test >= 0),
                precision=mean(prec),
                sta.precision=mean(prec.test >= 0),
                stability=mean(recovered)) %>%
      arrange(desc(prevalence.diff))
  } else {
    imp <- data.table(prevalence.diff=numeric(0),
                      sta.diff=numeric(0),
                      independence=numeric(0),
                      sta.independence=numeric(0),
                      precision=numeric(0),
                      sat.precision=numeric(0),
                      stability=numeric(0))
  }

  return(data.table(imp))
}

