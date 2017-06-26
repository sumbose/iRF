\name{iRF}
\alias{iRF}
\title{iteratively grows weighted random forests, finds stable feature interactions}
\description{
  Using repeated calls to \code{iRF::randomForest}, this function
  iteratively grows weighted ensembles of decision trees. Optionally,
  for every iteration, returns stable feature interactions by analyzing feature
  usage on decision paths of large leaf nodes. Implemented only for
  binary classification with numeric predictors and response taking values in {0,1}.
}
\usage{
  iRF(x, y, xtest=NULL, ytest=NULL, 
      n.iter=5, 
      ntree=500, 
      n.core=1,
      mtry.select.prob = rep(1/ncol(x), ncol(x)),
      keep.impvar.quantile=NULL, 
      interactions.return=NULL,
      wt.pred.accuracy=FALSE, 
      cutoff.unimp.feature = 0,
      rit.param=list(depth=5, ntree=100, nchild=2, 
                     class.cut=NULL, class.id=1), 
      varnames.grp=NULL,
      n.bootstrap=30,
      bootstrap.forest=TRUE, 
      verbose=TRUE, ...
     )
}
\arguments{
  \item{x, xtest}{numeric matrices of predictors}
  \item{y, ytest}{factor with two levels: 0, 1}
  \item{n.iter}{number of weighted random forest fits}
  \item{ntree}{number of trees to grow in each iteration}
  \item{n.core}{number of cores across which tree growing should be
    distributed}
  \item{mtry.select.prob}{initial weights specified for first random
    forest fit, defaults to equal weights}
  \item{keep.impvar.quantile}{a nonnegative fraction q. If provided, all
  the variables with Gini importance in the top 100*q percentile are
  retained during random splitting variable selection in the next iteration}
  \item{interactions.return}{a numeric vector specifying which iterations to 
    calculate interactions for. Note: interaction computation is
    computationally intensive particularly when \code{n.bootstrap} is large.}
  \item{wt.pred.accuracy}{Should leaf nodes be sampled proportional to both size
    and decrease in variabiliy of responses?}
  \item{cutoff.unimp.feature}{a non-negative fraction r. If provided,
    only features with Gini importance score in the top 100*(1-r)
    percentile are used to find feature interactions}
  \item{class.id}{which class of observations will be used to find
    class-specific interaction? Choose between 0 or 1. Default is set to 1.}
  \item{rit.param}{a named list, containing entries to specify 
    \code{depth}: depth of random intersection trees
    \code{ntree} number of random intersection trees
    \code{nchild}: number of children in each split of a random intersection tree
    \code{class.cut}: threshold to determine leaf nodes that are used to find
    interactions. Any leaf nodes with prediction greater than
    specified threshold will be used as input to RIT.  If \code{NULL}, all leaf 
    nodes from regression iRF will be used as input to RIT. Ignored if
    classification forest}
    \code{class.id}{which class of observations will be used to find
    class-specific interaction? Choose between 0 or 1. Default is set to 1. 
    Ignored if regression forest.}
  \item{varnames.grp}{If features can be grouped based on their
    demographics or correlation patterns, use the group of features or
    ``hyper-feature''s to conduct random intersection trees}
  \item{n.bootstrap}{Number of bootstraps replicates used to calculate
    stability scores of interactiosn obtained by RIT}
  \item{bootstrap.forest}{Should a new Random Forest be constructed for each
    bootstrap sample to evaluate stability? Setting to FALSE results in faster
    runtime.}
  \item{verbose}{Display progress messages and intermediate outputs on screen?}
  \item{...}{additional arguments passed to iRF::randomForest}
}
\value{ A list containing the following items:
  \item{rf.list}{A list of n.iter objects of the class randomForest}
  \item{interaction}{A list of length n.iter. Each element of the list
    contains a named numeric vector of stability scores, where the names
    are candidate interactions (feature names separated by "_"), defined
  as frequently appearing features and feature combinations on the
  decision paths of large leaf nodes}
}
%\details{
%}
%\references{
%}
\seealso{
  \code{randomForest}, \code{readForest}
}

%\author{Author Name \email{abcd\_xyz@domain.com}}
%\keyword{favourite keyword}
