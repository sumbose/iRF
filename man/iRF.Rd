\name{iRF}
\alias{iRF}
\title{iteratively grows weighted random forests, finds stable feature interactions}
\description{
  Using repeated calls to \code{iRF::randomForest}, this function
  iteratively grows weighted ensembles of decision trees. Optionally,
  for every iteration, returns stable feature interactions by analyzing feature
  usage on decision paths of large leaf nodes.
}
\usage{
  iRF(x, y, xtest=NULL, ytest=NULL, 
      n.iter=5, 
      ntree=500, 
      n.core=1,
      mtry.select.prob = rep(1, ncol(x)),
      interactions.return=NULL,
      rit.param=list(depth=5, ntree=500, 
                     nchild=2, class.id=1, 
                     min.nd=1, class.cut=NULL), 
      varnames.grp=colnames(x),
      n.bootstrap=1,
      select.iter=TRUE,
      signed=TRUE,
      verbose=TRUE, ...
     )
}
\arguments{
  \item{x, xtest}{numeric matrices of predictors}
  \item{y, ytest}{response vector. If factor classification is assumed.}
  \item{n.iter}{number of weighted random forest fits}
  \item{ntree}{number of trees to grow in each iteration}
  \item{n.core}{number of cores across which tree growing should be
    distributed. If -1, all available cores will be used.}
  \item{mtry.select.prob}{initial weights specified for first random
    forest fit, defaults to equal weights}
  \item{interactions.return}{a numeric vector specifying which iterations to 
    calculate interactions for. Note: interaction computation can be
    computationally intensive, particularly when \code{n.bootstrap} is large.}
  \item{rit.param}{a named list, containing entries to specify 
    \code{depth}: depth of random intersection trees
    \code{ntree}: number of random intersection trees
    \code{nchild}: number of children in each split of a random intersection tree
    \code{class.id}: which class of observations will be used to find
    class-specific interaction? Choose between 0 or 1. Default is set to 1. 
    Ignored if regression forest.
    \code{min.nd}: minimum leaf node size to run RIT over.
    \code{class.cut}: threshold to binarize leaf nodes for RIT. Any leaf nodes
    with prediction greater than specified threshold will be class-1 for RIT.
    Ignored if classification forest }
  \item{varnames.grp}{If features can be grouped based on their
    demographics or correlation patterns, use the group of features or
    ``hyper-feature''s to conduct random intersection trees}
  \item{n.bootstrap}{Number of bootstrap replicates used to calculate
    stability scores of interactiosn obtained by RIT}
  \item{select.iter}{If TRUE, interactions will be computed for iteration
    selected by OOB prediction accuracy}
  \item{signed}{If TRUE, signed interactions will be returned, indicating 
    the direction of the inequality associated with the first time a feature 
    is selected for splitting}
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
  \item{prevalence}{(if get.prevalence = TRUE) evaluate the prevalance of each
    recovered interaction along decision paths.}
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
