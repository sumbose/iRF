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
  iRF(x, y, xtest=NULL, ytest=NULL, n_iter=5, ntree=500, n_core=1,
      mtry_select_prob = rep(1/ncol(x), ncol(x)),
      keep_impvar_quantile=NULL, find_interaction=FALSE,
      cutoff_nodesize_prop = 0.1, cutoff_unimp_feature = 0,
      class_id=1, rit_param=c(5, 1000, 2), varnames_grp=NULL,
      n_bootstrap=30, verbose=TRUE, ...
     )
}
\arguments{
  \item{x, xtest}{numeric matrices of predictors}
  \item{y, ytest}{factor with two levels: 0, 1}
  \item{n_iter}{number of weighted random forest fits}
  \item{ntree}{number of trees to grow in each iteration}
  \item{n_core}{number of cores across which tree growing should be
    distributed}
  \item{mtry_select_prob}{initial weights specified for first random
    forest fit, defaults to equal weights}
  \item{keep_impvar_quantile}{a nonnegative fraction q. If provided, all
  the variables with Gini importance in the top 100*q percentile are
  retained during random splitting variable selection in the next iteration}
  \item{find_interaction}{find stable feature interaction using random
    intersection trees?}
  \item{cutoff_nodesize_prop}{a non-negative fraction p. Only leaf nodes
    containing at least a fraction p of the observations from clas 
    \code{class_id} are used to find feature interactions}
  \item{cutoff_unimp_feature}{a non-negative fraction r. If provided,
    only features with Gini importance score in the top 100*(1-r)
    percentile are used to find feature interactions}
  \item{class_id}{which class of observations will be used to find
    class-specific interaction? Choose between 0 or 1. Default is set to 1.}
  \item{rit_param}{a numeric triplet, containing (in this order) (a)
    depth of random intersection trees, (b) number of random
    intersection trees, (c) number of children in each split of a
    random intersection tree}
  \item{varnames_grp}{If features can be grouped based on their
    demographics or correlation patterns, use the group of features or
    ``hyper-feature''s to conduct random intersection trees}
  \item{n_bootstrap}{Number of bootstraps replicates used to calculate
    stability scores of interactiosn obtained by RIT}
  \item{verbose}{Display progress messages and intermediate outputs on screen?}
  \item{...}{additional arguments passed to iRF::randomForest}
}
\value{ A list containing the following items:
  \item{rf_list}{A list of n_iter objects of the class randomForest}
  \item{interaction}{A list of length n_iter. Each element of the list
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
