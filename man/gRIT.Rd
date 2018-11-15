\name{gRIT}
\alias{gRIT}
\title{Generalized Random Intersection Trees}
\description{
  Searches for prevalent interactions on random forest decision paths.
}
\usage{
  gRIT(x, y,
       rand.forest=NULL,
       read.forest=NULL,
       weights=rep(1, nrow(x)),
       varnames.grp=colnames(x),
       rit.param=list(depth=5, ntree=500,
                      nchild=2, class.id=1,
                      min.nd=1, class.cut=NULL),
       signed=TRUE,
       ints.full=NULL,
       n.core=1)
}
\arguments{
  \item{x}{numeric matrix of predictors}
  \item{y}{response vector. If factor classification is assumed.}
  \item{rand.forest}{An object of class randomForest. Required if 
    read.forest = NULL}
  \item{read.forest}{Ouput from the readForest function. Required if 
    rand.forest = NULL}
  \item{weights}{0, 1 indicator of weither an observation should be included
    when evaluating interactions}
  \item{varnames.grp}{If features can be grouped based on their
    demographics or correlation patterns, use the group of features or
    ``hyper-feature''s to conduct random intersection trees}
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
  \item{signed}{If TRUE, signed interactions will be returned, indicating 
    the direction of the inequality associated with the first time a feature 
    is selected for splitting}
  \item{ints.full}{A list of interactions to evaluate. If specified, prevalence
    will be evaluated for these interactions instead of those recovered by RIT.}
  \item{n.core}{number of cores across which tree growing should be
    distributed. If -1, all available cores will be used.}

}
\value{ A list containing the following items:
  \item{int}{A character vector indicating recovered interactions. Naming
    convention for an interaction is: <variable1(sign)>_ <variable2(sign)>_...}
  \item{prevalence}{A data table indicating: (1) the recovered interactions (2)
  the prevalence of each interaction in class-1 leaf nodes (3) the prevalence of
  each interaction in class-0 leaf nodes (4) the proportion of class-1
  observations in leaf nodes containing the interactions. If regression, class-1
  is determined by rit.param$class.cut}
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
