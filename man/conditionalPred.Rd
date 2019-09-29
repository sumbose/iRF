\name{conditionalPred}
\alias{conditionalPred}
\title{Evaluates interaction importance using conditional prediction}
\description{
  Computes the prediction accuracy of a fitted Random Forest using only the leaf
  nodes for which the specified interaction falls on the decision path. For 
  classification, accuracy is measured by AUROC. For regression, accuracy is
  measured by decrease in variance.
}
\usage{
  conditionalPred(rfobj, rd.forest, x, y, ints,
                  varnames.group=NULL,
                  n.cores=1,
  )
}
\arguments{
  \item{rfobj}{Fitted randomForest object}
  \item{rd.forest}{list as returned by \code{readForest}}
  \item{x}{numeric matrix of predictors}
  \item{y}{response vector}
  \item{ints}{a character vector specifying interactions, features separated by
    '_', as returned by \code{iRF}}
  \item{n.perms}{number of times to permute data matrix}
  \item{varnames.grp}{If features can be grouped based on their
    demographics or correlation patterns, use the group of features or
    ``hyper-feature''s to conduct random intersection trees}
  \item{n.cores}{number of cores to parallelize over}
}
\value{A numeric vector of the same length as \code{ints} giving the conditional
       prediction accuracy for each interaction}
%\details{
%}
%\references{
%}
\seealso{
  \code{randomForest}
}

%\author{Author Name \email{abcd\_xyz@domain.com}}
%\keyword{favourite keyword}
