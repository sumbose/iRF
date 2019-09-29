\name{permImportance}
\alias{permImportance}
\title{Evaluates interaction importance thhrough permutation}
\description{
  Computes the prediction accuracy of a fitted Random Forest evaluated on data
  for which all columns have been permuted except for the specified interaction.
  For classification, accuracy is measured by AUROC. For regression, accuracy is
  measured by decrease in variance.
}
\usage{
  permImportance(rfobj, x, y, ints,
                 n.perms=3,
                 varnames.group=NULL,
                 n.cores=1,
  )
}
\arguments{
  \item{rfobj}{Fitted randomForest object}
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
\value{A numeric vector of the same length as \code{ints} giving the prediction
       accuracy with all other features permuted}
%\details{
%}
%\references{
%}
\seealso{
  \code{randomForest}
}

%\author{Author Name \email{abcd\_xyz@domain.com}}
%\keyword{favourite keyword}
