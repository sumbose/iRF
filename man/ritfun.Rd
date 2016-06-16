\name{ritfun}
\alias{ritfun}
\title{run a basic version of random intersection tree (RIT)}
\description{
  This function runs a basic version of random intersection trees (Shah
  and Meinshausen, 2014) to
  find class-specific high-order interactions among binary features
}
\usage{
  ritfun(node_feature, tree_depth=5, n_child=2, n_ritree=500,
varnames=NULL, verbose=TRUE, node_wt=NULL)
}
\arguments{
  \item{node_feature}{a matrix of 0/1 or TRUE/FALSE, rows
    correspond to observations and columns correspond to
    features. Typically, output from \code{readForest()}}
  \item{tree_depth}{Depth of intersection trees in RIT. Default set at 5}
  \item{n_child}{numeric scalar. If set to n, RIT grows n-ary
    intersection trees. Default set at  2}
  \item{n_ritree}{Number of random intersection trees to grow}
  \item{varnames}{optionally, a vector containing feature names}
  \item{verbose}{display candidate interactions as the algorithm
    progresses}
  \item{node_wt}{a vector of non-negative weights. If specified,
    rows of the input matrix are sampled with probability proportional
    to the weights while growing RIT}
}
\value{ A character vector containing interactions, defined as prevalent
  features or feature combinations in the rows of the input
  matrix. Feature combinations are represented as feature names
  separated by "_"
}
%\details{
%}
\references{
Shah, R. D., & Meinshausen, N. (2014). Random intersection trees. The
Journal of Machine Learning Research, 15(1), 629-654.
}
%\examples{
%}
\author{Sumanta Basu \email{sumbose@berkeley.edu}}

