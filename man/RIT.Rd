\name{RIT}
\alias{RIT}
\title{Random Intersection Trees}
\description{
  Function to perform random intersection trees. When two binary data matrices \code{z} (class 1) and \code{z0} (class 0) are supplied, it searches for interactions. More precisely, since the data matrices are binary, each row of each matrix can be represented by the set of column indices with non-zero entries. The function searches for sets (interactions) that are more prevalent in class 1 than class 0, and then sets that are more prevalent in class 0 than class 1. When given a single binary matrix \code{z} with the argument \code{z0} omitted, the function simply finds sets with high prevalence. Prevalences of interactions returned are estimated using min-wise hashing.
}
\usage{
RIT(z, z0, weights=rep(1, nrow(z)), branch = 5, depth = 10L, n_trees = 100L,
  theta0 = 0.5, theta1 = theta0, min_inter_sz = 2L,
  L = 100L, n_cores = 1L, output_list = FALSE)
}
\arguments{
  \item{z}{data matrix where each row corresponds to an observation and columns correspond to variables. Can be in sparse matrix format (inherit from class "sparseMatrix" in the \pkg{Matrix} package).
  }
  \item{z0}{optional second data matrix with the same number of columns as \code{z}.
  }
  \item{weights}{weighting vector specifying the sampling probability for each
    observation in z}
  \item{branch}{average number of branches to use when creating each tree.
  }
  \item{depth}{maximum depth of trees.
  }
  \item{n_trees}{number of trees to be constructed.
  }
  \item{theta0}{when searching for sets of variables that are more prevalent in class 1 than class 0, the maximum threshold for prevalence in class 0.
  }
  \item{theta1}{as above but with class 1 and class 0 interchanged.
  }
  \item{min_inter_sz}{minimum size of the interactions to be returned
  }
  \item{L}{number of rows of the min-wise hash matrix used to estimate prevalences. A larger value will result in more accurate estimates, but computation time will increase linearly with \code{L}.
  }
  \item{n_cores}{number of cores for parallel processing. Only used when openMP is installed.
  }
  \item{output_list}{if \code{FALSE} returns each interaction set as a string with variable indices separated by spaces. If \code{TRUE} returns each interaction set as an integer vector.
  }
}
\details{
  There are two tasks which can be performed with this function depending on whether or not \code{z0} is supplied (note \code{z} must always be supplied).

  1. If \code{z0} is omitted, the function finds prevalent sets in \code{z} and \code{theta0} and \code{theta1} are ignored.

  2. If \code{z0} is supplied, it searches for sets that are prevalent in \code{z} but have prevalence at most \code{theta0} in \code{z0}. Next sets that are prevalent in \code{z0} but have prevalence in \code{z} at most \code{theta1} are found.
}
\value{
  If \code{output_list} is \code{FALSE} (the default), the output is either a data frame (if \code{z0} is omitted) or list of two data frames (if \code{z0} is supplied). The data frames have first column a character vector of interaction sets with the variables in the sets separated by spaces, and second column the estimated prevalences. When \code{z0} is supplied, the interactions in the first component of the list named \code{Class1} are those which are prevalent in \code{z} and their prevalences in \code{z} are reported. The second component named named \code{Class0} contains those interactions prevalent in \code{z0} and their prevalences in \code{z0}.

  When \code{output_list} is \code{TRUE}, each interaction is reported as an integer vector and so the collection of interactions is a list of such vectors.
}
\references{
  Shah, R. D. and Meinshausen, N. (2014) Random Intersection Trees. \emph{Journal of Machine Learning Research}, \bold{15}, 629--654.
}
\author{
  Hyun Jik Kim, Rajen D. Shah with slight modifications by Karl Kumbier
}
\examples{
## Generate two binary matrices
z <- matrix(rbinom(250*500, 1, 0.3), 250, 500)
z0 <- matrix(rbinom(250*500, 1, 0.3), 250, 500)

## Make the first and second cols of z identical
## so the set 1, 2 has prevalence roughly 0.3 compared
## to roughly 0.09 for any other pair of columns
z[, 1] <- z[, 2]

## Similarly for z0
z0[, 3] <- z0[, 4]

## Market basket analysis
out1 <- RIT(z)
out1[1:5, ]

## Finding interactions
out2 <- RIT(z, z0)
out2$Class1[1:5, ]
out2$Class0[1:5, ]

## Can also perform the above using sparse matrices
if (require(Matrix)) {
  S <- Matrix(z, sparse=TRUE)
  S0 <- Matrix(z0, sparse=TRUE)
  out3 <- RIT(S, S0)
}
}
