\name{partialPlot2var}
\alias{partialPlot2var}
\title{creates 3-dimensional surface plots for bivariate partial dependence functions}
\description{
  Given two feature vectors and a response vector, produces surface maps
  of partial dependence function. Uses R package \code{rgl}
}
\usage{
  partialDep2var(x1, x2, y, gridlength=NULL, x1_grid=NULL, x2_grid=NULL,
x1lab='v1', x2lab='v2', ylab=NA, range.color=NULL, col.palette=c('blue',
'yellow'), plot_quantile_scale = TRUE, plot.colorbar=TRUE, ...)
}
\arguments{
  \item{x1, x2, y}{numeric vectors containing the first, second and the response variable}
  \item{gridlength}{If either x1_grid or x2_grid is missin, a scalar
    indicating length of grids against which to make the plot}
  \item{x1_grid, x2_gris}{grids of points between which the surface map will be constant}
  \item{x1lab, x2lab, ylab}{Labels of x, y and z-axis respectively}
  \item{range.color}{range of values which should be colored. input in
    the form (min, max). If not provided, the program automatically
    sets it to the range of the partial dependence function calculated
    on grid points}
  \item{col.palette}{Color palette, typically a vector of length 2
    indicating the range of colors}
  \item{plot_quantile_scale}{If TRUE, plots against quantiles of x1 and
    x2 instead of their raw values}
  \item{plot.colorbar}{produce a separate colorbar for the partial dependence plot?}
}
\value{ 
A RGL plot of type persp3d() of the package \code{rgl}
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{rgl::persp3d}}
}
\examples{

}
\author{Sumanta Basu \email{sumbose@berkeley.edu}}
%\keyword{favourite keyword}
