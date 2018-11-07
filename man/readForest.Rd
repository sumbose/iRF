\name{readForest}
\alias{readForest}
\title{Pass data through a fitted forest, record node characteristics [works
  for continuous variables]}
\description{
  Passes a feature matrix (and optionally a response vector) through a fitted
  random forest object, records size of each node. Optionally,for every node,
  returns the features used to define the rule and the data points falling in
  that node. Uses \code{foreach} function to distribute computation across
  available cores.
}
\usage{
  function(rfobj, x, 
           return.node.feature=TRUE,
           return.node.obs=FALSE, 
           varnames.group=NULL,
           get.split=FALSE,
           first.split=TRUE,
           n.core=1)
}
\arguments{
  \item{rfobj}{a fitted \code{randomForest} object with the
    \code{forest} component in it}
  \item{x}{numeric matrix with the same number of predictors used in
    \code{rfobj} fit}
  \item{return.node.feature}{if TRUE, returns a sparse matrix that encodes the
    features and inequality directions used to define the decision rule associated 
    with a node}
  \item{return.node.obs}{if TRUE, returns a sparse matrix indicating the
    observations that fall in a leaf node}
  \item{varnames.grp}{a vector of length ncol(x) indicating feature groupings}
  \item{get.split}{if TRUE, entries of \code{node.feature} indicate the
    threshold associated with a feature split}
  \item{first.split}{if TRUE, only the threshold associated with the first time
    a feature is selected will be recorded}
  \item{n.core}{the number of cores to parallelize over}
}
\value{ A list containing the following items:
  \item{tree.info}{a data frame with number of rows equal to total number of
    nodes in the forest, giving node level attributes: 
    \code{prediction} (predicted response of leaf node), 
    \code{node.idx} (the forest level node index of the leaf), 
    \code{parent} (index of parent node), 
    \code{size.node} (number of data points falling in a node), 
    \code{tree} (index of the tree in the forest in which the node lives), 
    \code{dec.purity} (if wt.pred.accuracy=TRUE, the decrease in standard 
    deviation of responses relative to the full data).}
  \item{node.feature}{if return.node.feature = TRUE, returns a sparse matrix
    with ncol(x) columns, each row corresponding to a leaf node in the
    forest. The entries indicate which features were used to define the
    decision rule associated with a node}
  \item{node.feature}{if return.node.obs = TRUE, returns a sparse matrix
    with nrow(x) columns, each row corresponding to a leaf node in the
    forest. The entries indicate which observations fall in the corresponding
    leaf node.}
}

\seealso{
 \code{\link{getTree}}
}

\examples{
  library(doMC)
  registerDoMC()
  options(cores = 7)

  n = 200; p = 250
  X = array(rnorm(n*p), c(n, p))
  Y = (X[,1]>0.35 & X[,2]>0.35)|(X[,5]>0.35 & X[,7]>0.35)
  Y = as.factor(as.numeric(Y>0))

  train.id = 1:(n/2)
  test.id = setdiff(1:n, train.id)
  
  rf <- randomForest(x=X[train.id,], y=Y[train.id], xtest=X[test.id,],
                     ytest=Y[test.id], keep.forest=TRUE)

  rforest <- readForest(rfobj=rf, x=X[train.id,], y=as.numeric(Y[train.id])-1)

  head(rforest$tree_info)

  # only consider leaf nodes
  rforest <- readForest(rfobj=rf, x=X[train.id,], y = as.numeric(Y[train.id])-1)

  # count number of leaf nodes with at least 5 observations
  sum(rforest$tree.info$size.node > 5)

  # pass test data through the forest
  rforest <- readForest(rfobj=rf, X[test.id,], y=as.numeric(Y[test.id])-1)

  # count number of leaf nodes with at least 5 observations
  sum(rforest$tree.info$size.node > 5)

}

\author{Sumanta Basu \email{sumbose@berkeley.edu}}

