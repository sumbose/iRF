\name{readForest}
\alias{readForest}
\title{Pass data through a fitted forest, record node characteristics [works
  for binary classification with continuous variables]}
\description{
  Passes a feature matrix (and optionally a label vector) through a fitted random forest
  object, records size (and Gini impurity) of each node. Optionally,
  for every node, returns the features used to define the rule and the data
  points falling in that node. Uses \code{foreach} function to
  distribute computation across available cores.
}
\usage{
  function(rfobj, X, Y=NULL, return_node_feature=TRUE,
return_node_data=TRUE, leaf_node_only=TRUE, verbose=TRUE)
}
\arguments{
  \item{rfobj}{a fitted \code{randomForest} object with the
    \code{forest} component in it }
  \item{X}{numeric matrix with the same number of predictors used in
    \code{rfobj} fit}
  \item{Y}{Optionally, a numeric vector containing only 0 and 1}
  \item{return_node_feature}{if TRUE, returns a matrix containing features
    used to define the decision rule associated with a node}
  \item{return_node_data}{if TRUE, returns composition of nodes, i.e.,
    data points in X which fell in a node}
  \item{leaf_node_only}{if TRUE, only retain information on leaf nodes
    in the output}
  \item{verbose}{if TRUE, print progress of the algorithm on the screen}
}
\value{ A list containing the following items:
  \item{tree_info}{a data frame with number of rows equal to total
    number of nodes in the forest. The first six columns are as in 
    the output of \code{\link{getTree}}: \code{left daughter}, \code{right daughter},
    \code{split var}, \code{split point}, \code{status},
    \code{prediction}.

    If Y is not provided, the rest of the columns are:
    \code{parent} (parent of every node), \code{size_node} (number of data
    points falling in a node), \code{depth} (depth of the node in a
    tree) and \code{tree_id} (index of the tree in the forest in which the node lives)

    If Y is provided, contains three additional columns:
    \code{size_0} (number of observations of class \code{0} falling in a
      node), \code{Gini} (Gini impurity of a node) and
    \code{Decrease_Gini} (Decrease in Gini impurity when this node is split)
    }
  \item{node_feature}{if return_node_feature = TRUE, returns a TRUE/FALSE matrix
    with ncol(X) columns, each row corresponding to a node in the
    forest. The entries indicate which features were used to define the
    decision rule associated with a node}
  \item{node_data}{if return_node_data is TRUE, returns a TRUE/FALSE
    matrix with nrow(X) columns. Each row corresponds to a node in the
    forest. The entries indicate which data points in X fell in a particular node}
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

  rforest = readForest(rfobj = rf, X=X[train.id,]
                      , Y = as.numeric(Y[train.id])-1, leaf_node_only=FALSE)

  head(rforest$tree_info)

  # only consider leaf nodes
  rforest = readForest(rfobj = rf, X=X[train.id,]
                      , Y = as.numeric(Y[train.id])-1, leaf_node_only=TRUE)

  # count number of leaf nodes with at least 5 observations
  sum(rforest$tree_info$size_node > 5)

  # pass test data through the forest
  rforest = readForest(rfobj=rf, X[test.id,], Y=as.numeric(Y[test.id])-1)

  # count number of leaf nodes with at least 5 observations
  sum(rforest$tree_info$size_node > 5)

}

\author{Sumanta Basu \email{sumbose@berkeley.edu}}

