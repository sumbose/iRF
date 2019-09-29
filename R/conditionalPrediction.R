conditionalPred <- function(rfobj, rd.forest, x, y, ints, 
                            varnames.group=NULL, n.cores=1) {
  # Evaluate interaction based on prediction accuracy where predictions are 
  # made using only leaf nodes for which the given interaction falls on the 
  # decision path.
  require(parallel)
  if (is.null(varnames.group) && !is.null(colnames(x)))
    varnames.group <- colnames(x)
  
  y.hat <- mclapply(ints, predIntForest, rfobj=rfobj, rd.forest=rd.forest, 
                    x=x, y=y, varnames.group=varnames.group, mc.cores=n.cores)
  accuracy <- sapply(y.hat, predAccuracy, y=y)
  return(accuracy)
}

predIntForest <- function(rfobj, rd.forest, x, y, int, varnames.group) {
  # Predict responses from a RF using only leaf nodes for which a given 
  # interaction falls on the decision path.
  avg.response <- ifelse(is.factor(y), mean(as.numeric(y) - 1), mean(y))
  rd.forest$tree.info$forest.idx <- 1:nrow(rd.forest$tree.info)

  preds <- predict(rfobj, newdata=x, predict.all=TRUE, nodes=TRUE)
  node.mat <- attr(preds, 'nodes')
  interact.nodes <- getInteractNodes(nf=rd.forest$node.feature, 
                                     x.names=varnames.group, int=int)
  
  tree.preds <- sapply(1:rfobj$ntree, predIntTree, 
                          pred.tree=preds$individual, 
                          node.mat=node.mat, 
                          interact.nodes=interact.nodes, 
                          avg.response=avg.response,
                          rd.forest=rd.forest)
  interact.pred <- rowMeans(tree.preds)
} 

predIntTree <- function(pred.tree, node.mat, interact.nodes, 
                        avg.response, rd.forest, tree.idx) {
  # Predict responses from a decision tree using only leaf 
  # nodes for which a given interaction falls on the decision 
  # path.
  require(dplyr)
  # Get node indices for paths with full interaction
  tree.info <- filter(rd.forest$tree.info, tree == tree.idx)
  tree.interact <- interact.nodes[tree.info$forest.idx]
  tree.interact.nodes <- tree.info$node.idx[tree.interact]
  
  # Get predictions of observations that fall in interaction nodes
  tree.nodes <- node.mat[,tree.idx]
  is.interact <- tree.nodes %in% tree.interact.nodes
  tree.preds <- as.numeric(pred.tree[,tree.idx])
  tree.preds[!is.interact] <- avg.response
  return(tree.preds)
}

getInteractNodes <- function(nf, x.names, int) {
  # Determine which leaf nodes contain a given interactions along their 
  # decision paths
  int.split <- strsplit(int, '_')[[1]]
  if (!is.null(x.names)) {
    # group feature matrix by replicated variables
    grp.names <- unique(x.names)
    makeGroup <- function(x, g) apply(as.matrix(x[,x.names == g]), MAR=1, max)
    nf <- sapply(grp.names, makeGroup, x=nf)
    is.interact <- apply(nf[,int.split], MAR=1, sum) == length(int.split)
  } else {
    int.split <- as.numeric(int.split)
    is.interact <- apply(nf[,int.split], MAR=1, sum) == length(int.split)
  }
  return(is.interact)
}
