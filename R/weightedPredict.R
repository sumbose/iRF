weightedPredictForest <- function(rfobj, X, weighted=TRUE) {

  adj <- ifelse(fit$type == 'classification', 1, 0) # classification adjustment
  ntree <- rfobj$ntree
  trees <- lapply(1:ntree, getTree, rfobj=rfobj)
  pred.out <- lapply(trees, function(tt) weightedPredictTree(tt, X))
  pred.tree <- mapply(function(p, tt) tt[p$leaf, 'prediction'], 
                      pred.out, trees) - adj


  if (weighted) {
    weight <- sapply(pred.out, function(p) p$size)
  } else {
    weight <- matrix(1 / ntree, nrow=nrow(X), ncol=ntree)
  }
  weight <- t(apply(weight, MAR=1, function(z) z / sum(z)))
  pred <- rowSums(pred.tree * weight)
  return(pred)
}


weightedPredictTree <- function(rfTree, X){

  tt=as.data.frame(rfTree, check.names=FALSE, stringsAsFactors=FALSE) 
  n=nrow(X)
  p=ncol(X)
  n_node=nrow(tt)

  node_composition=array(FALSE, c(n_node, n))
  tt$size_node=0

  # add parent information
  leaf_id=(tt$status==-1)
  
  # initialize
  node_composition[1,]=TRUE
  tt$size_node[1]=n


  # run  over all non-leaf nodes 
  for (i in which(!leaf_id)){

    d_left=tt$"left daughter"[i]
    d_right=tt$"right daughter"[i]

    split_var=tt$"split var"[i]
    split_pt=tt$"split point"[i]

    parent_id=node_composition[i,]
    d_left_id=(X[,split_var] < split_pt) & parent_id
    d_right_id=(X[,split_var] >= split_pt) & parent_id

    node_composition[d_left,]=d_left_id
    node_composition[d_right,]=d_right_id

    tt$size_node[d_left]=sum(d_left_id)
    tt$size_node[d_right]=sum(d_right_id)

  }

  select_node = tt$status == -1
  node_composition[!select_node,] = FALSE
  leaf_id = apply(node_composition, MAR=2, which)
  size_node = tt$size_node[leaf_id]
  return(list(leaf=leaf_id, size=size_node))
}

