readForest <- function(rfobj, x, y=NULL, 
                       return.node.feature=TRUE,
                       wt.pred.accuracy=FALSE,
                       n.core=1){
 
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (wt.pred.accuracy & is.null(y))
    stop('y required to evaluate prediction accuracy')
 
  ntree <- rfobj$ntree
  p <- ncol(x)
  n <- nrow(x)
  out <- list()
  
  # read leaf node data from each tree in the forest 
  rd.forest <- mclapply(1:ntree, readTree, rfobj=rfobj, x=x, y=y,
                        return.node.feature=return.node.feature,
                        wt.pred.accuracy=wt.pred.accuracy,
                        mc.cores=n.core)
 
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  # aggregate sparse feature matrix across forest
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  nf <- aggregateNodeFeature(nf)
  out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], dims=c(max(nf[,1]), p))
  return(out)
  
}

readTree <- function(rfobj, k, x, y, return.node.feature, wt.pred.accuracy) {
  n <- nrow(x)
  p <- ncol(x)
  ntree <- rfobj$ntree
  
  # Read tree level data from RF
  out <- list()
  out$tree.info <- as.data.frame(getTree(rfobj, k))
  out$tree.info$node.idx <- 1:nrow(out$tree.info)
  parents <- getParent(out$tree.info)
  out$tree.info$parent <- as.integer(parents)
  out$tree.info$tree <- as.integer(k)
  out$tree.info$size.node <- 0L 
  
  # replicate each leaf node in node.feature based on specified sampling.
  select.node <- out$tree.info$status == -1
  rep.node <- rep(0, nrow(out$tree.info))
  out$tree.info <- select(out$tree.info, prediction, node.idx, parent, tree, size.node)

  if (is.null(rfobj$obs.nodes)) {
    # if nodes not tracked, pass data through forest to get leaf counts
    fit.data <- passData(rfobj, x, out$tree.info, k)
    leaf.counts <- rowSums(fit.data[select.node,])
    which.leaf <- apply(fit.data[select.node,], MAR=2, which)
    leaf.idx <- as.integer(which(select.node))
    if (wt.pred.accuracy) leaf.sd <- c(by(y, which.leaf, sdNode)) 
  } else {
    leaf.counts <- table(rfobj$obs.nodes[,k])
    leaf.idx <- as.integer(names(leaf.counts)) 
    if (wt.pred.accuracy) leaf.sd <- c(by(y, rfobj$obs.nodes[,k], sdNode))
  }

  out$tree.info$size.node[leaf.idx] <- as.integer(leaf.counts)
  if (wt.pred.accuracy) {
    out$tree.info$dec.purity <- 0
    out$tree.info$dec.purity[leaf.idx] <- pmax((sd(y) - leaf.sd) / sd(y), 0)
  }
  
  out$tree.info <- out$tree.info[select.node,]
  rep.node[select.node] <- 1
  
  # Extract decision paths from leaf nodes as binary sparse matrix
  if (return.node.feature) {
    row.offset <- 0
    var.nodes <- as.integer(rfobj$forest$bestvar[,k])
    total.rows <- n #sum(rep.node[select.node])
    sparse.idcs <- nodeVars(var.nodes, 
                            as.integer(length(select.node)), 
                            as.integer(p),
                            as.integer(parents),
                            as.integer(select.node),
                            as.integer(rep.node),
                            as.integer(row.offset),
                            matrix(0L, nrow=(total.rows * p), ncol=2))
    out$node.feature <- sparse.idcs[!sparse.idcs[,1] == 0,]
  }
  
  return(out)
}


getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent <- parent %% nrow(tree.info)
  parent[1] <- 0
  return(parent)
}


passData <- function(rfobj, x, tt, k) {
  # Pass data through rf object
  leaf.id <- tt$status == -1
  n <- nrow(x)
  n.node <- rfobj$forest$ndbigtree[k]
  node.composition <- matrix(FALSE, nrow=n.node, ncol=n)
  node.composition[1,] <- TRUE
  
  
  for (i in which(!leaf.id)){   
    # determine children and split point for current node
    d.left <- tt$"left daughter"[i]
    d.right <- tt$"right daughter"[i]
    split.var <- tt$"split var"[i]
    split.pt <- tt$"split point"[i]
    
    parent.id <- node.composition[i,]
    d.left.id <- (x[,split.var] <= split.pt) & parent.id
    d.right.id <- (x[,split.var] > split.pt) & parent.id
    
    node.composition[d.left,] <- d.left.id
    node.composition[d.right,] <- d.right.id
  }
  
  return(node.composition)
}

aggregateNodeFeature <- function(nf) {
  # aggregate list of node feature data returned from each tree

  ntree <- length(nf)
  row.offset <- c(0, cumsum(sapply(nf, function(z) max(z[,1])))[-ntree])
  n.rows <- sapply(nf, nrow)
  nf <- do.call(rbind, nf)
  nf[,1] <- nf[,1] + rep(row.offset, times=n.rows)
  return(nf)
}

sdNode <- function(x) {
  sd.node <- ifelse(length(x) == 1, 0, sd(x))
  return(sd.node)
}
