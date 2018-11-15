readForest <- function(rfobj, x, 
                       return.node.feature=TRUE, 
                       return.node.obs=FALSE,
                       varnames.grp=NULL,
                       get.split=FALSE,
                       first.split=TRUE,
                       n.core=-1){
  
  if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')
  if (is.null(varnames.grp)) varnames.grp <- 1:ncol(x)
 
  if (n.core == -1) n.core <- detectCores()
  registerDoParallel(n.core)
  ntree <- rfobj$ntree
  n <- nrow(x)
  p <- length(unique(varnames.grp))
  out <- list()
  
  # Determine leaf nodes for observations in x
  prf <- predict(rfobj, newdata=x, nodes=TRUE)
  nodes <- attr(prf, 'nodes')
  
  # read leaf node data from each tree in the forest 
  suppressWarnings(
  rd.forest <- foreach(tt=1:ntree) %dorng% {
    readTree(tt, rfobj=rfobj, x=x, nodes=nodes,varnames.grp=varnames.grp,
             return.node.feature=return.node.feature, 
             return.node.obs=return.node.obs, get.split=get.split,
             first.split=first.split)
  })
  
  # aggregate node level metadata across forest
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # aggregate sparse node level feature matrix across forest
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  nf <- aggregateNodeFeature(nf)
  
  if (get.split) 
    out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], x=nf[,3], 
                                     dims=c(max(nf[,1]), 2 * p))
  else
    out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2],
                                     dims=c(max(nf[,1]), 2 * p))

 
  # aggregate sparse node level observation matrix across forest
  if (return.node.obs) {
    nobs <- lapply(rd.forest, function(tt) tt$node.obs)
    nobs <- aggregateNodeFeature(nobs)
    out$node.obs <- sparseMatrix(i=nobs[,1], j=nobs[,2], dims=c(max(nf[,1]), n))
  }

  return(out)
  
}

readTree <- function(rfobj, k, x, nodes,
                     varnames.grp=1:ncol(x), 
                     return.node.feature=TRUE,
                     return.node.obs=FALSE,
                     get.split=FALSE,
                     first.split=TRUE) {
  
  n <- nrow(x) 
  ntree <- rfobj$ntree

  # Read tree metadata from forest
  tree.info <- as.data.frame(getTree(rfobj, k))
  n.node <- nrow(tree.info)
  tree.info$node.idx <- 1:nrow(tree.info)
  tree.info$parent <- getParent(tree.info) %% n.node
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  
  # replicate each leaf node in node.feature based on specified sampling.
  select.node <- tree.info$status == -1
  rep.node <- rep(0, nrow(tree.info))
  which.leaf <- nodes[,k]

  if (return.node.feature) {
    node.feature <- ancestorPath(tree.info, varnames.grp=varnames.grp, 
                                 split.pt=get.split, first.split=first.split)
    n.path <- sapply(node.feature, nrow)
    leaf.id <- rep(1:length(node.feature), times=n.path)
    node.feature <- cbind(leaf.id, do.call(rbind, node.feature))
  }

  # if specified, set node counts based on observation weights
  leaf.counts <- table(which.leaf)
  leaf.idx <- as.integer(names(leaf.counts))
  tree.info$size.node[leaf.idx] <- leaf.counts
  
  node.obs <- NULL
  if (return.node.obs) {
    id <- match(which.leaf, sort(unique(which.leaf)))
    node.obs <- cbind(id, 1:n)
    node.obs <- node.obs[order(node.obs[,1]),]
  }
  
  tree.info <- tree.info[select.node,]
  
  out <- list()
  col.remove <- c('left daughter', 'right daughter', 'split var',
                  'split point', 'status')
  out$tree.info <- tree.info[,!colnames(tree.info) %in% col.remove]
  out$node.feature <- node.feature
  out$node.obs <- node.obs
  return(out)
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

varNode <- function(x) {
  var.node <- ifelse(length(x) == 1, 0, var(x))
  return(var.node)
}

getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  parent <- match(1:nrow(tree.info), c(tree.info[,'left daughter'],
                                       tree.info[,'right daughter']))
  parent[1] <- 0
  return(parent)
}


ancestorPath <- function(tree.info, varnames.grp, split.pt=FALSE, 
                         first.split=TRUE) {
 
  # recursively extract path info for all nodes 
  paths <- getAncestorPath(tree.info, varnames.grp, split.pt=split.pt,
                           first.split=first.split)
  
  nlf <- sum(tree.info$status == -1)
  paths <- matrix(unlist(paths), ncol=nlf)
  nn <- paths[nrow(paths),]
  paths <- paths[1:(nrow(paths) - 1),]
  colnames(paths) <- nn

  # subset to only leaf nodes
  paths <- lapply(as.character(which(tree.info$status == -1)), 
                    function(z) {
                      idcs <- which(paths[,z] != 0)
                      return(cbind(id=idcs, att=paths[idcs, z]))
                    })
  return(paths)
}

getAncestorPath <- function(tree.info, varnames.grp, 
                            varnames.unq=unique(varnames.grp), 
                            node.idx=1, p=length(varnames.unq),
                            cur.path=NULL, depth=1L, split.pt=FALSE,
                            first.split=TRUE) {
 
  if (is.null(cur.path)) cur.path <- rep(0L, 2 * p + 1)
  id <- tree.info$`split var`[node.idx]
  sp <- tree.info$`split point`[node.idx]
  node.var <- which(varnames.grp[id] == varnames.unq)
  
  # Generate vector indicating depth/threshold for the first time a variable is 
  # selected on a decision paths
  left.set <- cur.path
  att <- ifelse(split.pt, sp, depth)
  if (first.split) {
    if (all(cur.path[node.var + p] == 0)) left.set[node.var] <- att
  } else {
    left.set[node.var] <- att
  }
  left.child <- tree.info$`left daughter`[node.idx]
  
  right.set <- cur.path
  if (first.split) {
    if (all(cur.path[node.var] == 0)) right.set[node.var + p] <- att
  } else {
    right.set[node.var + p] <- att
  }
  right.child <- tree.info$`right daughter`[node.idx]
  
  if (tree.info$status[node.idx] == -1) {
    cur.path[2*p + 1] <- node.idx
    return(cur.path)
  } else {
    l1 <- getAncestorPath(tree.info, varnames.grp, varnames.unq, 
                          left.child, p, left.set, depth+1, split.pt,
                          first.split)
    
    l2 <- getAncestorPath(tree.info, varnames.grp, varnames.unq,
                          right.child, p, right.set, depth+1, split.pt,
                          first.split)
    return(list(l1, l2))
  }
}

