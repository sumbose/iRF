#' Read forest
#'
#' Read out metadata from random forest decision paths
#'
#' @param rand.forest an object of class randomForest.
#' @param x numeric feature matrix
#' @param return.node.feature if True, will return sparse matrix indicating
#'  features used on each decision path of the rand.forest
#' @param return.node.obs if True, will return sparse matrix indicating
#'  observations in x that fall in each leaf node of rand.forest.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param get.split if True, node feature matrix will indicate splitting
#'  threshold for each selected feature.
#' @param first.split if True, splitting threshold will only be evaluated for
#'  the first time a feature is selected.
#' @param n.core number of cores to use. If -1, all available cores are used.
#'
#' @return a list containing the follosing entries
#' \itemize{
#'    \item{tree.info}{data frame of metadata for each leaf node in rand.forest}
#'    \item{node.feature}{optional sparse matrix indicating feature usage on
#'      each decision path}
#'    \item{node.obs}{optional sparse matrix indicating observations appearing
#'      in each leaf node}
#'  }
#'
#' @export
#'
#' @importFrom Matrix Matrix t sparseMatrix rowSums colSums
#' @importFrom data.table data.table rbindlist
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
readForest <- function(rand.forest, x, 
                       return.node.feature=TRUE, 
                       return.node.obs=FALSE,
                       varnames.grp=NULL,
                       get.split=FALSE,
                       first.split=TRUE,
                       n.core=1){
  
  if (is.null(rand.forest$forest))
    stop('No Forest component in the randomForest object')
  varnames.grp <- groupVars(varnames.grp, x)

  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)
  
  ntree <- rand.forest$ntree
  n <- nrow(x)
  p <- length(unique(varnames.grp))
  out <- list()
  
  # Determine leaf nodes for observations in x
  prf <- predict(rand.forest, newdata=x, nodes=TRUE)
  nodes <- attr(prf, 'nodes')
  
  # Read leaf node data from each tree in the forest 
  a <- floor(ntree / n.core)
  b <- ntree %% n.core
  ntree.core <- c(rep(a + 1, b), rep(a, n.core - b))
  core.id <- rep(1:n.core, times=ntree.core)

  suppressWarnings(
  rd.forest <- foreach(id=1:n.core) %dorng% {
    tree.id <- which(core.id == id)
    readTrees(k=tree.id, rand.forest=rand.forest, x=x, 
              nodes=nodes,varnames.grp=varnames.grp,
              return.node.feature=return.node.feature, 
              return.node.obs=return.node.obs, 
              get.split=get.split, first.split=first.split)
  })
  rd.forest <- unlist(rd.forest, recursive=FALSE)

  # Aggregate node level metadata
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))
  
  # Aggregate sparse node level feature matrix
  nf <- lapply(rd.forest, function(tt) tt$node.feature)
  nf <- aggregateNodeFeature(nf)
  
  if (get.split) 
    out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2], x=nf[,3], 
                                     dims=c(max(nf[,1]), 2 * p))
  else
    out$node.feature <- sparseMatrix(i=nf[,1], j=nf[,2],
                                     dims=c(max(nf[,1]), 2 * p))

 
  # Aggregate sparse node level observation matrix
  if (return.node.obs) {
    nobs <- lapply(rd.forest, function(tt) tt$node.obs)
    nobs <- aggregateNodeFeature(nobs)
    out$node.obs <- sparseMatrix(i=nobs[,1], j=nobs[,2], 
                                 dims=c(max(nf[,1]), n))
  }

  stopImplicitCluster()
  return(out)
  
}

readTrees <- function(rand.forest, k, x, nodes,
                      varnames.grp=1:ncol(x),
                      return.node.feature=TRUE,
                      return.node.obs=FALSE,
                      get.split=FALSE,
                      first.split=TRUE) {

  out <- lapply(k, readTree, rand.forest=rand.forest, x=x, 
                nodes=nodes,varnames.grp=varnames.grp,
                return.node.feature=return.node.feature, 
                return.node.obs=return.node.obs, 
                get.split=get.split, first.split=first.split)
  return(out)
}

readTree <- function(rand.forest, k, x, nodes,
                     varnames.grp=1:ncol(x), 
                     return.node.feature=TRUE,
                     return.node.obs=FALSE,
                     get.split=FALSE,
                     first.split=TRUE) {
  
  n <- nrow(x) 
  ntree <- rand.forest$ntree

  # Read tree metadata from forest
  tree.info <- as.data.frame(getTree(rand.forest, k))
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

