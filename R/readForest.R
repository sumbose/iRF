#' Read forest
#'
#' Read out metadata from random forest decision paths
#'
#' @param rand.forest an object of class randomForest.
#' @param x numeric feature matrix.
#' @param return.node.feature if True, will return sparse matrix indicating
#'  features used on each decision path of the rand.forest.
#' @param return.node.obs if True, will return sparse matrix indicating
#'  observations in x that fall in each leaf node of rand.forest.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param oob.importance if TRUE, importance measures are evaluated on OOB
#'  samples.
#' @param first.split if True, splitting threshold will only be evaluated for
#'  the first time a feature is selected.
#' @param n.core number of cores to use. If -1, all available cores are used.
#'
#' @return a list containing the following entries
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
#' @importFrom data.table rbindlist
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG "%dorng%"
#' @importFrom parallel detectCores
#' @importFrom fastmatch fmatch
readForest <- function(rand.forest, x,
                       return.node.feature=TRUE,
                       return.node.obs=TRUE,
                       varnames.grp=NULL,
                       oob.importance=TRUE,
                       first.split=TRUE,
                       n.core=1) {

  # Check for valid input RF
  if (is.null(rand.forest$forest))
    stop('No Forest component in the random forest object')

  # Get variable names from fitted RF
  varnames <- readVariableNames(rand.forest)

  # Check that RF variables match x
  if (!is.null(colnames(x))) {
    if (any(colnames(x) != varnames))
      stop('variable names of x do not match RF variables')
  } else {
    colnames(x) <- varnames
  }

  # Set variable grouping if not specified
  if (is.null(varnames.grp))
    varnames.grp <- colnames(x)

  # Register cores for parallelization
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)

  ntree <- readNtree(rand.forest)

  n <- nrow(x)
  p <- length(unique(varnames.grp))
  out <- list()

  # Pass observations through RF to determine leaf node membership
  nodes <- readNodes(rand.forest, return.node.obs, x)

  if (n.core == 1) {
    rd.forest <- readTrees(1:ntree,
                           rand.forest=rand.forest,
                           x=x, nodes=nodes,
                           varnames.grp=varnames.grp,
                           oob.importance=oob.importance,
                           return.node.feature=return.node.feature,
                           return.node.obs=return.node.obs,
                           first.split=first.split)
  } else {
    # Split trees across cores to read forest in parallel
    a <- floor(ntree / n.core)
    b <- ntree %% n.core
    ntree.core <- c(rep(a + 1, b), rep(a, n.core - b))
    core.id <- rep(1:n.core, times=ntree.core)

    # Read decision paths across each tree in the forest
    suppressWarnings(
    rd.forest <- foreach(id=1:n.core) %dorng% {
      tree.id <- which(core.id == id)
      readTrees(k=tree.id, rand.forest=rand.forest, x=x, nodes=nodes,
                varnames.grp=varnames.grp,
                oob.importance=oob.importance,
                return.node.feature=return.node.feature,
                return.node.obs=return.node.obs,
                first.split=first.split)
    })
    rd.forest <- unlist(rd.forest, recursive=FALSE)
  }

  # Aggregate node level metadata
  offset <- cumsum(sapply(rd.forest, function(tt) nrow(tt$tree.info)))
  offset <- c(0L, offset[-length(offset)])
  out$tree.info <- rbindlist(lapply(rd.forest, function(tt) tt$tree.info))

  # Aggregate sparse node level feature matrix
  if (return.node.feature) {
    nf <- nfSparse(rd.forest, offset, p)
    dims <- c(nrow(out$tree.info), 2 * p)
    out$node.feature <- sparseMatrix(i=nf$i, j=nf$j, x=nf$x, dims=dims)
  }

  # Aggregate sparse node level observation matrix
  if (return.node.obs) {
    nobs <- nobsSparse(rd.forest, offset, out$tree.info)
    dims <- c(n, nrow(out$tree.info))
    out$node.obs <- sparseMatrix(i=nobs$i, j=nobs$j, dims=dims)
  }

  # Adjust predicted value for randomForest classification
  if (all(class(rand.forest) != 'ranger') &&
      rand.forest$type == 'classification') {
    out$tree.info$prediction <- out$tree.info$prediction - 1
  }

  stopImplicitCluster()
  return(out)
}

readTrees <- function(rand.forest, k, x, nodes,
                      varnames.grp=1:ncol(x),
                      oob.importance=TRUE,
                      return.node.feature=TRUE,
                      return.node.obs=FALSE,
                      first.split=TRUE) {

  out <- lapply(k, readTree, rand.forest=rand.forest, x=x,
                nodes=nodes,varnames.grp=varnames.grp,
                oob.importance=oob.importance,
                return.node.feature=return.node.feature,
                return.node.obs=return.node.obs,
                first.split=first.split)
  return(out)
}

readTree <- function(rand.forest, k, x, nodes,
                     varnames.grp=1:ncol(x),
                     oob.importance=TRUE,
                     return.node.feature=TRUE,
                     return.node.obs=FALSE,
                     first.split=TRUE) {

  n <- nrow(x)
  ntree <- readNtree(rand.forest)

  # Read metadata for current tree
  tree.info <- getTree(rand.forest, k)

  # Set additional metadata features for current tree
  tree.info$node.idx <- 1:nrow(tree.info)
  tree.info$parent <- getParent(tree.info) %% nrow(tree.info)
  tree.info$tree <- as.integer(k)
  tree.info$size.node <- 0L
  select.node <- tree.info$status

  # Read active features for each decision path
  if (return.node.feature) {
    node.feature <- readFeatures(tree.info, varnames.grp=varnames.grp,
                                 first.split=first.split)
  }

  tree.info <- tree.info[select.node, ]

  # Read leaf node membership for each observation
  node.obs <- NULL
  if (return.node.obs) {
    id.leaf <- nodes[, k]

    # Get OOB observations if evaluating importance meausres on OOB
    if (!is.null(rand.forest$inbag) && oob.importance) {
      oob.id <- readOOB(rand.forest, k)
    } else {
      oob.id <- rep(TRUE, nrow(nodes))
      if (oob.importance) warning('keep.inbag = FALSE, using all observations')
    }

    # Group oservations by leaf node membership
    id.obs <- (1:length(id.leaf))[oob.id]
    id.leaf <- id.leaf[oob.id]
    unq.leaf <- sort(unique(id.leaf))
    if (any(class(rand.forest) == 'ranger')) unq.leaf <- unq.leaf + 1
    id <- fmatch(unq.leaf, tree.info$node.idx)
    node.obs <- c(by(id.obs, id.leaf, list))
    names(node.obs) <- id

    tree.info$size.node[id] <- sapply(node.obs, length)
  }

  out <- list()
  out$tree.info <- tree.info[, -(1:5)]
  out$node.feature <- node.feature
  out$node.obs <- node.obs

  return(out)
}


readVariableNames <- function(x, ...) UseMethod("readVariableNames")
readVariableNames.ranger <- function(rand.forest)
    names(rand.forest$variable.importance) 
readVariableNames.randomForest <- function(rand.forest)
    rownames(rand.forest$importance)

readNtree <- function(x, ...) UseMethod("readNtree")
readNtree.ranger <- function(rand.forest) rand.forest$num.trees
readNtree.randomForest <- function(rand.forest) rand.forest$ntree

readNodes <- function(x, ...) UseMethod("readNodes")
readNodes.ranger <- function(rand.forest, return.node.obs, x) {
  pred <- predict(rand.forest, data=x, type='terminalNodes')
  nodes <- pred$predictions
  return(nodes)
}
readNodes.randomForest <- function(rand.forest, return.node.obs, x) {
  if (return.node.obs) {
    pred <- predict(rand.forest, newdata=x, nodes=TRUE)
    nodes <- attr(pred, 'nodes')
    return(nodes)
  } else {
    return(NULL)
  }
}   

readOOB <- function(x, ...) UseMethod("readOOB")
readOOB.ranger <- function(rand.forest, k) rand.forest$inbag[[k]] == 0
readOOB.randomForest <- function(rand.forest, k) rand.forest$inbag[,k] == 0


getParent <- function(tree.info) {
  # Generate a vector of parent node indices from output of getTree
  children <- c(tree.info$`left daughter`, tree.info$`right daughter`)
  parent <- fmatch(1:nrow(tree.info), children)
  parent[1] <- 0
  return(parent)
}

readFeatures <- function(tree.info, varnames.grp,
                         first.split=TRUE) {

  # Pre-allocate variables for path ancestry
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)
  nlf <- sum(tree.info$status)
  nodeVarIndices <- fmatch(varnames.grp, varnames.unq)

  paths <- ancestorPath(tree.info, nodeVarIndices,
                           p, nlf, first.split)

  # Generate sparse matrix of decision path feature selection
  paths <- Matrix(paths, nrow=nlf, byrow=TRUE, sparse=TRUE)
  rownames(paths) <- paths[,ncol(paths)]
  paths <- paths[, 1:(2 * p)]

  # Reorder leaf nodes according to tree.info
  idlf <- tree.info$node.idx[tree.info$status]
  paths <- paths[fmatch(idlf, rownames(paths)),]
  return(paths)
}

nfSparse <- function(rd.forest, offset, p) {
  # Row indices by tree
  nfRow <- function(rf, offset) rf$node.feature@i + 1L + offset
  nf.i <- mapply(nfRow, rd.forest, offset, SIMPLIFY=FALSE)

  # Col indices by tree
  nfCol <- function(rf, p) rep(1:(2 * p), times=diff(rf$node.feature@p))
  nf.j <- lapply(rd.forest, nfCol, p=p)

  # X values y tree
  nfX <- function(rf) rf$node.feature@x
  nf.x <- lapply(rd.forest, nfX)

  return(list(i=.Internal(unlist(nf.i, FALSE, FALSE)),
              j=.Internal(unlist(nf.j, FALSE, FALSE)),
              x=.Internal(unlist(nf.x, FALSE, FALSE))))
}

nobsSparse <- function(rd.forest, offset, tree.info) {
  # Col indices by tree
  nobs.j <- mapply(function(rf, oo) {
    id <- as.numeric(names(rf$node.obs)) + oo
    nrep <- tree.info$size.node[id]
    return(rep(id, times=nrep))
  }, rd.forest, offset)

  # Row indices by tree
  nobs.i <- lapply(rd.forest, function(rf) rf$node.obs)
  return(list(i=unlist(nobs.i), j=unlist(nobs.j)))
}
