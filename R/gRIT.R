#' Generalized random intersection trees
#'
#' Run RIT across decision paths of a fitted random forest.
#'
#' @param x numeric feature matrix
#' @param y response vector. If factor, classification is assumed.
#' @param rand.forest an object of class randomForest. Required if read.forest
#'  is NULL.
#' @param read.forest output of readForest. Required if rand.forest is NULL.
#' @param rit.param named list specifying RIT parameters. Entries include
#'  \code{depth}: depths of RITs, \code{ntree}: number of RITs, \code{nchild}:
#'  number of child nodes for each RIT, \code{class.id}: 0-1 indicating which
#'  leaf nodes RIT should be run over, \code{min.nd}: minimum node size to run
#'  RIT over, \code{class.cut}: threshold for converting leaf nodes in
#'  regression to binary classes.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param weights numeric weight for each observation. Leaf nodes will be
#'  sampled for RIT with probability proprtional to the total weight of
#'  observations they contain.
#' @param signed if TRUE, signed interactions will be returned
#' @param oob.importance if TRUE, importance measures are evaluated on OOB
#'  samples.
#' @param ints.eval interactions to evaluate. If specified, importance metrics
#'  will be evaluated for these interactions instead of those recovered by RIT.
#' @param ints.idx.eval like \code{ints.eval}, but specifies the indice of the
#'  interactions instead of their names. Intended for internal use only.
#' @param n.core number of cores to use. If -1, all available cores are used.
#'
#' @return a data table containing the recovered interactions and importance
#'  scores.
#'
#' @export
#'
#' @importFrom dplyr mutate right_join
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom data.table data.table as.data.table
gRIT <- function(x, y,
                 rand.forest=NULL,
                 read.forest=NULL,
                 rit.param=list(depth=5,
                         ntree=500,
                         nchild=2,
                         class.id=1,
                         min.nd=1,
                         class.cut=NULL),
                 varnames.grp=colnames(x),
                 weights=rep(1, nrow(x)),
                 signed=TRUE,
                 oob.importance=TRUE,
                 ints.idx.eval = NULL,
                 ints.eval=NULL,
                 n.core=1) {

  class.irf <- is.factor(y)
  if (n.core == -1) n.core <- detectCores()
  if (n.core > 1) registerDoParallel(n.core)

  # Check rit parameters and set default values if not specified
  if (is.null(rand.forest) && is.null(read.forest))
    stop('Supply random forest or read forest output')
  if (is.null(rit.param$depth))
    rit.param$depth <- 5
  if (is.null(rit.param$ntree))
    rit.param$ntree <- 500
  if (is.null(rit.param$nchild))
    rit.param$nchild <- 2
  if (is.null(rit.param$class.id) && class.irf)
    rit.param$class.id <- 1
  if (is.null(rit.param$min.nd))
    rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) && !class.irf)
    rit.param$class.cut <- median(y)
  if (!class.irf)
    y <- as.numeric(y >= rit.param$class.cut)

  if (!is.null(ints.eval) && !is.null(ints.idx.eval))
    warning(paste('Both ints.eval and ints.idx.eval are specified,',
                  'ignoring the latter'))

  # Set feature names for grouping interactions
  if (is.null(varnames.grp) && !is.null(colnames(x)))
    varnames.grp <- colnames(x)
  else if (is.null(varnames.grp))
    varnames.grp <- as.character(1:ncol(x))

  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)

  # Read RF object to extract decision path metadata
  if (is.null(read.forest)) {
    read.forest <- readForest(rand.forest, x=x,
                              oob.importance=oob.importance,
                              return.node.feature=TRUE,
                              return.node.obs=TRUE,
                              varnames.grp=varnames.grp,
                              n.core=n.core)
  }

  # Collapse node feature matrix for unsigned iRF
  if (!signed) read.forest$node.feature <- collapseNF(read.forest$node.feature)

  # Evaluate leaf node size and subset forest based on minimum node size
  count <- read.forest$tree.info$size.node
  idcnt <- count >= rit.param$min.nd
  read.forest <- subsetReadForest(read.forest, idcnt)
  count <- count[idcnt]

  # Evaluate leaf node precision
  precision <- nodePrecision(read.forest, y, count, weights)

  # Select class specific leaf nodes
  if (class.irf)
    idcl <- read.forest$tree.info$prediction == rit.param$class.id
  else
    idcl <- read.forest$tree.info$prediction > rit.param$class.cut

  if (sum(idcl) < 2) {
    return(nullReturnGRIT())
  } else {

    # Run RIT on leaf nodes of selected class
    ints <- runRIT(subsetReadForest(read.forest, idcl),
                   weights=count[idcl] * precision[idcl],
                   rit.param=rit.param, n.core=n.core)

    if (length(ints) == 0) return(nullReturnGRIT())

    # Set recovered interactions or convert to indices if supplied
    if (!is.null(ints.eval)) {
      ints.eval <- lapply(ints.eval, int2Id, signed=signed,
                          varnames.grp=varnames.grp)
    } else if (!is.null(ints.idx.eval)) {
      ints.eval <- ints.idx.eval
    } else {
      ints.eval <- ints
    }

    # Evaluate importance metrics for interactions and lower order subsets.
    ints.eval <- lapply(ints.eval, unname)
    ints.sub <- lapply(ints.eval, intSubsets)
    ints.sub <- unique(unlist(ints.sub, recursive=FALSE))

    # Convert node feature matrix to list of active features for fast lookup
    node.feature <- as(read.forest$node.feature, 'dgTMatrix')
    nf.list <- split(node.feature@i + 1L, node.feature@j + 1L)
    names(nf.list) <- unique(node.feature@j) + 1L

    ximp <- lapply(ints.sub, intImportance, nf=nf.list, weight=count,
                   precision=precision)
    ximp <- rbindlist(ximp)

    imp.test <- lapply(ints.eval, subsetTest, importance=ximp, ints=ints.sub)
    imp.test <- rbindlist(imp.test)
  }

  # Aggregate evaluated interations for return
  ints.recovered <- nameInts(ints, varnames.unq, signed=signed)

  name.full <- nameInts(ints.eval, varnames.unq, signed=signed)
  imp.test <- mutate(imp.test, int=name.full)

  name.sub <- nameInts(ints.sub, varnames.unq, signed=signed)
  id.recovered <- name.sub %in% ints.recovered
  ximp <- mutate(ximp, int.idx=ints.sub, int=name.sub,
                 recovered=id.recovered) %>%
          right_join(imp.test, by='int')

  stopImplicitCluster()

  return(as.data.table(ximp))
}

runRIT <- function(read.forest, weights, rit.param, n.core=1) {
  # Run a weighted version of RIT across RF decision paths
  xrit <- read.forest$node.feature[weights > 0,]
  interactions <- RIT(xrit,
                      weights=weights[weights > 0],
                      depth=rit.param$depth,
                      n_trees=rit.param$ntree,
                      branch=rit.param$nchild,
                      output_list=TRUE,
                      n_cores=n.core)$Interaction

  return(interactions)
}

subsetReadForest <- function(read.forest, subset.idcs) {
  # Subset nodes from readforest output
  if (!is.null(read.forest$node.feature))
    read.forest$node.feature <- read.forest$node.feature[subset.idcs,]

  if (!is.null(read.forest$tree.info))
    read.forest$tree.info <- read.forest$tree.info[subset.idcs,]

  if (!is.null(read.forest$node.obs))
    read.forest$node.obs <- read.forest$node.obs[,subset.idcs]

  return(read.forest)
}

collapseNF <- function(x) {
  # Convert from signed to unsigned node feature matrix
  p <- ncol(x) / 2
  x <- x[,1:p] + x[,(p + 1):(2 * p)]
  return(x)
}

nullReturnGRIT <- function() {
  # Return empty interaction and importance
  out <- data.table(prev1=numeric(0), prev0=numeric(0),
                    prec=numeric(0), int=character(0),
                    recovered=logical(0), prev.test=numeric(0),
                    prec.test=numeric(0))
  return(out)
}
