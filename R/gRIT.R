gRIT <- function(x, y, 
                 rand.forest=NULL, 
                 read.forest=NULL,
                 weights=rep(1, nrow(x)),
                 varnames.grp=colnames(x),
                 rit.param=list(depth=5,
                         ntree=500,
                         nchild=2,
                         class.id=1,
                         min.nd=1,
                         class.cut=NULL),
                 signed=TRUE,
                 ints.full=NULL,
                 n.core=-1) {

  out <- list()
  class.irf <- is.factor(y)
  if (ncore == -1) n.core <- detectCores()  
  registerDoMC(n.core)

  # Check rit parameters and set default values if not specified
  if (is.null(rand.forest) & is.null(read.forest))
    stop('Supply random forest or read forest output')
  if (!all(weights %in% 0:1))
    stop('Only 0-1 weighting supported')
  if (is.null(rit.param$depth)) 
    rit.param$depth <- 5
  if (is.null(rit.param$ntree)) 
    rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) 
    rit.param$nchild <- 2
  if (is.null(rit.param$class.id) & class.irf) 
    rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) 
    rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & !class.irf) 
    rit.param$class.cut <- median(y)

  # Set feature names for grouping interactions
  if (is.null(varnames.grp) & !is.null(colnames(x)))
    varnames.grp <- colnames(x)
  else if (is.null(varnames.grp))
    varnames.grp <- as.character(1:ncol(x))
  varnames.unq <- unique(varnames.grp)
  p <- length(varnames.unq)
  
  # Read RF object to extract decision path metadata
  if (is.null(read.forest)) {
    read.forest <- readForest(rand.forest, x=x,
                              return.node.feature=TRUE,
                              return.node.obs=TRUE,
                              varnames.grp=varnames.grp,
                              get.split=TRUE,
                              n.core=n.core)
  }

  # Collapse node feature matrix for unsigned iRF
  if (!signed) 
    read.forest$node.feature <- collapseNF(read.forest$node.feature)

  # Evaluate leaf node attributes: number of observations in node, proportion of
  # class-1 observations in node
  yprec <- precision(read.forest, y, weights)
  ndcnt <- Matrix::colSums(t(read.forest$node.obs) * weights)
  rit.param$min.nd <- min(rit.param$min.nd, quantile(ndcnt, prob=0.9))
  idcnt <- ndcnt >= rit.param$min.nd
  ndcnt <- ndcnt[idcnt]
  read.forest <- subsetReadForest(read.forest, idcnt)
  
  # Select class specific leaf nodes
  if (class.irf)
    idcl <- read.forest$tree.info$prediction == rit.param$class.id + 1
  else
    idcl <- read.forest$tree.info$prediction > rit.param$class.cut

  if (sum(idcl) < 2) {
    return(nullReturn())
  } else {
  
    # Run RIT on leaf nodes of selected class  
    out$int <- runRIT(subsetReadForest(read.forest, idcl),
                      weights=ndcnt[idcl],
                      rit.param=rit.param,
                      n.core=n.core)

    # Evaluate prevalence of recovered/supplied interactions 
    if (!is.null(out$int)) {

      if (is.null(ints.full)) {
        ints.full <- out$int
      } else {
        ints.full <- lapply(ints.full, int2Id, 
                            varnames.grp=varnames.grp, 
                            signed=signed)
      }
     
      ints.sub <- lapply(ints.full, intSubsets)
      ints.sub <- unique(unlist(ints.sub, recursive=FALSE))
      suppressWarnings(
      imp <- foreach(int=ints.sub) %dorng% {
        intImportance(int, nf=read.forest$node.feature,
                      yprec=yprec, select.id=idcl, 
                      weight=ndcnt)
      })
      out$imp <- rbindlist(imp) 
      imp.test <- lapply(ints.full, subsetTest, importance=out$imp, ints=ints.sub)
      imp.test <- rbindlist(imp.test)
    }

    if (is.null(out$int)) return(nullReturn())
    out$int <- nameInts(out$int, varnames.unq, signed=signed)
    imp.test <- imp.test %>%
      mutate(int=nameInts(ints.full, varnames.unq, signed=signed))
    out$imp <- out$imp %>%
      mutate(int=nameInts(ints.sub, varnames.unq, signed=signed)) %>%
      right_join(imp.test, by='int')
    out$int <- out$int[out$int %in% out$imp$int]
  }
  return(out)
}

runRIT <- function(read.forest, weights, rit.param, n.core=1) {
  # Run a weighted version of RIT across RF decision paths
  xrit <- cbind(read.forest$node.feature[weights > 0,])
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
    read.forest$node.obs <- read.forest$node.obs[subset.idcs,]

  return(read.forest)
}

collapseNF <- function(x) {
  # Convert from signed to unsigned node feature matrix
  p <- ncol(x) / 2
  x <- x[,1:p] + x[,(p + 1):(2 * p)]
  return(x)
}

nullReturn <- function() {
  # Return empty interaction and importance
  out <- list()
  out$int <- character(0)
  out$imp <- data.table(prev1=numeric(0), prev0=numeric(0),
                        prec=numeric(0), prev.test=numeric(0),
                        prec.test=numeric(0), int=character(0))
  return(out)
}
