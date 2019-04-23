#' @importFrom data.table data.table
#' @importFrom Matrix colSums rowSums t
#' @importFrom utils combn
intImportance <- function(int, nf, prec.nd, select.id, weight) {
  # Calculate importance metrics for an interaction across selected elaf nodes
  # of a fitted random forest.
  #   prev: the prevalence of an interaction across all selected leaf nodes,
  #   weighted by observations in each leaf node.
  #   prec: the proportion of class-1 observations in leaf nodes containing the
  #   interaction.
  
  # Remove all 0-weighted leaf nodes from further analysis
  id.rm <- weight == 0
  select.id <- select.id[!id.rm]
  weight <- weight[!id.rm]
  nf <- nf[!id.rm,]
  prec.nd <- prec.nd[!id.rm]

  intord <- length(int)
  if (intord == 1)
    int.id <- nf[, int] != 0
  else
    int.id <- Matrix::rowSums(nf[, int] != 0) == intord

  if (sum(int.id) == 0)
    return(data.table(prev1=0, prev0=0, prec=0))

  # Compute prevalence and precision for given interaction
  prev <- prevalence(weight, int.id, select.id)
  prec <- mean(prec.nd[int.id & select.id])
  return(data.table(prev1=prev[1], prev0=prev[2], prec=prec))
}

prevalence <- function(weight, idint, idcl) {
  # Computes the prevalence of an interaction among class-1 and class-0 leaf 
  # nodes in the fitted RF.
  sint1 <- sum(weight[idint & idcl])
  sint0 <- sum(weight[idint & !idcl])
  s1 <- sum(weight[idcl])
  s0 <- sum(weight[!idcl])

  prev1 <- sint1 / s1
  prev0 <- sint0 / s0
  return(c(prev1, prev0))
}

nodeAttr <- function(read.forest, y, weight=rep(1, length(y))) {
  # Evaluate class proportion of class-1 observations in each leaf node of  a
  # fitted RF.
  if (is.factor(y)) y <- as.numeric(y) - 1

  ndcnt <- t(read.forest$node.obs)
  ndcntY <- Matrix::colSums(ndcnt * y * weight)
  ndcnt <- Matrix::colSums(ndcnt * weight)
  prec.nd <- ndcntY / ndcnt
  prec.nd[ndcnt == 0] <- 0
  return(list(precision=prec.nd, ndcnt=ndcnt))
}

subsetTest <- function(int, ints, importance) {
  # Compare interaction importance metrics to simple baselines.
  # Prevalence: expectation under independent selection.
  # Precision: precision of lower-order interactions.

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  
  # Determine all s-1 order interactions to evaluate
  ss <- combn(int, length(int) - 1, simplify=FALSE)
  setEq <- function(x, y) all(x %in% y) & all(y %in% x)
  setsEq <- function(x, y) sapply(y, setEq, x=x)
  getPair <- function(z) setsEq(z, ints) | setsEq(setdiff(int, z), ints)
  pairs <- lapply(ss, getPair)

  id <- setsEq(int, ints)
  prev <- importance$prev1[id]
  prev.null <- sapply(pairs, function(z) prod(importance$prev1[z]))
  prev.test <- max(mean(prev - prev.null) / prev, -1)
  prec <- importance$prec[id]
  prec.null <- importance$prec[unlist(lapply(pairs, which))]
  prec.test <- max(mean(prec - prec.null) / prec, -1)
  return(data.table(prev.test=prev.test, prec.test=prec.test))
}
