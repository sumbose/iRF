#' @importFrom data.table data.table
#' @importFrom Matrix colSums rowSums t
#' @importFrom utils combn
#' @importFrom fastmatch fmatch
#' @importFrom memoise memoise
intImportance <- function(int, nf, precision, weight) {
  # Calculate importance metrics for an interaction across selected leaf nodes
  # of a fitted random forest.
  #   prev: the prevalence of an interaction across all selected leaf nodes,
  #   weighted by observations in each leaf node.
  #   prec: the proportion of class-1 observations in leaf nodes containing the
  #   interaction.
  
  # Determine which leaf nodes contain the given interaction
  int.id <- Reduce(fast.intersect, nf[as.character(int)])
  if (length(int.id) == 0)
    return(data.frame(prev1=0, prev0=0, prec=0))

  # Compute prevalence and precision for given interaction
  weight1 <- weight * precision
  sint1 <- sum(weight1[int.id])
  s1 <- sum(weight1)
  prev1 <- sint1 / s1

  weight0 <- weight - weight1
  sint0 <- sum(weight0[int.id])
  s0 <- sum(weight0)
  prev0 <- sint0 / s0

  prec <- sint1 / (sint0 + sint1)

  return(data.frame(prev1=prev1, prev0=prev0, prec=prec))
}

# Faster than base::intersect, assuming x and y don't contain duplications.
fast.intersect <- memoise(function(x, y) {
  if (is.null(x) || is.null(y)) return(numeric(0))
  y[fmatch(x, y, 0L)]
})

prevalence <- function(weight, idint) {
  # Computes the prevalence in selected nodes 
  sint <- sum(weight[idint])
  s <- sum(weight)
  return(sint / s)
}

nodePrecision <- function(read.forest, y, count, weights=rep(1, length(y))) {
  # Evaluate class proportion of class-1 observations in each leaf node.
  if (is.factor(y)) y <- as.numeric(y) - 1
  count.y <- Matrix::colSums(read.forest$node.obs * y * weights)
  precision <- count.y / count
  precision[count == 0] <- 0
  return(precision)
}

subsetTest <- function(int, ints, importance) {
  # Compare interaction importance metrics to simple baselines.
  # Prevalence: expectation under independent selection.
  # Precision: precision of lower-order interactions.

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  id.int <- fmatch(list(int), ints)

  # Evaluate prevalence relative to independent selection
  id <- sapply(int, fmatch, ints)
  prev.null <- prod(importance$prev1[id])
  prev <- importance$prev1[id.int]
  prev.test <- prev - prev.null

  # Determine all s-1 order interactions to evaluate
  ss <- combn(int, length(int) - 1, simplify=FALSE)
  id <- sapply(ss, function(ii) fmatch(list(ii), ints))
  prec.null <- importance$prec[id]
  prec <- importance$prec[id.int]
  prec.test <- min(prec - prec.null)

  return(data.table(prev.test=prev.test, prec.test=prec.test))
}
