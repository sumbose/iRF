intImportance <- function(int, nf, yprec, select.id, weight) {
  # Calculate the prevalence of an interaction across selected leaf nodes of a
  # random forest
  id.rm <- weight == 0
  select.id <- select.id[!id.rm]
  weight <- weight[!id.rm]
  nf <- nf[!id.rm,]
  yprec <- yprec[!id.rm]

  intord <- length(int)
  if (intord == 1)
    int.id <- nf[, int] != 0
  else
    int.id <- Matrix::rowSums(nf[, int] != 0) == intord

  # Compute conditional probabilities of interaction given class and class given
  # interaction
  if (sum(int.id) == 0)
    return(data.table(prev1=0, prev0=0, prec=0))

  prev <- prevalence(weight, int.id, select.id)
  prec <- mean(yprec[int.id & select.id], na.rm=TRUE)
  return(data.table(prev1=prev[1], prev0=prev[2], prec=prec))
}

prevalence <- function(weight, idint, idcl) {
  sint1 <- sum(weight[idint & idcl])
  sint0 <- sum(weight[idint & !idcl])
  s1 <- sum(weight[idcl])
  s0 <- sum(weight[!idcl])

  prev1 <- sint1 / s1
  prev0 <- sint0 / s0
  return(c(prev1, prev0))
}

precision <- function(read.forest, y, weights) {
  # Evaluate class proportion in each leaf node
  class.irf <- is.factor(y)
  if (class.irf) y <- as.numeric(y) - 1

  stopifnot(all(weights %in% 0:1))
  if (!all(weights == 1)) weights <- 1 - weights

  # TODO: implement for regression
  ndcnt <- t(read.forest$node.obs) * weights
  ndcntY <- Matrix::colSums(ndcnt * y)
  ndcnt <- Matrix::colSums(ndcnt)
  yprec <- ndcntY / ndcnt
  yprec[is.nan(yprec)] <- 0
  return(yprec)
}

subsetTest <- function(int, ints, importance) {
  # Compare prevalence of interaction on decision paths to expectation under
  # independent selection

  if (length(int) == 1) return(c(prev.test=0, prec.test=0))
  ss <- combn(int, length(int) - 1, simplify=FALSE)

  setEq <- function(x, y) all(x %in% y) & all(y %in% x)
  setsEq <- function(x, y) sapply(y, setEq, x=x)
  getPair <- function(z) setsEq(z, ints) | setsEq(setdiff(int, z), ints)
  pairs <- lapply(ss, getPair)

  id <- setsEq(int, ints)
  prev <- importance$prev1[id]
  prev.null <- sapply(pairs, function(z) prod(importance$prev1[z]))
  prec <- importance$prec[id]
  prec.null <- importance$prec[unlist(lapply(pairs, which))]
  return(data.table(prev.test=min(prev - prev.null),
                    prec.test=min(prec - prec.null)))
}
