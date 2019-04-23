permImportance <- function(rfobj, x, y, ints, n.perms=3, varnames.group=NULL,
                             n.cores=1) {
  # Evaluate the importance of an interaction by permuting all other features
  # args:
  #   rfobj: fitted random forest
  #   x: raw data matrix
  #   int: a character vector indicating the features for which to 
  #     use raw data, separated by '_'
  #   n.perm: number of permutations
  #   varnames.group: feature names corresponding to columns of x
  #   collapse, should predictions be averaged over permutations?
  #   n.cores: number of cores to parallelize over

  if (is.null(varnames.group) & !is.null(colnames(x)))
    varnames.group <- colnames(x)
  if (is.null(colnames(x)))
    varnames.group <- as.character(1:ncol(x))
  
  preds <- permPredict(rfobj, x, ints, n.perms, varnames.group, 
                       collapse=TRUE, n.cores=n.cores)
  interact.score <- apply(preds, MAR=2, predAccuracy, y=y)
  return(interact.score)
}

permPredict <- function(rfobj, x, int, n.perms=3, varnames.group=NULL, 
                        collapse=TRUE, n.cores=1) {
  # Wrapper function for predictInt, to evaluate predictions over multiple 
  # interactions
  # args:
  #   rfobj: fitted random forest
  #   x: raw data matrix
  #   int: a character vector indicating the features for which to 
  #     use raw data, separated by '_'
  #   n.perm: number of permutations
  #   varnames.group: feature names corresponding to columns of x
  #   collapse, should predictions be averaged over permutations?
  #   n.cores: number of cores to parallelize over
  
  n.cores <- min(n.cores, length(int))
  if (is.null(varnames.group) & !is.null(colnames(x)))
    varnames.group <- colnames(x)
  
  pred.ints <- mclapply(int, predictInt, rfobj=rfobj, x=x, n.perms=n.perms, 
                        varnames.group=varnames.group, mc.cores=n.cores)
  
  if (collapse) {
    pred.ints <- sapply(pred.ints, rowMeans)
    colnames(pred.ints) <- int
  }
  
 return(pred.ints)
}

predictInt <- function(rfobj, x, int, n.perms=3, varnames.group=NULL) {
  # Predict responses with a subset of feautres fixed to raw data values and 
  # the remaining features permuted
  # args:
  #   rfobj: fitted random forest
  #   x: raw data matrix
  #   int: a length 1 character vector indicating the features for which to 
  #     use raw data, separated by '_'
  #   n.perm: number of permutations
  #   varnames.group: feature names corresponding to columns of x
  
  pred.type <- ifelse(rfobj$type == 'regression', 'response', 'prob')
  x.perm <- replicate(n.perms, permuteVars(x, 1:ncol(x)), simplify=FALSE)
  x.fixed <- lapply(x.perm, fixInteract, x=x, int=int, x.names=varnames.group)
  pred.fixed <- lapply(x.fixed, predict, object=rfobj, type=pred.type)
  if (pred.type == 'prob') pred.fixed <- lapply(pred.fixed, function(p) p[,2])
  pred.fixed <- do.call(cbind, pred.fixed)
  return(pred.fixed)
}

predAccuracy <- function(y.hat, y) {
  # Evaluate prediction accuracy: AUROC for classification and decrease 
  # variance for regression.
  # args:
  #   y.hat: predicted response
  #   y: true response
  require(AUC)
  if (is.factor(y)) {
    accuracy <- auc(roc(y.hat, y))
  } else {
    accuracy <- 1 - mean((y - y.hat) ^ 2) / var(y)
    if (accuracy < 0) accuracy <- 0
  }
  return(accuracy)
}

fixInteract <- function(x, x.perm, int, x.names=NULL) {
  # Returns feature matrix where all variables in int have been fixed to their 
  # raw values and all other variables have been permuted.
  # args:
  #   x: raw data matrix
  #   x.perm: data matrix with all columns permuted
  #   int: a length 1 character vector indicating the features for which to 
  #     use raw data, separated by '_'
  #   x.names: character vector indicating feature names
  stopifnot(length(int) == 1)
  
  int.split <- unlist(strsplit(int, '_'))
  if (!is.null(x.names)) {
    stopifnot(length(x.names) == ncol(x))
    int.split <- which(x.names %in% int.split)
  } else if (!is.null(colnames(x))) {
    int.split <- which(colnames(x) %in% int.split)
  } else {
    int.split <- as.numeric(int.split)
  }
  
  x.perm[,int.split] <- x[,int.split]
  return(x.perm)
}

permuteVars <- function(x, vars) {
  # Permutes the columns of x corresponding to vars
  x[,vars] <- apply(x[,vars], MAR=2, sample)
  return(x)
}
