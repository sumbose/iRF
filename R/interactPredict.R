#' Predict interaction
#'
#' Generate predictions from random forest decision rules corresponding to a
#' signed interaction.
#'
#' @param x numeric feature matrix
#' @param int a signed interaction. Formatted as 'X1+_X2+_X3-_...'
#' @param read.forest output of readForest.
#' @param varnames grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param min.nd minimum leaf node size to extract decision rules from.
#'
#' @return a numeric vector of length nrow(x), entries indicating a predicted
#'  response for the corresponding observation. Predictions are generated from
#'  random forest decision rules using only the features in int.
#'  
#' @export
interactPredict <- function(x, int, read.forest, varnames=NULL, min.nd=1) {
  
  p <- ncol(x) 
  stopifnot(p == ncol(read.forest$node.feature) / 2)
    
  # Set feature names and check for replicates
  if (is.null(colnames(x))) colnames(x) <- paste0('X', 1:ncol(x))
  if (is.null(varnames)) varnames <- colnames(x)

  if (any(duplicated(varnames)))
    stop('Replicate features not supported')
  
  class.rf <- all(read.forest$tree.info$prediction %in% 0:1)

  # Filter out small leaf nodes, and class-0 nodes if classification
  id.keep <- read.forest$tree.info$size.node >= min.nd
  if (class.rf) id.keep <- id.keep & read.forest$tree.info$prediction == 1
  read.forest <- subsetReadForest(read.forest, id.keep)
  
  # Get indices for interaction features in <x> and <nf>
  int <- strsplit(int, '_')[[1]]
  stopifnot(length(int) > 1)
  int.clean <- str_remove_all(int, '[-\\+]')
  int.nf <- int2Id(int, varnames, split=TRUE, signed=TRUE)
  int.x <- int.nf %% p + p * (int.nf %% p == 0)
  int.pos <- int.nf > p

  # Subset node feature matrix and data matrix based on interacting features 
  tree.info <- read.forest$tree.info
  nf <- read.forest$node.feature[,int.nf]
  x <- x[,int.x]

  # Determine leaf nodes corresponding to the specified interaction 
  int.nds <- Matrix::rowMeans(nf != 0) == 1
  
  if (sum(int.nds) < 2) {
    warning('interaction does not appear on RF paths')
    return(rep(0, nrow(x)))
  }

  # Subset to nodes containing the specified interaction
  nf <- nf[int.nds,]
  tree.info <- tree.info[int.nds,]
  
  # Set response values for each region proportional to node size
  y <- tree.info$prediction
  size <- tree.info$size.node
 
  # Sample one decision rule per tree containing the specified interaction
  trees <- unique(tree.info$tree)
  ss <- sapply(trees, sampleTree, tree=tree.info$tree, size=size)
  nrule <- length(ss)

  # Flip sign for negative components of signed interaction
  if (any(!int.pos)) {
    x[, !int.pos] <- -x[, !int.pos]
    nf[, !int.pos] <- -nf[, !int.pos] 
  }

  # Iterate over sampled leaf nodes and generate predictions.
  tx <- t(x)
  preds <- numeric(nrow(x))
  for (s in ss) {
    int.active <- Matrix::colSums(tx > nf[s,]) == length(int)
    preds <- preds + int.active * y[s]
  }
  out <- preds / nrule
  return(out)
}

sampleTree <- function(k, tree, size) { 
  # Sample a leaf node from specified tree
  tree.id <- which(tree == k)
  if (length(tree.id) == 1) return(tree.id)
  sample.id <- sample(tree.id, 1, prob=size[tree.id])
  return(sample.id)
}

isPositive <- function(x) {
  # Generate T/F vector indicating which signed interactions are +
  out <- rep(FALSE, length(x))
  out[grep('+', x, fixed=TRUE)] <- TRUE
  return(out)
}

intUnsign <- function(x) gsub('(-|\\+)$', '', x)


intId <- function(int, varnames, adj) { 
  # Evaluate 1:2p index of interaction term
  int.adj <- which(varnames == int) + adj * length(varnames)
  return(int.adj)
}
