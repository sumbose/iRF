#' Predict interaction
#'
#' Generate predictions from random forest decision rules corresponding to a
#' signed interaction.
#'
#' @param x numeric feature matrix
#' @param int a signed interaction. Formatted as 'X1+_X2+_X3-_...'
#' @param read.forest output of readForest.
#' @param varnames.grp grouping "hyper-features" for RIT search. Features with
#'  the same name will be treated as identical for interaction search.
#' @param min.nd minimum leaf node size to extract decision rules from.
#' @param mask if True, only leaf nodes containing the full interaction will be
#'  used to generate predictions. If false, any leaf node containing a subset of
#'  interacting features will be used to generate predictions.
#'
#' @return a numeric vector of length nrow(x), entries indicating a predicted
#'  response for the corresponding observation. Predictions are generated from
#'  random forest decision rules using only the features in int.
#'  
#' @export
interactPredict <- function(x, int, read.forest, varnames.grp=NULL,
                            min.nd=1, mask=TRUE) {
  
  p <- ncol(x) 
  stopifnot(p == ncol(read.forest$node.feature) / 2)
  varnames.grp <- groupVars(varnames.grp, x)

  # filter out small leaf nodes
  id.keep <- read.forest$tree.info$size.node >= min.nd
  read.forest <- subsetReadForest(read.forest, id.keep)

  # Get indices for interaction features in <x> and <nf>
  int <- strsplit(int, '_')[[1]]
  stopifnot(length(int) > 1)
  int.adj <- isPositive(int)
  int.unsgn <- intUnsign(int) 
  
  id.nf <- mapply(function(i, a) {
    intId(int=i, varnames=varnames.grp, adj=a)
  }, int.unsgn, int.adj, SIMPLIFY=TRUE)
  id.x <- id.nf %% p  + p * (id.nf == p | id.nf == 2 * p) 
  id.pos <- id.nf > p

  # Subset node feature matrix and data matrix based on interacting features 
  tree.info <- read.forest$tree.info
  nf <- read.forest$node.feature[,id.nf]
  x <- x[,id.x]

  # if classification, subset to class 1 leaf nodes
  if (all(tree.info$prediction %in% 1:2)) {
    nf <- nf[tree.info$prediction == 2,]
    tree.info <- tree.info[tree.info$prediction == 2,]
    tree.info$prediction <- tree.info$prediction - 1
  }
  
  # Determine leaf nodes corresponding to the specified interaction 
  nint <- Matrix::rowSums(nf != 0)
  if (mask) {
    # Only evaluate nodes containing the specified interaction
    int.nds <- nint == length(id.x)
  } else {
    # Evaluate nodes containing any interacting features
    int.nds <- nint > 0
  }

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

  # Iterate over sampled leaf nodes and generate predictions using rule
  # associated with specified interaction.
  preds <- numeric(nrow(x))
  for (s in ss) {
    id.active <- nf[s, ] != 0
    tlow <- t(x[,!id.pos & id.active]) <= nf[s, !id.pos & id.active]
    thigh <-  t(x[,id.pos & id.active]) > nf[s, id.pos & id.active]
    
    int.active <- (colSums(tlow) + colSums(thigh)) == sum(id.active)
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


intId <- function(int, varnames.grp, adj) { 
  # Evaluate 1:2p index of interaction term
  int.adj <- which(varnames.grp == int) + adj * length(varnames.grp)
  return(int.adj)
}
