interactPredict <- function(x, int, read.forest, varnames.grp=1:ncol(x), 
                            nrule=1000, min.node=1, mask='low', wt=TRUE, 
                            is.split=FALSE, aggregate=TRUE) {

  # Generate RF predictions for given interactions, using only information
  # from interacting features.
  p <- ncol(x)
  stopifnot(p == length(varnames.grp))
  stopifnot(p == ncol(read.forest$node.feature) / 2)

  nf <- read.forest$node.feature[read.forest$tree.info$size >= min.node,]
  tree.info <- read.forest$tree.info[read.forest$tree.info$size >= min.node,]

  # Get signed and raw feature indices for interaction terms and determine which
  # interaction features are positive
  if (!is.split) int <- strsplit(int, '_')[[1]]
  int.adj <- isPositive(int)
  int.unsgn <- intUnsign(int) 
  
  id.sgn <- mapply(function(i, a) {
    intId(int=i, varnames=varnames.grp, adj=a)
  }, int.unsgn, int.adj, SIMPLIFY=TRUE)
  id.raw <- id.sgn %% p  + p * (id.sgn == p | id.sgn == 2 * p) 
  id.pos <- id.sgn > p
  
  # if classification, subset to class 1 leaf nodes
  if (all(tree.info$prediction %in% 1:2)) {
    nf <- nf[tree.info$prediction == 2,]
    tree.info <- tree.info[tree.info$prediction == 2,]
    tree.info$prediction <- tree.info$prediction - 1
  }
  
  # Subset node feature matrix to and data matrix based on interacting features
  nf <- nf[,id.sgn]
  x <- x[,id.raw]

  if (length(id.sgn) == 1) {
    # order = 1 interaction
    int.nds <- nf != 0
    nf <- as.matrix(nf[int.nds])
    x <- as.matrix(x)
  } else {
    # order > 1 interaction
    nint <- Matrix::rowSums(nf != 0)
    
    # determine which nodes contain desired interaction based on indicated
    # masking: 
    #   low -- only order = s nodes
    #   high -- only order < s nodes
    #   none -- both order = s and order < s nodes
    if (mask == 'low') {
      int.nds <- nint == length(id)
    } else if (mask == 'high') {
      int.nds <- nint < length(id) & nint > 0
    } else if (mask == 'none') {
      int.nds <- nint > 0
    }
    nf <- nf[int.nds,]
    nint <- nint[int.nds]
  }
  
  tree.info <- tree.info[int.nds,]
  if (sum(int.nds) == 0) {
    warning('interaction does not appear on RF paths')
    return(rep(0, nrow(x)))
  }

  # Set response values for each region proportional to node size
  y <- tree.info$prediction
  if (wt) size <- tree.info$size.node
  else size <- rep(1, nrow(tree.info))
 
  # evaluate predictions over subsample of active rules
  nrule <- min(nrule, nrow(nf))
  ss <- sample(nrow(nf), nrule, prob=size)
  
  preds <- matrix(0, nrow=nrow(x), ncol=nrule)
  for (i in 1:nrule) {
    s <- ss[i]
    id.active <- nf[s, ] != 0
    tlow <- t(x[,!id.pos & id.active]) <= nf[s, !id.pos & id.active]
    if (is.null(ncol(tlow))) tlow <- matrix(0, nrow=0, ncol=nrow(x))
    
    thigh <-  t(x[,id.pos & id.active]) > nf[s, id.pos & id.active]
    if (is.null(ncol(thigh))) thigh <- matrix(0, nrow=0, ncol=nrow(x))
    
    int.active <- (colSums(tlow) + colSums(thigh)) == sum(id.active)
    preds[,i] <- int.active * size[s] * y[s]
  }

  if (aggregate) {
    out <- rowSums(preds) / sum(size[ss])
  } else if (length(id.sgn) == 1) {
    out <- rowSums(preds) / sum(size[ss])
    out <- cbind(out, out)
  } else {
    ilow <- nint[ss] < length(id.sgn) 
    if (any(ilow)) 
      out.low <- rowSums(as.matrix(preds[,ilow])) / sum(size[ss[ilow]])
    else
      out.low <- rep(0, nrow(x))
    
    if (!all(ilow))
      out.high <- rowSums(as.matrix(preds[,!ilow])) / sum(size[ss[!ilow]])
    else
      out.high <- rep(0, nrow(x))

    out <- cbind(out.low, out.high)
  }

  return(out)
}

isPositive <- function(x) {
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
