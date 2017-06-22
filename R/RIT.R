RIT <- function(z, z0, weights=rep(1, nrow(z)), branch=5, depth=10L, n_trees=100L, theta0=0.5, theta1=theta0,
                min_inter_sz=2L, L=100L, n_cores=1L, output_list=FALSE) {
  ## check L, branch, depth, t, min_inter_sz, n_cores, output_list
  L <- as.integer(L)
  branch <- as.double(branch)
  depth <- as.integer(depth)
  n_trees <- as.integer(n_trees)
  theta0 <- as.double(theta0)
  theta1 <- as.double(theta1)
  min_inter_sz <- as.integer(min_inter_sz)
  n_cores <- as.integer(n_cores)
  output_list <- as.logical(output_list)
  weights <- weights / sum(weights)
  
  if (L < 1L)
    stop("L must be >= 1")
  if (branch < 1)
    stop("branch must be >= 1")
  if (depth <= 1L)
    stop("depth must be >= 2")
  if (n_trees <= 0L)
    stop("n_trees must be >= 0")
  if (min_inter_sz < 2L)
    stop("min_inter_sz must be >= 2")
  if (n_cores <1L)
    stop("n_cores must be >= 1")
  if(theta0<0 || theta0>1)
    stop("theta0 must be between 0 and 1")
  if(theta1<0 || theta1>1)
    stop("theta1 must be between 0 and 1")
  if (any(weights < 0))
    stop("weights must be >= 0")
  if (length(weights) != nrow(z))
    stop("length weights must equal nrow z")

  
  ## check z,z0 and convert to suitable data types for 'ExportedcppFunctions.cpp', then carry out RIT
  is_2_class <- !missing(z0)
  is_sparse <- is(z,"Matrix")
  
  # If z or z0 is sparse then also make the other one sparse, so that they are of the same type
  if (is_2_class) {
    if (is_sparse && !is(z0,"Matrix")) {
      z0 <- Matrix(z0, sparse=TRUE)
    }
    else if (is(z0,"Matrix") && !is_sparse) {
      z <- Matrix(z, sparse=TRUE)
      is_sparse <- TRUE
    }
  }
  
  if (is_sparse) {
    if (is_2_class) {
      if (nrow(z0) == 0) stop("z0 must have at least one row")
      if (ncol(z0) != ncol(z)) stop("z and z0 must have the same number of columns")
      z0 <- Matrix::t(z0)
      z0 <- list(z0@i, z0@p)
    }
    if (nrow(z) == 0) stop("z must have at least one row")
    z <- t(z)
    z <- list(z@i, z@p)
  }
  
  ## if z, z0 are not sparse matrices, check that they are matrices
  if (!is_sparse) {
    if (!is.matrix(z)) stop("z must be a matrix")
    if (nrow(z) == 0) stop("z must have more than 0 rows")
    if (is_2_class) {
      if (!is.matrix(z0)) stop("z0 must be a matrix")
      if (nrow(z0) == 0) stop("z0 must have more than 0 rows")
      if (ncol(z) != ncol(z0)) stop("z and z0 must have the same number of columns")
    }
  }

  # carry out RIT_2class
  if (is_2_class) {
    output <- RIT_2class(z, z0, L, branch, depth, n_trees, theta0, theta1, min_inter_sz, n_cores, is_sparse)
    
    # reorder output in decreasing prevalence
    prev_order <- order(output$Class1$Prevalence, decreasing=TRUE)
    prev_order0 <- order(output$Class0$Prevalence, decreasing=TRUE)
    output$Class1$Prevalence <- output$Class1$Prevalence[prev_order]
    output$Class1$Interactions <- output$Class1$Interactions[prev_order]
    output$Class0$Prevalence <- output$Class0$Prevalence[prev_order0]
    output$Class0$Interactions <- output$Class0$Interactions[prev_order0]
    
    # check whether output should be a list or a data.frame
    if(!output_list) {
      data1 <- convert.to.data.frame(output$Class1)
      data0 <- convert.to.data.frame(output$Class0)
      output <- list("Class1"=data1, "Class0"=data0)
    }
    return(output)
  }
  
  # carry out RIT_1class
  output<-RIT_1class(z, weights, L, branch, depth, n_trees, min_inter_sz, n_cores, is_sparse)
  
  # reorder output in decreasing prevalence
  prev_order <- order(output$Prevalence, decreasing=TRUE)
  output$Prevalence <- output$Prevalence[prev_order]
  output$Interactions <- output$Interactions[prev_order]
  
  # check whether output should be a list or a data.frame
  if (!output_list) {
    output<-convert.to.data.frame(output)
  }
  return(output)
}

# Converts the list output of RIT to a data frame
# Not exported
convert.to.data.frame <- function(Inter_Prev_List) {
  str_inter <- character(length(Inter_Prev_List$Interactions))
  for (i in seq_along(Inter_Prev_List$Interactions)) {
    str_inter[i] <- paste(Inter_Prev_List$Interactions[[i]],collapse=" ")
  }
  data <- data.frame(str_inter, Inter_Prev_List$Prevalence, stringsAsFactors = FALSE)
  colnames(data) <- c("Interaction", "Prevalence")
  return(data)
}
