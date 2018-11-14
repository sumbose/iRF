# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=-1, 
                mtry.select.prob=rep(1, ncol(x)),
                interactions.return=NULL, 
                rit.param=list(depth=5, ntree=500, 
                               nchild=2, class.id=1, 
                               min.nd=1, class.cut=NULL), 
                varnames.grp=colnames(x), 
                n.bootstrap=1,
                select.iter=TRUE,
                signed=TRUE,
                verbose=TRUE,
                ...) {
 
  # Check input format
  if (!class(x) %in% c('data.frame', 'matrix'))
    stop('x must be matrix or data frame')
  if (nrow(x) != length(y))
    stop('x and y must contain the same number of observations')
  if (ncol(x) < 2 & (!is.null(interactions.return) | select.iter))
    stop('cannot find interaction - X has less than two columns!')
  if (any(interactions.return > n.iter))
    stop('selected iteration greater than n.iter')
  if (length(mtry.select.prob) != ncol(x))
    stop('length mtry.select.prop must equal number of features')
  if (!is.null(xtest)) {
    if (ncol(xtest) != ncol(x)) 
      stop('training/test data must have same number of features')
  }
  if (!is.null(xtest) & is.null(ytest))
    stop('test set responses not indicated')

  # Check all RIT and set to defaul if missing
  if (is.null(rit.param$depth)) rit.param$depth <- 5
  if (is.null(rit.param$ntree)) rit.param$ntree <- 500
  if (is.null(rit.param$nchild)) rit.param$nchild <- 2
  if (is.null(rit.param$class.id) & is.factor(y)) rit.param$class.id <- 1
  if (is.null(rit.param$min.nd)) rit.param$min.nd <- 1
  if (is.null(rit.param$class.cut) & is.numeric(y)) 
    rit.param$class.cut <- median(y)
 
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  importance <- ifelse(class.irf, 'MeanDecreaseGini', 'IncNodePurity')
   
  rf.list <- list()
  if (!is.null(interactions.return) | select.iter) {
    stability.score <- list()
    importance.score <- list()
  }

  # Set number of trees to grow in each core 
  if (ncore == -1) n.core <- detectCores()
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))

  registerDoMC(n.core)
  for (iter in 1:n.iter) {
    
    # Grow Random Forest on full data
    print(paste('iteration = ', iter))
    suppressWarnings(
    rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dorng% {
                                 randomForest(x, y, 
                                              xtest, ytest, 
                                              ntree=ntree.id[i], 
                                              mtry.select.prob=mtry.select.prob, 
                                              keep.forest=TRUE,
                                              ...)
                               }
    )
    
    # Update feature selection probabilities
    mtry.select.prob <- rf.list[[iter]]$importance

    # Evaluate test set error if supplied
    if (!is.null(xtest) & class.irf & verbose) {
      auroc <- auc(roc(rf.list[[iter]]$test$votes[,2], ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    } else if (!is.null(xtest) & verbose) {
      pct.var <- 1 - mean((rf.list[[iter]]$test$predicted - ytest) ^ 2) / var(ytest)
      pct.var <- max(pct.var, 0)
      print(paste('% var explained:', pct.var * 100))
    }
  }

  # Select iteration for which to return interactions based on minimizing 
  # prediction error on OOB samples
  if (select.iter) {
    interactions.return <- selectIter(rf.list, y=y)
    if (verbose) print(paste('selected iter:', interactions.return))
  }

  # Combine training and test set for RIT, and weight test set to 0. Leaf node
  # class proportions will be evaluated on test set if supplied
  if (!is.null(xtest)) {
    xx <- rbind(x, xtest); yy <- c(y, ytest)
    if (class.irf) yy <- as.factor(yy - 1)
    weights <- c(rep(1, nrow(x)), rep(0, nrow(xtest)))
  } else {
    xx <- x; yy <- y
    weights <- rep(1, nrow(x))
  }

  
  for (iter in interactions.return) {
    # Evaluate interactions in full data random forest
    ints.full <- gRIT(rand.forest=rf.list[[iter]], x=xx, y=yy,
                                weights=weights,
                                varnames.grp=varnames.grp,
                                rit.param=rit.param,
                                signed=signed,
                                n.core=n.core)

    # Find interactions across bootstrap replicates
    if (verbose) cat('finding interactions ... ')
    
    interact.list <- list()
    imp.list <- list()

    for (i.b in 1:n.bootstrap) { 

      sample.id <- bootstrapSample(y)
      # Use feature weights from current iteraction of full data RF
      if (iter == 1) 
        mtry.select.prob <- rep(1, ncol(x))
      else
        mtry.select.prob <- rf.list[[iter - 1]]$importance

      # Fit random forest on bootstrap sample
      suppressWarnings(
      rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                      .multicombine=TRUE, .packages='iRF') %dorng% {
                        randomForest(x[sample.id,], y[sample.id], 
                                     xtest, ytest, 
                                     ntree=ntree.id[i], 
                                     mtry.select.prob=mtry.select.prob, 
                                     keep.forest=TRUE, 
                                     ...)
                      }
      )

      # Run generalized RIT on rf.b to learn interactions
      ints <- gRIT(rand.forest=rf.b, x=xx, y=yy,
                   weights=weights,
                   varnames.grp=varnames.grp,
                   rit.param=rit.param,
                   signed=signed,
                   ints.full=ints.full$int,
                   n.core=n.core)
      
      interact.list[[i.b]] <- ints$int
      imp.list[[i.b]] <- ints$imp
      rm(rf.b)       
    }
    
    # Calculate stability scores of interactions
    stability.score[[iter]] <- summarizeInteract(interact.list)
    importance.score[[iter]] <- summarizeImp(imp.list)
  } # end for (iter in ... )
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)) {
    out$interaction <- stability.score
    out$importance <- importance.score
  }

  if (length(interactions.return) == 1) {
    out$rf.list <- out$rf.list[[interactions.return]]
    out$interaction <- out$interaction[[interactions.return]]
    out$selected.iter <- interactions.return
    out$importance <- out$importance[[interactions.return]]
  }

  return(out)
}


summarizeInteract <- function(x){
  # Aggregate interactions across bootstrap samples

  n.bootstrap <- length(x)
  x <- unlist(x)
  
  if (length(x) >= 1){
    int.tbl <- sort(c(table(x)), decreasing=TRUE)
    int.tbl <- int.tbl / n.bootstrap
    return(int.tbl)
  } else {
    return(c(interaction=numeric(0)))
  }
}

summarizeImp <- function(imp) {
  # Summarize interaction prevalence across bootstrap samples
  nbs <- length(imp)
  
  imp <- rbindlist(imp)
  if (nrow(imp) > 0) {
    imp <- mutate(imp, diff=(prev1-prev0)) %>%
      group_by(int) %>%
      summarize(sta.diff=mean(diff > 0),
                diff=mean(diff),
                sta.prev=mean(prev.test > 0),
                prev1=mean(prev1), 
                prev0=mean(prev0),
                sta.prec=mean(prec.test > 0), 
                prec=mean(prec)) %>%
      arrange(desc(diff))
    return(data.table(imp))
  } else {
    # If no interactions recovered return empty data table
    imp <- data.table(sta.diff=numeric(0), diff=numeric(0),
                       prev1=numeric(0), prev0=numeric(0), 
                       prec=numeric(0))
    return(imp)
  }
}

sampleClass <- function(y, cl, n) {
  # Sample indices specific to a given class
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}

selectIter <- function(rf.list, y) {
  # Evaluate optimal iteration based on prediction error in OOB samples.
  # For classification, error is given by weighted missclassification rate.
  # For regression, error is given by MSE
  predicted <- lapply(rf.list, function(z) 
                      as.numeric(z$predicted) - is.factor(y))
  
  if (is.factor(y)){
    y <- as.numeric(y) - 1
    errFun <- function(y, yhat) {
      wt1 <- mean(y == 1)
      err1 <- sum(y == 0 & yhat == 1)
      wt0 <- mean(y == 0)
      err0 <- sum(y == 1 & yhat == 0)
      out <- wt1 * err1 + wt0 * err0
     return(out)
    }
  } else {
    errFun <- function(y, yhat) mean((yhat - y) ^ 2, na.rm=TRUE)
  }
  
  error <- sapply(predicted, errFun, y=y)
  min.err <- min(error)
  id.select <- max(which(error == min.err))
  return(id.select)
}

bootstrapSample <- function(y) {
  n <- length(y)
  id <- sample(n, n, replace=TRUE)
  return(id)
}
