# Iteratively grows random forests, finds case specific feature interactions
iRF <- function(x, y, 
                xtest=NULL, ytest=NULL, 
                n.iter=5, 
                ntree=500, 
                n.core=1, 
                mtry.select.prob=rep(1/ncol(x), ncol(x)), 
                keep.impvar.quantile=NULL, 
                interactions.return=NULL, 
                wt.pred.accuracy=FALSE, 
                cutoff.unimp.feature=0,  
                rit.param=list(depth=5, ntree=100, nchild=2, class.id=1, class.cut=NULL), 
                varnames.grp=NULL, 
                n.bootstrap=30,
                bootstrap.forest=TRUE, 
                verbose=TRUE,
               ...) {
 

  if (!is.matrix(x) | (!is.null(xtest) & !is.matrix(xtest)))
    stop('either x or xtest is not a matrix !')
  if (!is.numeric(x) | (!is.null(xtest) & !is.numeric(xtest)))
    stop('either x or xtest is not a numeric matrix!')
  if (ncol(x) < 2 & !is.null(interactions.return))
    stop('cannot find interaction - X has less than two columns!')
  if (any(interactions.return > n.iter))
    stop('interaction iteration to return greater than n.iter')
 
  n <- nrow(x)
  p <- ncol(x)
  class.irf <- is.factor(y)
  if (n.core > 1) registerDoParallel(cores=n.core)  

  rf.list <- list()
  if (!is.null(interactions.return)) {
    interact.list <- list()
    stability.score <- list()
    prevalence <- list()
  }
  
  # Set number of trees to grow in each core
  a <- floor(ntree / n.core) 
  b <- ntree %% n.core
  ntree.id <- c(rep(a + 1, b), rep(a, n.core - b))
  
  for (iter in 1:n.iter) {
    
    ## 1: Grow Random Forest on full data
    print(paste('iteration = ', iter))
    rf.list[[iter]] <- foreach(i=1:length(ntree.id), .combine=combine, 
                               .multicombine=TRUE, .packages='iRF') %dopar% {
                                 randomForest(x, y, 
                                              xtest, ytest, 
                                              ntree=ntree.id[i], 
                                              mtry.select.prob=mtry.select.prob, 
                                              keep.forest=TRUE,
                                              track.nodes=TRUE, 
                                              ...)
                               }
    
    ## 2.1: Find interactions across bootstrap replicates
    if (iter %in% interactions.return){
      if (verbose){cat('finding interactions ... ')}
      
      interact.list.b <- list()      
      for (i.b in 1:n.bootstrap) { 

        if (class.irf) {
          n.class <- table(y)
          sample.id <- mapply(function(cc, nn) sampleClass(y, cc, nn),
                              as.factor(names(n.class)), n.class)
          sample.id <- unlist(sample.id)
        } else {
          sample.id <- sample(n, n, replace=TRUE)
        }
        
        if (bootstrap.forest) {
          i <- NULL 
          #2.1.1: fit random forest on bootstrap sample
          rf.b <- foreach(i=1:length(ntree.id), .combine=combine, 
                          .multicombine=TRUE, .packages='iRF') %dopar% {
                            randomForest(x[sample.id,], y[sample.id], 
                                         xtest, ytest, 
                                         ntree=ntree.id[i], 
                                         mtry.select.prob=mtry.select.prob, 
                                         keep.forest=TRUE, 
                                         track.nodes=TRUE, 
                                         ...)
                        }
        } else {
          rf.b <- rf.list[[iter]]
        }
        
        
        #2.1.2: run generalized RIT on rf.b to learn interactions
        ints <- generalizedRIT(rf=rf.b, 
                               x=x[sample.id,], y=y[sample.id],  
                               wt.pred.accuracy=wt.pred.accuracy,
                               class.irf=class.irf, 
                               varnames.grp=varnames.grp,
                               cutoff.unimp.feature=cutoff.unimp.feature,
                               rit.param=rit.param,
                               n.core=n.core) 
        interact.list.b[[i.b]] <- ints
        rm(rf.b)       
        
      }
     
      interact.list[[iter]] <- interact.list.b 
      # 2.2: calculate stability scores of interactions
      if (!is.null(varnames.grp))
        varnames.new <- unique(varnames.grp)
      else if (!is.null(colnames(x)))
        varnames.new <- colnames(x)
      else
        varnames.new <- 1:ncol(x)
      
      summary.interact <- summarizeInteract(interact.list[[iter]], varnames=varnames.new)
      stability.score[[iter]] <- summary.interact$interaction
    
    } # end if (find_interaction)
    
    ## 3: update mtry.select.prob 
    if (!class.irf) 
      mtry.select.prob <- rf.list[[iter]]$importance[,'IncNodePurity']
    else
      mtry.select.prob <- rf.list[[iter]]$importance[,'MeanDecreaseGini']
    
    
    if (!is.null(xtest) & class.irf){
      auroc <- auc(roc(rf.list[[iter]]$test$votes[,2], ytest))
      print(paste('AUROC: ', round(auroc, 2)))
    }
    
  } # end for (iter in ... )
  
  
  out <- list()
  out$rf.list <- rf.list
  if (!is.null(interactions.return)){
    out$interaction <- stability.score
  }
  return(out)
}


generalizedRIT <- function(rf, x, y, wt.pred.accuracy, class.irf, varnames.grp,
                           cutoff.unimp.feature, rit.param, n.core) {

  # Extract decision paths from rf as sparse binary matrix to be passed to RIT
  rforest <- readForest(rf, x=x, y=y, 
                        return.node.feature=TRUE,
                        wt.pred.accuracy=wt.pred.accuracy, 
                        n.core=n.core)
  class.id <- rit.param$class.id 

  # Select class specific leaf nodes
  select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  if (class.irf) {
    select.leaf.id <- rforest$tree.info$prediction == as.numeric(class.id) + 1
  } else if (is.null(rit.param$class.cut)) {
    select.leaf.id <- rep(TRUE, nrow(rforest$tree.info))
  } else {
    select.leaf.id <- rforest$tree.info$prediction > rit.param$class.cut
  }
  
  rforest <- subsetReadForest(rforest, select.leaf.id)
  nf <- rforest$node.feature
  if (wt.pred.accuracy) {
    wt <- rforest$tree.info$size.node * rforest$tree.info$dec.purity
  } else {
    wt <- rforest$tree.info$size.node
  }
  rm(rforest)
  
  if (sum(select.leaf.id) < 2){
    return(character(0))
  } else {
    # group features if specified
    if (!is.null(varnames.grp)) nf <- groupFeature(nf, grp=varnames.grp)
    
    # drop feature if cutoff.unimp.feature is specified
    if (cutoff.unimp.feature > 0){
      if (!class.irf)
        rfimp <- rf$importance[,'IncNodePurity']
      else
        rfimp <- rf$importance[,'MeanDecreaseGini']   
      drop.id <- which(rfimp < quantile(rfimp, prob=cutoff.unimp.feature))
      nf[,drop.id] <- FALSE
    }
    
    interactions <- RIT(nf, weights=wt, depth=rit.param$depth, 
                        n_trees=rit.param$ntree, branch=rit.param$nchild, 
                        n_cores=n.core)
    interactions$Interaction <- gsub(' ', '_', interactions$Interaction)
    return(interactions)
  }
}

subsetReadForest <- function(rforest, subset.idcs) {
  # Subset nodes from readforest output 
  if (!is.null(rforest$node.feature)) 
    rforest$node.feature <- rforest$node.feature[subset.idcs,]
  if(!is.null(rforest$tree.info))
    rforest$tree.info <- rforest$tree.info[subset.idcs,]
  return(rforest)
}

groupFeature <- function(node.feature, grp){
  # Group feature level data in node.feature 
  sparse.mat <- is(node.feature, 'Matrix')
  
  grp.names <- unique(grp)
  makeGroup <- function(x, g) apply(as.matrix(x[,grp == g]), MARGIN=1, max) 
  node.feature.new <- sapply(grp.names, makeGroup, x=node.feature)
  if (sparse.mat) node.feature.new <- Matrix(node.feature.new, sparse=TRUE)
  
  colnames(node.feature.new) <- grp.names
  
  return(node.feature.new)
}



summarizeInteract <- function(store.out, varnames=NULL){
  # Aggregate interactions across bootstrap samples
  n.bootstrap <- length(store.out)
  store <- do.call(rbind, store.out)

  if (length(store) >= 1){
    int.tbl <- sort(table(store$Interaction), decreasing = TRUE)
    int.tbl <- int.tbl / n.bootstrap

    prev.tbl <- c(by(store$Prevalence, store$Interaction, sum))
    prev.tbl <- prev.tbl / n.bootstrap
    prev.tbl <- prev.tbl[names(int.tbl)]
  } else {
    return(list(interaction=numeric(0), prevalence=numeric(0)))
  }
  
  if (!is.null(varnames)) {
    stopifnot (names(int.tbl) == names(prev.tbl))
    names.int <- lapply(names(int.tbl), strsplit, split='_')
    names.int <- lapply(names.int, unlist)
    names.int <- sapply(names.int, function(n) {
      nn <- as.numeric(n)
      return(paste(varnames[nn], collapse='_'))
    })
    names(int.tbl) <- names.int
    names(prev.tbl) <- names.int

  }
  out <- list(interaction=int.tbl, prevalence=prev.tbl)
  return(out)
}

sampleClass <- function(y, cl, n) {
  # Sample indices specific to a given class
  sampled <- sample(which(y == cl), n, replace=TRUE)
  return(sampled)
}
