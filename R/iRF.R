# iteratively grows random forests, finds case specific feature interactions
# in binary classification problems

iRF <- function(x
              , y   # assume y, ytest are factors with level 0-1
              , xtest = NULL
              , ytest = NULL
              , n_iter = 5
              , ntree = 500
              , n_core = 1
              , mtry_select_prob = rep(1/ncol(x), ncol(x))
              , keep_impvar_quantile = NULL 
              , find_interaction = FALSE
              , cutoff_nodesize_prop = 0.1
              , cutoff_unimp_feature = 0
              , class_id = 1
              , rit_param = c(5, 1000, 2) # c(rit depth (D), n_rit (M), n_child)
              , varnames_grp = NULL
              , n_bootstrap = 30
              , verbose = TRUE
              , ...
               ){

n = nrow(X)
p = ncol(X)
keep_subset_var = NULL

# check if y, ytest are factors with levels 0,1
if (!is.factor(y))
    stop('y is not a factor')
if (!setequal(as.numeric(levels(y)), c(0,1)))
    stop('y must have levels 0, 1')

if ((!is.null(ytest)) & (!is.factor(ytest)))
    stop('ytest is not a factor')
if ((!is.null(ytest)) & (!setequal(as.numeric(levels(ytest)), c(0,1))) )
    stop('ytest must have levels 0, 1')

# check if x is a numeric matrix with two or more columns if find_interaction = TRUE
if (!is.matrix(x) | ((!is.null(xtest)) & (!is.matrix(xtest))))
    stop('either x or xtest is not a matrix !')
if ((!is.numeric(x)) | ((!is.null(xtest) & (!is.numeric(xtest)))))
    stop('either x or xtest is not a numeric matrix!')
if ((ncol(x) < 2) & find_interaction == TRUE)
    stop('cannot find interaction - X has less than two columns!')

    
# initialize outputs
rf_list = list()
if (find_interaction){
    interact_list = list()
    Stability_Score = list()
}

# set number of trees to grow in each core
a = floor(ntree/n_core); b = ntree %% n_core
ntree_id = c(rep(a+1, b), rep(a, n_core-b))


for (iter in 1:n_iter){
## 1: run random forest
   print(paste('iteration = ', iter))

   rf <- foreach(nt=ntree_id, .combine=combine
               , .multicombine=TRUE, .packages='iRF')%dopar%{
       randomForest(x
                  , y
                  , xtest
                  , ytest
                  , ntree = nt
                  , mtry_select_prob = mtry_select_prob
                  , keep_subset_var = keep_subset_var
                  , keep.forest=TRUE
                  , ...
                   )
   }

   rf_list[[iter]] = rf


## 2: find interactions stable over multiple bootstrap replicates
   if (find_interaction == TRUE){
      if (verbose){cat('finding interactions ... ')}

      interact_list[[iter]] = list()
      Stability_Score[[iter]] = list()
      
      n1 = sum(y==1); n0 = sum(y==0)

      # 2.1: find interactions in  multiple bootstrap samples to assess stability
      for (i_b in 1:n_bootstrap){
         if (verbose){cat(paste('b = ', i_b, ';  ', '\n', sep=''))}

         sample_id = c(sample(which(y==0), n0, replace=TRUE)
                     , sample(which(y==1), n1, replace=TRUE)
                      )
         #2.1.1: fit random forest
         rf_b <- foreach(nt=ntree_id, .combine=combine
                     , .multicombine=TRUE, .packages='iRF')%dopar%{
               randomForest(x[sample_id,]
                          , y[sample_id]
                          , xtest
                          , ytest
                          , ntree=nt
                          , mtry_select_prob = mtry_select_prob
                          , keep_subset_var = keep_subset_var
                          , keep.forest=TRUE
                          , ...
                           )
          }
        # 2.1.2: select large leaf nodes
        rforest_b = readForest(rf_b, X = x[sample_id,]
                             , return_node_feature = TRUE
                             , return_node_data = FALSE
                             , leaf_node_only = TRUE
                             , verbose = verbose
                              )

        cutoff_nodesize = cutoff_nodesize_prop*ifelse(class_id == 1, n1, n0)
        select_leaf_id = rforest_b$tree_info$size_node > cutoff_nodesize &
                         rforest_b$tree_info$prediction == (as.numeric(class_id)+1)
                         # since tree predictions are saved as 1, 2 for binary classification

        # 2.1.3: apply random intersection trees
        nf = rforest_b$node_feature[select_leaf_id,]

        if (!sum(select_leaf_id)){
            if (verbose) print('no large node - not running RIT!')
            interact_list[[iter]][[i_b]] = character(0)
        }
        else if (sum(select_leaf_id) == 1){
            if (verbose) print('only one large node - not running RIT!')
            interact_list[[iter]][[i_b]] = paste(which(as.vector(nf)!=0), collapse='_')
        }
        else{
        
        # grp features, if specified
        if (!is.null(varnames_grp))
            nf = groupFeature(nf, grp = varnames_grp)

        # drop features, if cutoff_unimp_feature is specified
        if (cutoff_unimp_feature > 0){
            if (ncol(rf_b$importance) == 1)
                rfimp = rf_b$importance
            else
                rfimp = rf_b$importance[,4]
            
            drop_id = which(rfimp < quantile(rfimp, prob=cutoff_unimp_feature))
            nf[,drop_id] = FALSE
        }
            
        interact_list[[iter]][[i_b]] =
                ritfun(node_feature = nf
                     , tree_depth = rit_param[1]
                     , n_ritree = rit_param[2]
                     , n_child = rit_param[3]
                     , varnames = NULL
                     , verbose = verbose
                      )
        }


      } # end for (i_b in ... )

   # 2.2: calculate stability scores of interactions
   if (!is.null(varnames_grp))
       vnames_new = unique(varnames_grp)
   else if (!is.null(colnames(x)))
       vnames_new = colnames(x)
   else
       vnames_new = seq(ncol(x))

   Stability_Score[[iter]] =
   summarize_interact(interact_list[[iter]]
                    , varnames = vnames_new)
      
   } # end if (find_interaction)
   
## 3: change mtry_select_prob (and keep_subset_var, if applicable) for next iteration
   if (ncol(rf$importance) == 1)
       mtry_select_prob = rf$importance
   else
       mtry_select_prob = rf$importance[,4]

   if (!is.null(keep_impvar_quantile)){
       keep_subset_var = which(mtry_select_prob >
                               quantile(mtry_select_prob, keep_impvar_quantile)
                              )
   }
    
   if (!is.null(xtest)){
       auroc = auc(roc(rf$test$votes[,2], ytest))
       print(paste('AUROC: ', round(auroc, 2)))
   }

} # end for (iter in ... )


out=list()
out$rf_list = rf_list
if (find_interaction == TRUE){
    out$interaction = Stability_Score
}
    
return(out)
}
