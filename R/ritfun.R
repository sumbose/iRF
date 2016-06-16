ritfun=function(node_feature
               , tree_depth=5
               , n_child=2
               , n_ritree=500
               , varnames=NULL
               , verbose=TRUE
               , node_wt = NULL
              ){
   if (!is.matrix(node_feature))
       stop('in ritfun(), input node_feature is not a matrix !')

   if (nrow(node_feature)==1)
       stop('node_feature has only one row -- not running RIT !')
   
   S=set2int(!!node_feature)
   v=random_intersect_trees(S, D=tree_depth
                                , M=n_ritree
                                , n_child=n_child
                                , verbose=FALSE
                                , node_wt = node_wt
                                 )
   if (is.null(v)){
       if (verbose) print('no interaction detected!')
       return(character(0))
   }
   else{
   v=int2set(v, p=ncol(node_feature), label=TRUE
              , var.names=varnames)
   tt=table(v)
   tt=sort(tt, decreasing=TRUE)
   return(names(tt))
   }
}

# input: matrix with sets (binary seq) in rows
# ouput: matrix of integers with nrow(input) rows, ceiling(ncol(output)/32) columns
set2int=function(S){
   if (length(dim(S))!=2)
       stop('input to set2int() is not a matrix!')
    
   n=nrow(S); p=ncol(S)
   pad=ifelse(p%%32, 32-p%%32, 0); n_block=(p+pad)/32

   if (pad != 0)
      S=cbind(S, array(FALSE, c(n, pad)))
   store=array(as.integer(0), c(n, n_block))
   temp<-foreach(block=1:n_block, .combine=cbind, .multicombine=TRUE)%dopar%{
      tmp=packBits(t(S[,((block-1)*32+seq(32))])
                           , type='integer')
      tmp
   }
   store=as.matrix(temp)
   return(store)
}

# input: matrix with integers 
int2set=function(v, p, label=FALSE, var.names=NULL){

if (length(dim(v))!=2)
    stop('input to int2set() is not a matrix!')
    
if (is.null(var.names))
   var.names=seq(p)

   n_block=ncol(v)
   S<- foreach(block=1:n_block, .combine=cbind, .multicombine=TRUE)%dopar%{
      tmp=matrix(as.integer(intToBits(v[,block])), ncol=32, byrow=TRUE)
      tmp
   }

   # remove trailing zeros
   S=S[,1:p]

   S=(S==1)
   if (!label)
      return(S)
   else{
      if (is.null(dim(S)))
      S=as.matrix(t(S), nrow=1)
      return(apply(S, 1
                 , FUN=function(vv){
                         return(paste(sort((var.names[which(vv)])), collapse='_'))
                       }
                  )
            )
   }
}

random_intersect_trees=function(S, D=5, M=100, n_child=3, verbose=TRUE, node_wt = NULL){

   if (length(dim(S))!=2)
       stop('input to random_intersection_trees() is not a matrix!')
   
   n=nrow(S)
   n_block=ncol(S)

   if (!is.null(node_wt)){
      if (length(node_wt)!= nrow(S)){
         warning('length of node weight vector does not equal number of nodes -- importance sampling ignored!!')
         node_wt = NULL
      }
   }

   
   out=foreach(m=1:M, .combine=rbind, .multicombine=TRUE)%dopar%{
      tmp=S[sample(n, 1, prob=node_wt),]
      tmp=as.matrix(t(tmp), nrow=1)

      for (d in 1:D){
         if (is.null(tmp)|(nrow(tmp)==0))
            next
         j_d=nrow(tmp)
         store=matrix(0, nrow=1, ncol=n_block)
         for (j in 1:j_d){
            select_samples=as.matrix(S[sample(n, n_child, replace=FALSE, prob = node_wt),]); 
            for (block in 1:n_block){
               select_samples[,block]=bitwAnd(tmp[j,block], select_samples[,block])
            }
               store=rbind(store, select_samples)
         }# end for (j in 1:j_d)
         store=as.matrix(store[-1,])

         na_id=apply(is.na(store), 1, FUN=max)
         blank_id=!apply(!!store, 1, FUN=max)
         keep_id=!(na_id|blank_id)
         tmp=matrix(store[keep_id,], nrow=sum(keep_id))

      }# end for (d in 1:D)
      if (is.null(tmp)|(nrow(tmp)==0)){
         if (verbose){
            print(paste('tree ', m, ': no set survived'))
         }
         tmp=matrix(NA, nrow=1, ncol=n_block)
      }
      tmp
   }# end foreach(m=1:M)
   
   na_id=apply(is.na(out), 1, FUN=max)
   blank_id=!apply(!!out, 1, FUN=max)

   idx = sum(!(na_id|blank_id))
   if (idx == 0){
       out=NULL
   }
   else if (idx ==1){
       out = matrix(out[!(na_id|blank_id),], nrow=1)
   }
   else{
       out=as.matrix(out[!(na_id|blank_id),])
   }
return(out)
}

groupFeature=function(node_feature
                    , grp
                     ){
   grp_names = unique(grp)
   node_feature_new = array(FALSE, c(nrow(node_feature), length(grp_names)))
   for (g in grp_names){
     grp_id = which(grp== g)
     if (length(grp_id)==1){
         node_feature_new[,grp_names==g] = node_feature[,grp_id]
     }
     else{
         node_feature_new[,grp_names==g] = apply(node_feature[,grp_id], 1, FUN=max)
     }
   }
node_feature_new = !!node_feature_new
colnames(node_feature_new) <- grp_names

return(node_feature_new)
}



summarize_interact = function(store_out, varnames = NULL){
n_bootstrap = length(store_out)

store = character()

for (b in 1:n_bootstrap){
    tmp = store_out[[b]]
    if ((length(tmp)>0))
       store = c(store, tmp)
    
}

if (length(store) >= 1){
    store_tbl= sort(table(store), decreasing = TRUE)
    tt = as.vector(store_tbl)
    names(tt) = names(store_tbl)
    tt = tt/n_bootstrap
}
else{
    return(character(0))
}

if (!is.null(varnames)){
    for (j in 1:length(tt)){
        tmp = as.numeric(strsplit(names(tt)[j], split='_')[[1]])
        names(tt)[j] = paste(varnames[tmp], collapse='_')
    }
}
return(tt)
}
