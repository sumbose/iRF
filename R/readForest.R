#   PASS DATA THROUGH FOREST
readForest <- function(rfobj  # a randomForest object with forest component in it
                     , X   # n x p data matrix 
                     , Y=NULL   # label vector (numeric, 0/1)
                     , return_node_feature=TRUE
                     , return_node_data=TRUE
                     , leaf_node_only = TRUE
                     , verbose = TRUE
                       ){
if (is.null(rfobj$forest))
    stop('No Forest component in the randomForest object')

if (!is.null(Y)){
  if (!is.numeric(Y)){
    stop("Y must be a numeric vector")
  }
  else{
     tmp = unique(Y)
     if (!setequal(tmp, c(0,1)))
       stop('Y must be 0-1 valued')
  }

  if (length(Y)!=nrow(X))
      stop("Number of observations in X, Y do not match!")
}

n = nrow(X)
p = ncol(X)
ncol_tree_info = ifelse(is.null(Y), 9, 12)

store<-foreach(k=1:rfobj$ntree, .combine=rbind, .multicombine=TRUE, .packages='iRF')%dopar%{
   if (verbose==TRUE & (k%%50==0)) {print(paste('tree:', k))}
   tt = getTree(rfobj, k=k, labelVar=FALSE); rownames(tt) <- NULL
   tt_out=readTree(rfTree = tt, X = X, Y = Y
                  , return_node_feature=return_node_feature
                  , return_node_data=return_node_data
                  , leaf_node_only=leaf_node_only)
   
   tt_out$tree_info$tree_id=k

   tmp=as.matrix(tt_out$tree_info)
   if (return_node_feature){
   colnames(tt_out$node_feature)<- seq(ncol(tt_out$node_feature))
   tmp=cbind(tmp, tt_out$node_feature)
   }
   if (return_node_data){
   colnames(tt_out$node_data) <- seq(ncol(tt_out$node_data))
   tmp=cbind(tmp, tt_out$node_data)
   }
   tmp
}
rownames(store) <- NULL

# account for new column: tree_id
ncol_tree_info = ncol_tree_info+1

out=list()
out$tree_info =
as.data.frame(store[,1:ncol_tree_info]
            , stringsAsFactors=FALSE
            , check.names=FALSE
             )

if (return_node_feature)
out$node_feature = store[,(ncol_tree_info+seq(p))]
if (return_node_data)
out$node_data = store[,(ncol_tree_info+p*(return_node_feature)+seq(n))]

return(out)
}



#######################################
#  PASS DATA THROUGH A TREE
#######################################
readTree <- function(rfTree # a single tree in the forest, as in the output from getTree() function
                 , X  # n x p matrix of features
                 , Y=NULL  # numeric (0/1) p-vector of labels (NOT as factors)
                 , return_node_feature = TRUE
                 , return_node_data = TRUE
                 , leaf_node_only = TRUE
                  ){
   tt=as.data.frame(rfTree, check.names=FALSE, stringsAsFactors=FALSE)
   
   n=nrow(X)
   p=ncol(X)
   n_node=nrow(tt)
   
   node_face=array(FALSE, c(n_node, p))
   node_composition=array(FALSE, c(n_node, n))
      
   tt$parent=0
   tt$size_node=0
   tt$depth=0

   if (!is.null(Y)){
     tt$size_0=0
     tt$Gini=0
     tt$Decrease_Gini=0
   }
   
   # add parent information
   leaf_id=(tt$status==-1)
   for (i in which(!leaf_id)){
      tt$parent[tt$"left daughter"[i]]=i
      tt$parent[tt$"right daughter"[i]]=i
   }
   
   # initialize
   node_composition[1,]=TRUE
   tt$size_node[1]=n
   tt$depth[1]=0   

   if (return_node_feature)
   node_face[1,]=FALSE
   
   if (!is.null(Y)){
   tt$size_0[1]=sum(Y==0)
   tt$Gini[1]=2*tt$size_0[1]*(n-tt$size_0[1])/(n^2)
   }

   # run  over all non-leaf nodes 
   for (i in which(!leaf_id)){

      d_left=tt$"left daughter"[i]
      d_right=tt$"right daughter"[i]
      
      split_var=tt$"split var"[i]
      split_pt=tt$"split point"[i]

      parent_id=node_composition[i,]
      d_left_id=(X[,split_var] < split_pt) & parent_id
      d_right_id=(X[,split_var] >= split_pt) & parent_id

      if (return_node_feature){
      node_face[d_left,]=node_face[i,]
      node_face[d_right,]=node_face[i,]
      node_face[c(d_left, d_right), split_var]=TRUE
      }
      

      node_composition[d_left,]=d_left_id
      node_composition[d_right,]=d_right_id

      
      tt$size_node[d_left]=sum(d_left_id)
      tt$size_node[d_right]=sum(d_right_id)
      
      tt$depth[d_left]=tt$depth[i]+1
      tt$depth[d_right]=tt$depth[i]+1
      
  if (!is.null(Y)){
      tt$size_0[d_left]=sum(Y[d_left_id]==0)
      tt$size_0[d_right]=sum(Y[d_right_id]==0)
      
      tt$Gini[d_left]=2*tt$size_0[d_left]*
         (tt$size_node[d_left]-tt$size_0[d_left])/
         (tt$size_node[d_left])^2
      tt$Gini[d_right]=2*tt$size_0[d_right]*
         (tt$size_node[d_right]-tt$size_0[d_right])/
         (tt$size_node[d_right])^2
      
      tt$Decrease_Gini[i]=tt$Gini[i]-
         (tt$size_node[d_left]*tt$Gini[d_left]+
         tt$size_node[d_right]*tt$Gini[d_right])/
         tt$size_node[i]
  }

  }

select_node = 1:nrow(tt)
if (leaf_node_only)
  select_node = tt$status == -1

out=list()
out$tree_info=tt[select_node,]
if (return_node_feature)
    out$node_feature = node_face[select_node,]
if (return_node_data)
    out$node_data = node_composition[select_node,]

return(out)
}

