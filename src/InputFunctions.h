#ifndef InputFunctions_inc
#define InputFunctions_inc

#include<Rcpp.h>
#include "RaggedArray.h"

using namespace Rcpp;

RaggedArray InputLogicalMatrix(LogicalMatrix z) {
  RaggedArray x;
  for (int i=0; i<z.nrow(); i++) {
    x.new_row();
    for (int j=0; j<z.ncol(); j++) {
      if (z(i,j)) {
        x.push_back(j);
      }
    }
  }
  return x;
}

RaggedArray InputSparseMatrix(IntegerVector i_vec, IntegerVector p_vec){
  //i and p are @i and @p of the transpose of the sparseMatrix
  RaggedArray x;
  for (int p_it=1; p_it<p_vec.size(); p_it++) {
    x.new_row();
    for (int j=p_vec[p_it-1]; j<p_vec[p_it]; j++) {
      x.push_back(i_vec[j]);
    }
  }
  return x;
}
#endif
