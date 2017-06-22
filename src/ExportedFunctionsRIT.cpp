#include<vector>
#include<algorithm>
#include<set>
#include<Rcpp.h>
#include "RaggedArray.h"
#include "InputFunctions.h"
#include "RITmain.h"
#include "RITaux.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP RIT_1class(SEXP z, NumericVector weights, int L, int branch, int depth, int n_trees, int min_inter_sz, int n_cores, bool is_sparse) {
    //construct RaggedArray
    RaggedArray x;
    if (is_sparse) {
      List iplist(z);
      IntegerVector i_sparse(iplist[0]);
      IntegerVector p_sparse(iplist[1]);
      x = InputSparseMatrix(i_sparse, p_sparse);
    }
    else {
      x = InputLogicalMatrix(z);
    }
    
    //n and p are the number of rows and cols of z respectively
    int n=x.nrow();
    int p=x.ncol();
    
    //perform RIT to find total_candidate_interactions
    set<vector<int> > total_candidate_interactions=RIT_basic(x, weights, L,branch,depth,n_trees,min_inter_sz, n_cores,n);
    
    //Get vector of prevalences giving the estimated prevalence of each candidate interaction
    //For this we need to create the minhash matrix Ht
    int** Ht;
    Ht = new int* [p];
    for (int k=0; k<p; k++) {
      Ht[k] = new int [L];
    }

    CreateHt(x, L, Ht);
    
    const double n_plus1_over_n=(n+1)/n; 
    const double recip_n_plus_1=1/(n+1);
    vector<double> prevalences=PrevEst_inter(total_candidate_interactions, Ht, L,n_plus1_over_n,recip_n_plus_1);
    
    //return list 'Output', containing both total_candidate_interactions and prevalences
    List output;
    output["Interactions"]=AddOne(total_candidate_interactions);
    output["Prevalence"]=prevalences;
    
    // Delete Ht
    for (int k=0; k<p; k++) {
		  delete[] Ht[k];
	  }
	  delete[] Ht;
    
    //return output
    return output;
}

// [[Rcpp::export]]
SEXP RIT_2class(SEXP z, SEXP z0, int L, int branch, int depth, int n_trees, double theta0, double theta1,
  int min_inter_sz, int n_cores, bool is_sparse) {
  //construct RaggedArrays
  RaggedArray x; RaggedArray x0;
  if (is_sparse) {
    // extract i_Sparse,p_Sparse & i_Sparse0,p_Sparse0 components from z,z0 respectively
    List iplist(z); List iplist0(z0);
    IntegerVector i_sparse(iplist[0]); IntegerVector i_sparse0(iplist0[0]);
    IntegerVector p_sparse(iplist[1]); IntegerVector p_sparse0(iplist0[1]);
        
    x = InputSparseMatrix(i_sparse, p_sparse);
    x0 = InputSparseMatrix(i_sparse0, p_sparse0);
  } 
  else {
    x = InputLogicalMatrix(z);
    x0 = InputLogicalMatrix(z0);
  }
  //n,n0 and p are the number of rows and cols of z,z0 respectively
  int n=x.nrow(); int n0=x0.nrow();
  int p=max(x.ncol(),x0.ncol());
  
  // Setup
  const double n_plus1_over_n=(n+1)/n; const double n0_plus1_over_n0=(n0+1)/n0;
  const double recip_n_plus_1=1/(n+1); const double recip_n0_plus_1=1/(n0+1);
  
  // Create H matrices
  int** Ht; int** H0t;
  Ht = new int* [p]; H0t=new int* [p];
  for (int k=0; k<p; k++) {
      Ht[k] = new int [L];
      H0t[k] = new int [L];
	}
  CreateHt(x, L, Ht); CreateHt(x0,L,H0t);
  
  // Perform RIT
  // Sets prevalent in z but not in z0
  set<vector<int> > total_candidate_interactions=RIT_minhash(x, L, branch, depth, n_trees, 
    theta0, theta1,min_inter_sz, n_cores, n, H0t, n0_plus1_over_n0, recip_n0_plus_1);
  // Sets prevalent in in z0 but not in z
  set<vector<int> > total_candidate_interactions0=RIT_minhash(x0, L, branch, depth, n_trees, 
    theta0, theta1, min_inter_sz, n_cores, n, Ht, n_plus1_over_n, recip_n_plus_1);

  vector<double> prevalences=PrevEst_inter(total_candidate_interactions, Ht, L,n_plus1_over_n, recip_n_plus_1);
  vector<double> prevalences0=PrevEst_inter(total_candidate_interactions0, H0t, L,n0_plus1_over_n0, recip_n0_plus_1);
  
  //create lists 'output','output0' containing both total_candidate_interactions and prevalences for each class
  List output; List output0;
  output["Interactions"]=AddOne(total_candidate_interactions);
  output0["Interactions"]=AddOne(total_candidate_interactions0);
  output["Prevalence"]=prevalences; output0["Prevalence"]=prevalences0;
  List combined_output;
  combined_output["Class1"]=output;
  combined_output["Class0"]=output0;
    
  // Delete H Matrices
  for (int k=0; k<p; k++) {
    delete[] Ht[k]; delete[] H0t[k];
	 }
	 delete[] Ht; delete[] H0t;
  
  //return combined_output
  return combined_output;
}
