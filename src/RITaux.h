#ifndef RITaux
#define RITaux

#include<vector>
#include<algorithm>
#include<random>
#include<set>
#include<math.h>
#include<Rcpp.h>
#include "RaggedArray.h"

using namespace std;

typedef vector<int>::iterator v_iterator;

// compute intersection between two vectors a and b where size(a)>>size(b)
// for each element of b, we see if it is in a using binary search
// we assume elements of a,b are sorted in ascending order
vector<int> binary_intersect(v_iterator large_begin, v_iterator large_end,
    v_iterator small_begin, v_iterator small_end) {
  vector<int> intersection;
  v_iterator it;
  for (it = small_begin; it != small_end; it++) {
    if (binary_search(large_begin, large_end, *it)) {
      intersection.push_back(*it);
    }
  }
  return intersection;
}

void CreateHt(RaggedArray &x, const int L, int** Ht) {
  // Ht is p by L
  random_device rd; //seed for Random Number Generator(RNG)
  mt19937_64 mt(rd()); //Use Mersenne Twister as RNG
  
  const int n=x.nrow();
  const int p=x.ncol();
  vector<int> perm(n); //vector perm_i=i
  for (int k=0; k<n; k++) {
    perm[k] = k;
  }
  for (int l=0; l<L; l++) {
    shuffle(perm.begin(), perm.end(), mt);
    for (int k=0;k<p;k++){
      bool T=false;
      int i=0;
      while (!T && (i<n)) {
        T=binary_search(x.begin(perm[i]),x.end(perm[i]),k);
        i++;
      }
      if (i==n) {
        Ht[k][l]=0;
      }
      else {
        Ht[k][l]=i;
      }
    }
  }
}

// Computes P1(L,intersection,H)=1/L{#l st H[l][k]=H[l][k'] for all k,k' in intersection}
// Computes P2(L,intersection,H)={(n+1)/n}*{1/(average over l(min(H[l][k]))) - 1/(n+1)}
double PrevEst(const vector<int> &intersection, int** Ht, const int L, const double n_plus_1_over_n, const double recip_n_plus_1) {
  if (intersection.size() <= 1) return 1; // if the interaction has size 1, return 1
  int P1 = 0;
  double P2 = 0;
  int T;
  int m;
  for (int l=0; l<L; l++){
  	T = 1;
    m = Ht[intersection[0]][l];
		for (unsigned int q=1; q<intersection.size(); q++) {
			if (m != Ht[intersection[q]][l]) {
				T = 0;
        m = min(m, Ht[intersection[q]][l]);
			}
      if ((m == 1) && (T == 0)) break;
		}    
		P1 += T;
    P2 += double(m);
	}
  return n_plus_1_over_n * (double(L)/P2-recip_n_plus_1) * double(P1)/double(L);
}

// Given a set of interactions (integer vectors) copmutes their prevalence using the minwise hash matrix
// Returns a vector of prevalences
vector<double> PrevEst_inter(const set<vector<int> > &interactions, int** Ht, const int L, const double n_plus_1_over_n, const double recip_n_plus_1) {
  vector<double> prevalences(interactions.size());
  set<vector<int> >::const_iterator it;
  int j=0;
  for (it=interactions.begin(); it!=interactions.end(); it++) {
    prevalences[j]=PrevEst(*it, Ht, L, n_plus_1_over_n, recip_n_plus_1);
    j++;
  }
  return prevalences;
}

// Given a set of interactions (integer vectors) this produces a set with 1 added to each entry
// Note function gives an error when we pass &interactions instead and try to modify interactions
Rcpp::List AddOne(const set<vector<int> > &interactions) {
  Rcpp::List interactions_plus1(interactions.size());
  set<vector<int> >::iterator it;
  int i=0;
  for (it=interactions.begin(); it!=interactions.end(); it++) {
    Rcpp::IntegerVector interaction = wrap(*it);
    interactions_plus1[i] = interaction + 1; // uses Rcpp sugar
    i++;
  }
  return interactions_plus1;
}

#endif
