#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector nodeVars(IntegerVector varnodes,
                       int nrnodes,
                       int p,
                       IntegerVector parents,
                       IntegerVector idcskeep,                  
                       IntegerVector nodect,
                       int rowoffset,
                       IntegerVector nodevars) {
  
  int idx = 0; //current index of output vector
  int nobsnode = 0; //number of observations in current node
  int annd = 0; //ancestor nodes
  int pp = 0; //track current col of sparse matrix

  
  for (int nd = 0; nd < nrnodes; nd++) {
    

    if (idcskeep[nd] == 1) {
      // Keeping the current node. Determine variables on path and store
      // entries in output vector based on counts for associated node.
      
      nobsnode = nodect[nd];
      annd = nd;
      IntegerVector pathvars(p);
      while (annd > 0) {  
        annd = parents[annd] - 1; //ancestor node
        pp = varnodes[annd]; // selected variable for ancestor
        
        if (pathvars[pp - 1] == 0) {
          for (int rw = 0; rw < nobsnode; rw ++) {
            nodevars[idx] = rw + rowoffset + 1;
            nodevars[nodevars.size() / 2 + idx] = pp;
            idx += 1;
          }
          pathvars[pp - 1] = 1;
        }
      }
      rowoffset += nobsnode;
    }
    
  }
  
  return nodevars;
}


