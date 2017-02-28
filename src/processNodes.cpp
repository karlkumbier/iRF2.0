#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector nodeObs(IntegerVector obsnodes, 
                      int n, int ntree, IntegerVector nrnodes, 
                      IntegerVector nodeobs) {

/* Read data from RF with node tracking and return nleaf x n binary
 * matrix indicating which leaf nodes each observation falls in */  
  int cumnodes = 0; //number of processed nodes
  
  for (int t = 0; t < ntree; t++) {
    
    int nrnodest = nrnodes[t];
    
    for (int nd = 0; nd < nrnodest; nd ++) {
      
      for (int i = 0; i < n; i ++) {
        
        if (obsnodes[i + t * n] == -1) continue;       
          nodeobs[i + (obsnodes[i + t * n]) * n + n * cumnodes] = 1;
      }
    }
    cumnodes += nrnodest;
  }
  
  return nodeobs;
}


// [[Rcpp::export]]
IntegerVector nodeVars(IntegerVector varnodes, 
                       int p, int ntree, int maxnodes, int nr,
                       IntegerVector nrnodes, 
                       IntegerVector nodevars,
                       IntegerVector idcskeep) {
 
/* Read data from RF with node tracking and return length(idcskeep) x p binary
 * matrix indicating decision path features at each node, the rows
 * of this matrix will be subset based on nodes specified by idcskeep*/ 
  int cumnodes = 0; //number of processed nodes
  int currentnd = 0; //current node across all trees
  
  for (int t = 0; t < ntree; t++) {
    
    int nrnodest = nrnodes[t]; 
    int j = 0;
    
    // create dense matrix indicating which observations are at each node for
    //  the specified subset of nodes
    for (int nd = 0; nd < nrnodest; nd ++) {
      for (int i = 0; i < p; i ++) {
        if (varnodes[j + t * nr] == (i + p * nd) &&
            idcskeep[currentnd] - 1 == nd + cumnodes) {
          nodevars[i + p * currentnd] = 1;
          j += 1;
        } else if (varnodes[j + t * nr] == (i + p * nd)) {
          j += 1;
        }
      }
      if (idcskeep[currentnd] - 1 == nd + cumnodes) currentnd += 1;
    }
    cumnodes += nrnodest;
  }
  
  return nodevars;
}

