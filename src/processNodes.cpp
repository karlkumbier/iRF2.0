#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector nodeObs(IntegerVector obsnodes, 
    int n, int ntree, 
    IntegerVector nrnodes, 
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
    int ntree, 
    int nrnodes,
    int p,
    IntegerVector parents,
    IntegerVector idcskeep,                  
    IntegerVector nodect,
    IntegerVector nnodest,
    IntegerVector nodevars) {

  int cumnodes = 0; //cumulative number of nodes processed over forest
  int idx = 0; //current index of output vector
  int nobsnode = 0; //number of observations in current node
  int annd = 0; //ancestor nodes
  int pp = 0; //track current col of sparse matrix
  int rowoffset = 0; //track current row of sparse matrix

  for (int t = 0; t < ntree; t++) {

    for (int nd = 0; nd < nnodest[t]; nd++) {

      if (idcskeep[nd + cumnodes] == 1) {
        // Keeping the current node. Determine variables on path and store
        // entries in output vector based on counts for associated node.

        nobsnode = nodect[nd + cumnodes];
        annd = nd;
        IntegerVector pathvars(p);
        while (annd > 0) {  
          annd = parents[annd + cumnodes] - 1;
          pp = varnodes[annd + nrnodes * t];

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
    cumnodes += nnodest[t];
  }

  return nodevars;
}


/*
// [[Rcpp::export]]
IntegerVector nodeVars(IntegerVector varnodes,
                       int p, int ntree, int nnodes,
                       IntegerVector idcskeep,
                       IntegerVector nodect,
                       IntegerVector nnodest,
                       IntegerVector nodevars) {


  int rowoffset = 0;
  int idx = 0;
  int cumnodes = 0;
  int rw, cl, nxtrw, nxtcl;
  int nobsnode;
  for (int t = 0; t < ntree; t++) {
      for (int i = 0; i < nnodes; i ++) {

        rw = varnodes[i + nnodes * t] % p;
        cl = (varnodes[i + nnodes * t] - rw) / p;

        if (idcskeep[cl + cumnodes] == 1) {

          nobsnode = nodect[cl + cumnodes];
          for (int j = 0; j < nobsnode; j ++) {
            nodevars[idx] = j + rowoffset + 1;
            nodevars[nodevars.size() / 2 + idx] = rw + 1;
            idx += 1;
          }

          nxtrw = varnodes[1 + i + nnodes * t] % p;
          nxtcl = (varnodes[1 + i + nnodes * t] - nxtrw) / p;
          if (nxtcl != cl) rowoffset += nobsnode;
        }
      }
      cumnodes += nnodest[t];
  }

  return nodevars;
  }

// [[Rcpp::export]]
IntegerVector nodeVars(IntegerVector varnodes, 
                       int p, int ntree, int nr,
                       IntegerVector nodect,
                       IntegerVector nrnodes, 
                       IntegerVector nodevars,
                       IntegerVector idcskeep) {
 
 Read data from RF with node tracking and return length(idcskeep) x p binary
  matrix indicating decision path features at each node, the rows
 of this matrix will be subset based on nodes specified by idcskeep
  int cumnodes = 0; //number of processed nodes
  int currentnd = 0; //current node across all trees
  
  for (int t = 0; t < ntree; t++) {
    Rprintf("tree %d\n", t);     
    int ii = 0;
    int nrnodest = nrnodes[t];
    IntegerVector nodectt(nrnodes);
    for (int i = cumnodes; i < (cumnodes + nrnodest); i++) {
      nodectt[ii] = nodect[i];
      ii += 1;
    } 

    int j = 0;
    
    // create dense matrix indicating which observations are at each node for
    //  the specified subset of nodes
    for (int nd = 0; nd < nrnodest; nd ++) {
      for (int i = 0; i < p; i ++) {
        if (varnodes[j + t * nr] == (i + p * nd) &&
            idcskeep[currentnd] - 1 == nd + cumnodes) {
          nodevars[i + p * currentnd] = 1;
          // replicate based on counts in current node
          //Rprintf("current node %d\n", currentnd);
          //Rprintf("node count %d\n", nodectt[nd]);
          //for (int colidx = currentnd; colidx < (currentnd + nodectt[nd]); colidx++) {
          //  nodevars[i + p * colidx] = 1;   
          //} 
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
*/
