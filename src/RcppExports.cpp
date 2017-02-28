// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nodeObs
IntegerVector nodeObs(IntegerVector obsnodes, int n, int ntree, IntegerVector nrnodes, IntegerVector nodeobs);
RcppExport SEXP iRF_nodeObs(SEXP obsnodesSEXP, SEXP nSEXP, SEXP ntreeSEXP, SEXP nrnodesSEXP, SEXP nodeobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type obsnodes(obsnodesSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nrnodes(nrnodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nodeobs(nodeobsSEXP);
    rcpp_result_gen = Rcpp::wrap(nodeObs(obsnodes, n, ntree, nrnodes, nodeobs));
    return rcpp_result_gen;
END_RCPP
}
// nodeVars
IntegerVector nodeVars(IntegerVector varnodes, int p, int ntree, int maxnodes, int nr, IntegerVector nrnodes, IntegerVector nodevars, IntegerVector idcskeep);
RcppExport SEXP iRF_nodeVars(SEXP varnodesSEXP, SEXP pSEXP, SEXP ntreeSEXP, SEXP maxnodesSEXP, SEXP nrSEXP, SEXP nrnodesSEXP, SEXP nodevarsSEXP, SEXP idcskeepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type varnodes(varnodesSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    Rcpp::traits::input_parameter< int >::type maxnodes(maxnodesSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nrnodes(nrnodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nodevars(nodevarsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idcskeep(idcskeepSEXP);
    rcpp_result_gen = Rcpp::wrap(nodeVars(varnodes, p, ntree, maxnodes, nr, nrnodes, nodevars, idcskeep));
    return rcpp_result_gen;
END_RCPP
}