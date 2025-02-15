// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pbinsumRow
NumericVector pbinsumRow(NumericVector y, double N, NumericVector p);
RcppExport SEXP _nmixgof_pbinsumRow(SEXP ySEXP, SEXP NSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbinsumRow(y, N, p));
    return rcpp_result_gen;
END_RCPP
}
// pbinsum
NumericMatrix pbinsum(NumericMatrix y, NumericVector N, NumericMatrix p);
RcppExport SEXP _nmixgof_pbinsum(SEXP ySEXP, SEXP NSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbinsum(y, N, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nmixgof_pbinsumRow", (DL_FUNC) &_nmixgof_pbinsumRow, 3},
    {"_nmixgof_pbinsum", (DL_FUNC) &_nmixgof_pbinsum, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_nmixgof(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
