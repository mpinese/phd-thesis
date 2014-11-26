// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// snmfl_atomic
Rcpp::List snmfl_atomic(const arma::mat& A, arma::uword k, double eta, double alpha, size_t max_iter, double conv_delta, size_t conv_interval);
RcppExport SEXP snmfl_snmfl_atomic(SEXP ASEXP, SEXP kSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP max_iterSEXP, SEXP conv_deltaSEXP, SEXP conv_intervalSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP );
        Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP );
        Rcpp::traits::input_parameter< double >::type eta(etaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< size_t >::type max_iter(max_iterSEXP );
        Rcpp::traits::input_parameter< double >::type conv_delta(conv_deltaSEXP );
        Rcpp::traits::input_parameter< size_t >::type conv_interval(conv_intervalSEXP );
        Rcpp::List __result = snmfl_atomic(A, k, eta, alpha, max_iter, conv_delta, conv_interval);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
