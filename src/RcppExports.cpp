#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

arma::mat P3D_cpp(const arma::mat& Y, const arma::mat& M, const arma::mat& X, const arma::mat& Vinv, const int ncores);
RcppExport SEXP _ASRgwas_P3D_cpp(SEXP YSEXP, SEXP MSEXP, SEXP XSEXP, SEXP VinvSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vinv(VinvSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(P3D_cpp(Y, M, X, Vinv, ncores));
    return rcpp_result_gen;
END_RCPP
}

arma::mat schulz_cpp(const arma::mat& X, const arma::mat& Xinv_init, const arma::uvec& na, const int niter);
RcppExport SEXP _ASRgwas_schulz_cpp(SEXP XSEXP, SEXP Xinv_initSEXP, SEXP naSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xinv_init(Xinv_initSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type na(naSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(schulz_cpp(X, Xinv_init, na, niter));
    return rcpp_result_gen;
END_RCPP
}

arma::mat P3D_schulz_cpp(const arma::mat& Y, const arma::mat& M, const arma::mat& X, const arma::mat& V, const arma::mat& Vinv, const int niter, const int ncores);
RcppExport SEXP _ASRgwas_P3D_schulz_cpp(SEXP YSEXP, SEXP MSEXP, SEXP XSEXP, SEXP VSEXP, SEXP VinvSEXP, SEXP niterSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vinv(VinvSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(P3D_schulz_cpp(Y, M, X, V, Vinv, niter, ncores));
    return rcpp_result_gen;
END_RCPP
}

arma::mat woodbury_cpp(const arma::mat& X, const arma::uvec& na, const arma::uvec& nona);
RcppExport SEXP _ASRgwas_woodbury_cpp(SEXP XSEXP, SEXP naSEXP, SEXP nonaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type na(naSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type nona(nonaSEXP);
    rcpp_result_gen = Rcpp::wrap(woodbury_cpp(X, na, nona));
    return rcpp_result_gen;
END_RCPP
}

arma::mat P3D_woodbury_cpp(const arma::mat& Y, const arma::mat& M, const arma::mat& X, const arma::mat& Vinv, const int ncores);
RcppExport SEXP _ASRgwas_P3D_woodbury_cpp(SEXP YSEXP, SEXP MSEXP, SEXP XSEXP, SEXP VinvSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vinv(VinvSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(P3D_woodbury_cpp(Y, M, X, Vinv, ncores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ASRgwas_P3D_cpp", (DL_FUNC) &_ASRgwas_P3D_cpp, 5},
    {"_ASRgwas_schulz_cpp", (DL_FUNC) &_ASRgwas_schulz_cpp, 4},
    {"_ASRgwas_P3D_schulz_cpp", (DL_FUNC) &_ASRgwas_P3D_schulz_cpp, 7},
    {"_ASRgwas_woodbury_cpp", (DL_FUNC) &_ASRgwas_woodbury_cpp, 3},
    {"_ASRgwas_P3D_woodbury_cpp", (DL_FUNC) &_ASRgwas_P3D_woodbury_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_ASRgwas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
