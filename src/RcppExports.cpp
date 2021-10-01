// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getPairDist
arma::mat getPairDist(const arma::mat x);
RcppExport SEXP _DR_SC_getPairDist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(getPairDist(x));
    return rcpp_result_gen;
END_RCPP
}
// calYenergy2D_sp
arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta);
RcppExport SEXP _DR_SC_calYenergy2D_sp(SEXP ySEXP, SEXP AdjSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(calYenergy2D_sp(y, Adj, K, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// obj_beta
double obj_beta(const arma::ivec& y, const arma::mat& R, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta);
RcppExport SEXP _DR_SC_obj_beta(SEXP ySEXP, SEXP RSEXP, SEXP AdjSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_beta(y, R, Adj, K, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// runICM_sp
Rcpp::List runICM_sp(const arma::mat& X, arma::ivec& y, const arma::mat& W0, const arma::vec& Lam_vec0, const arma::mat& Mu0, const arma::cube& Sigma0, const arma::sp_mat& Adj, const arma::vec& alpha, const arma::vec& beta_grid, double beta, int maxIter_ICM);
RcppExport SEXP _DR_SC_runICM_sp(SEXP XSEXP, SEXP ySEXP, SEXP W0SEXP, SEXP Lam_vec0SEXP, SEXP Mu0SEXP, SEXP Sigma0SEXP, SEXP AdjSEXP, SEXP alphaSEXP, SEXP beta_gridSEXP, SEXP betaSEXP, SEXP maxIter_ICMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W0(W0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lam_vec0(Lam_vec0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu0(Mu0SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Sigma0(Sigma0SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_grid(beta_gridSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter_ICM(maxIter_ICMSEXP);
    rcpp_result_gen = Rcpp::wrap(runICM_sp(X, y, W0, Lam_vec0, Mu0, Sigma0, Adj, alpha, beta_grid, beta, maxIter_ICM));
    return rcpp_result_gen;
END_RCPP
}
// icmem_heterCpp
Rcpp:: List icmem_heterCpp(const arma::mat& X, const arma::sp_mat& Adj, const arma::ivec& y_int, const arma::mat& Mu_int, const arma::mat& W_int, const arma::cube& Sigma_int, const arma::vec& Lam_vec_int, const arma::vec& alpha, const double& beta_int, const arma::vec& beta_grid, const int& maxIter_ICM, const int& maxIter, const double& epsLogLik, const int& verbose, const bool& homo, const bool& diagSigmak);
RcppExport SEXP _DR_SC_icmem_heterCpp(SEXP XSEXP, SEXP AdjSEXP, SEXP y_intSEXP, SEXP Mu_intSEXP, SEXP W_intSEXP, SEXP Sigma_intSEXP, SEXP Lam_vec_intSEXP, SEXP alphaSEXP, SEXP beta_intSEXP, SEXP beta_gridSEXP, SEXP maxIter_ICMSEXP, SEXP maxIterSEXP, SEXP epsLogLikSEXP, SEXP verboseSEXP, SEXP homoSEXP, SEXP diagSigmakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type y_int(y_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu_int(Mu_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W_int(W_intSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Sigma_int(Sigma_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lam_vec_int(Lam_vec_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type beta_int(beta_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_grid(beta_gridSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter_ICM(maxIter_ICMSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsLogLik(epsLogLikSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type homo(homoSEXP);
    Rcpp::traits::input_parameter< const bool& >::type diagSigmak(diagSigmakSEXP);
    rcpp_result_gen = Rcpp::wrap(icmem_heterCpp(X, Adj, y_int, Mu_int, W_int, Sigma_int, Lam_vec_int, alpha, beta_int, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, homo, diagSigmak));
    return rcpp_result_gen;
END_RCPP
}
// EMmPCpp_heter
Rcpp:: List EMmPCpp_heter(const arma::mat& X, const arma::vec& Pi_int, const arma::mat& Mu_int, const arma::mat& W_int, const arma::cube& Sigma_int, const arma::vec& Lam_vec_int, const int& maxIter, const double& epsLogLik, const bool& verbose, const bool& homo, const bool& diagSigmak);
RcppExport SEXP _DR_SC_EMmPCpp_heter(SEXP XSEXP, SEXP Pi_intSEXP, SEXP Mu_intSEXP, SEXP W_intSEXP, SEXP Sigma_intSEXP, SEXP Lam_vec_intSEXP, SEXP maxIterSEXP, SEXP epsLogLikSEXP, SEXP verboseSEXP, SEXP homoSEXP, SEXP diagSigmakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi_int(Pi_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu_int(Mu_intSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W_int(W_intSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Sigma_int(Sigma_intSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lam_vec_int(Lam_vec_intSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsLogLik(epsLogLikSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type homo(homoSEXP);
    Rcpp::traits::input_parameter< const bool& >::type diagSigmak(diagSigmakSEXP);
    rcpp_result_gen = Rcpp::wrap(EMmPCpp_heter(X, Pi_int, Mu_int, W_int, Sigma_int, Lam_vec_int, maxIter, epsLogLik, verbose, homo, diagSigmak));
    return rcpp_result_gen;
END_RCPP
}
// calculateWeight
arma::vec calculateWeight(const arma::mat& X, const int& nPCs);
RcppExport SEXP _DR_SC_calculateWeight(SEXP XSEXP, SEXP nPCsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int& >::type nPCs(nPCsSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateWeight(X, nPCs));
    return rcpp_result_gen;
END_RCPP
}
// wpcaCpp
Rcpp::List wpcaCpp(const arma::mat& X, const int& nPCs, const bool& weighted);
RcppExport SEXP _DR_SC_wpcaCpp(SEXP XSEXP, SEXP nPCsSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int& >::type nPCs(nPCsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(wpcaCpp(X, nPCs, weighted));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DR_SC_getPairDist", (DL_FUNC) &_DR_SC_getPairDist, 1},
    {"_DR_SC_calYenergy2D_sp", (DL_FUNC) &_DR_SC_calYenergy2D_sp, 5},
    {"_DR_SC_obj_beta", (DL_FUNC) &_DR_SC_obj_beta, 6},
    {"_DR_SC_runICM_sp", (DL_FUNC) &_DR_SC_runICM_sp, 11},
    {"_DR_SC_icmem_heterCpp", (DL_FUNC) &_DR_SC_icmem_heterCpp, 16},
    {"_DR_SC_EMmPCpp_heter", (DL_FUNC) &_DR_SC_EMmPCpp_heter, 11},
    {"_DR_SC_calculateWeight", (DL_FUNC) &_DR_SC_calculateWeight, 2},
    {"_DR_SC_wpcaCpp", (DL_FUNC) &_DR_SC_wpcaCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_DR_SC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
