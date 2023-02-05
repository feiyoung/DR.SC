#ifndef idrsc2_hpp
#define idrsc2_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>


using namespace Rcpp;
using namespace arma;
using namespace std;


struct Objdrsc{
    ivec y;
    mat Ezz;
    vec beta0;
    mat Mu0;
    cube Sigma0;
    mat W0;
    vec Lam_vec0;
    double loglik;
    vec loglik_seq;
};

Objdrsc drsc(const arma::mat& X, const arma::sp_mat& Adj, const arma::ivec& y_int, const mat& Mu_int,
             const cube& Sigma_int, const arma::mat& W_int, const arma::vec& Lam_vec_int,
             const arma::vec& alpha, const double& beta_int, const arma::vec& beta_grid, const int& maxIter_ICM, 
             const int& maxIter, const double& epsLogLik, const bool& verbose, const bool& homo,
             const bool& diagSigmak);


struct Objdrsc_nonspa{
    ivec y;
    mat Ezz;
    vec Pi0;
    mat Mu0;
    cube Sigma0;
    mat W0;
    vec Lam_vec0;
    double loglik;
    vec loglik_seq;
};



Objdrsc_nonspa drsc_nonspa(const arma::mat& X, const arma::vec& Pi_int, const arma::mat& Mu_int, const arma::mat&
  W_int, const arma::cube& Sigma_int, const  arma::vec& Lam_vec_int,
  const int& maxIter, const double& epsLogLik, const bool& verbose, 
  const bool& homo, const bool& diagSigmak);







#endif /* CoMM_covar_pxem_ptr_hpp */
