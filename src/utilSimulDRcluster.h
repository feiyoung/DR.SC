#ifndef utilSimulDRcluster_hpp

#define utilSimulDRcluster_hpp

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <string.h>
// #include <time.h>
// #include <iostream>


using namespace Rcpp;
using namespace arma;
using namespace std;
void multi_det_SkCpp2(const arma::mat& X, const arma::vec& Lam_vec0,
                      const arma::mat& W0, const arma::mat& Ck, 
                      const arma::rowvec Muk, const mat& Sigmak,
                      double& logdSk, arma::vec& mSk);

void multi_det_SkCpp(const arma::mat& X, const arma::vec& Lam_vec0, const arma::mat& W0, const arma::mat& Ck, 
                     const arma::rowvec Muk, 
                     double& logdSk, arma::vec& mSk);

vec decomp(const arma::mat& Cki, const arma::mat& W0);

mat update_W0(const mat& X, const mat& R, const cube& Ez, const cube& Ci,  const vec& N);

cube update_Sigma0(const mat& R, const cube& Ez, const cube& Ci, const mat& Mu, 
                   const vec&  N, const bool& diagSigmak);

vec update_Lam(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci,const bool& homo);

// vec update_Lam2(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci);

double Q_fun(const mat& X, const mat& R,  const cube& Ez, const arma::cube& Ci, 
             const mat& W0, const mat& Mu0, const cube& Sigma0, const vec& Pi0, 
             const vec& Lam_vec0);

#endif
