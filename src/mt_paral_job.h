#ifndef mt_paral_job_hpp
#define mt_paral_job_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "simulDRcluster.h"

using namespace std;
using namespace arma;
using namespace Rcpp;


class par_DRSC{
public:
    
    int current_idx = 0;
    
	mat X;
    sp_mat Adj;
    imat y;
    field<mat> Mu0;
    field<cube> Sigma0;
    mat W0;
    vec Lam_vec0;
    field<vec> alpha0;
    vec beta0;
    vec beta_grid;
    int maxIter_ICM;
    int maxIter;
    double epsLogLik;
    bool verbose;
    bool homo;
    bool diagSigmak;
    int maxK, minK;
    int g;    

    mat out_param;
    struct Objdrsc output[50];

    

	par_DRSC(const mat& X, const sp_mat& Adj, const imat& y,
                      const field<mat>& Mu0, field<cube> Sigma0, const arma::mat& W0,
                      const vec& Lam_vec0, const field<vec>& alpha0, vec beta0, const arma::vec& beta_grid,
                      const int& maxIter_ICM, const int& maxIter, const double & epsLogLik, const bool& verbose,
                      const bool& homo, const bool& diagSigmak, 
                      const int maxK, const int minK, mat& out_param){

		this->X = X;
		this->Adj = Adj;
		this->y = y;
		this->Mu0 = Mu0;
		this->Sigma0 = Sigma0;
		this->W0 = W0;
		this->Lam_vec0 = Lam_vec0;
        this->alpha0 = alpha0;
        this->beta0 = beta0;
        this->beta_grid = beta_grid;
        this->maxIter_ICM = maxIter_ICM;
        this->maxIter = maxIter;
        this->epsLogLik = epsLogLik;
        this->verbose = verbose;
        this->homo = homo;
        this->diagSigmak = diagSigmak;
        this->maxK = maxK;
        this->minK = minK;
        this->out_param = out_param;
        
	}

	void loop_by_K_drsc(int g);
	void update_by_thread_drsc(int thread_id);
	int  next_drsc();

};


#endif 

