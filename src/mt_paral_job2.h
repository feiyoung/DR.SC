#ifndef mt_paral_job2_hpp
#define mt_paral_job2_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "simulDRcluster.h"

using namespace std;
using namespace arma;
using namespace Rcpp;


class par_DRSC_nonspa{
public:
    
	mat X;
    field<mat> Mu0;
    field<cube> Sigma0;
    arma::mat W0;
    vec Lam_vec0;
    field<vec> Pi0;
    int maxIter_ICM;
    int maxIter;
    double epsLogLik;
    bool verbose;
    bool homo;
    bool diagSigmak;
    int maxK, minK;
    int g;    
    int current_idx;
    mat out_param;
    struct Objdrsc_nonspa output[50];

    

	par_DRSC_nonspa(const mat& X, const field<mat>& Mu0, field<cube> Sigma0, const arma::mat& W0,
                      const vec& Lam_vec0, const field<vec>& Pi0,
                      const int& maxIter, const double & epsLogLik, const bool& verbose,
                      const bool& homo, const bool& diagSigmak, 
                      const int maxK, const int minK, mat& out_param){

		this->X = X;
		this->Mu0 = Mu0;
		this->Sigma0 = Sigma0;
		this->W0 = W0;
		this->Lam_vec0 = Lam_vec0;
      this->Pi0 = Pi0;
        // this->maxIter_ICM = maxIter_ICM;
        this->maxIter = maxIter;
        this->epsLogLik = epsLogLik;
        this->verbose = verbose;
        this->homo = homo;
        this->diagSigmak = diagSigmak;
        this->maxK = maxK;
        this->minK = minK;
        this->current_idx = minK;
        this->out_param = out_param;
        
	}

	void loop_by_K_drsc(int g);
	void update_by_thread_drsc(int thread_id);
	int  next_drsc();

};


#endif 

