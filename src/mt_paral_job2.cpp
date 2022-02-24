#include "mt_paral_job2.h"
#include "simulDRcluster.h"



using namespace std;
using namespace arma;
using namespace Rcpp;







  


void par_DRSC_nonspa::loop_by_K_drsc(int g){
    
  mat Mu_int = Mu0(g-minK);
  vec Lam_vec_int = Lam_vec0;  
  cube Sigma_int = Sigma0(g-minK);
  vec Pi_int = Pi0(g-minK);
  int K = Mu_int.n_rows;
  mat W_int = W0;

    output[g-minK] = drsc_nonspa(X, Pi_int, Mu_int,  W_int, Sigma_int, 
                        Lam_vec_int, maxIter, epsLogLik, verbose,
                        homo, diagSigmak);
    
    int p = X.n_cols;
    int q = W_int.n_cols;
    int n = X.n_rows;
    
    // caluculate degree of freedom

    int pen_const = 1;
    
    // calucate AIC and BIC
    double aic = -2.0* output[g-minK].loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen_const; // adjusted  bic and aic for high dimension
    double bic =  -2.0* output[g-minK].loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n))*pen_const; 
    

    // output the matrix of paramters AIC BIC n p and degree of freedom
    out_param(g-minK, 0) = aic;
    out_param(g-minK, 1) = bic;
    out_param(g-minK, 2) = n;
    out_param(g-minK, 3) = p;
    out_param(g-minK, 4) = q;
    out_param(g-minK, 5) = 1+p*(q+1) + K*(q+q*(q+1)/2.0);


}

std::mutex _mtx23;
int par_DRSC_nonspa::next_drsc(){
	std::lock_guard<std::mutex> lockGuard(_mtx23);
	if (current_idx >= maxK + 1){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void par_DRSC_nonspa::update_by_thread_drsc(int thread_id){
	while (true){
		int idx = next_drsc();
		if (idx == -1){
			break;
		}
		loop_by_K_drsc(idx);
	}
}





