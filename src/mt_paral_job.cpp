#include "mt_paral_job.h"
#include "simulDRcluster.h"



using namespace std;
using namespace arma;
using namespace Rcpp;







  


void par_DRSC::loop_by_K_drsc(int g){
    

    
  mat Mu_int = Mu0(g);
  vec Lam_vec_int = Lam_vec0; 
  mat W_int = W0;
  vec alpha_int = alpha0(g);
  cube Sigma_int = Sigma0(g);
  ivec y_int = y.col(g);
  double beta_int = beta0(g);
  int K = Mu_int.n_rows;
    
    /*
    cout << "accuX-" << g << " " << accu(X) << endl;
    cout << "accuAdj-" << g << " " << accu(Adj) << endl;
    cout << "accuy_int-" << g << " " << accu(y_int) << endl;
    cout << "accuMu_int-" << g << " " << accu(Mu_int) << endl;
    cout << "accuSigma_int-" << g << " " << accu(Sigma_int) << endl;
    cout << "accuW_int-" << g << " " << accu(W0) << endl;
    cout << W_int(span(0,4),span(0,4)) << endl;
    cout << "accuLam_vec_int-" << g << " " << accu(Lam_vec_int) << endl;
    cout << "accualpha_int-" << g << " " << accu(alpha_int) << endl;
    */
    
    output[g] = drsc(X, Adj, y_int, Mu_int, Sigma_int,  W_int,
                        Lam_vec_int, alpha_int,  beta_int,  beta_grid,
                        maxIter_ICM,  maxIter, epsLogLik, verbose,
                        homo, diagSigmak);
    
    int p = X.n_cols;
    int q = W_int.n_cols;
    int n = X.n_rows;
    
    // caluculate degree of freedom

    int pen_const = 1;
    
    // calucate AIC and BIC
    double aic = -2.0* output[g].loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* 2* log(log(p+n))*pen_const; // adjusted  bic and aic for high dimension
    double bic =  -2.0* output[g].loglik + (1+p*(q+1) + K*(q+q*(q+1)/2.0))* log(n)* log(log(p+n))*pen_const; 
    

    // output the matrix of paramters AIC BIC n p and degree of freedom
    out_param(g, 0) = aic;
    out_param(g, 1) = bic;
    out_param(g, 2) = n;
    out_param(g, 3) = p;
    out_param(g, 4) = q;
    out_param(g, 5) = 1+p*(q+1) + K*(q+q*(q+1)/2.0);


}

std::mutex _mtx22;
int par_DRSC::next_drsc(){
	std::lock_guard<std::mutex> lockGuard(_mtx22);
	if (current_idx >= maxK - minK + 1){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void par_DRSC::update_by_thread_drsc(int thread_id){
	while (true){
		int idx = next_drsc();
		if (idx == -1){
			break;
		}
		loop_by_K_drsc(idx);
	}
}





