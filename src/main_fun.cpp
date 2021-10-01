// // This script considers the case with the homogenous error and without spatial information.
// 
// #include "RcppArmadillo.h"
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// #define INT_MIN (-INT_MAX - 1)
// 
// using namespace Rcpp;
// using namespace arma;
// using namespace std;
// 
// 
// /*
//  * Auxiliary
//  */
// //' @keywords internal
// //' @noRd
// //' 
// //' Data X has been centered.
// 
// //' Evaluate the determinant of the covariance matrix and energy function values on each sample
// //' point, by SVD decomposition to achieve efficeint computation.
// //' @param X: a n x p  matrix, such as gene expression logcount matrix
// //' @param W0: a p x q matrix, loading matrix
// //' @param Ck: a q x q matrix,
// //' @param Mu0: a q x K matrix, gaussian mixture mean component
// //' @param dSk: a real number, the returned determinant of covariance matrix
// //' @param mSk: a colunmn vector, the returned energy funciton values. 
// //' 
// //' @return No return
// //' 
// 
// void multi_det_Sk_homo(const arma::mat& X, const arma::mat& W0, const arma::mat& Ck, 
//                      const arma::mat& Mu0, 
//                       int k, double& dSk, arma::vec& mSk){
//   //int p = X.n_cols;
//   int n = X.n_rows;
//   mSk = zeros(n);
//   dSk = 0;
//   
//   mat WC12,  tmp2;
//   vec tmp1, s, tmp3;
//   mat U, V, X_tk;
//   
//   svd(U, s, V, Ck.i());
//   
//   WC12 = W0 * (U * diagmat(sqrt(s)));
//   vec d = svd(WC12);
//   dSk = arma::as_scalar(prod(1- d % d));
//   X_tk = X - repmat(Mu0.row(k)* W0.t(), n, 1);
//   tmp1 = sum(X_tk % X_tk, 1);
//   tmp2 = X_tk * WC12;
//   tmp3 = sum(tmp2 % tmp2, 1);
//   mSk = tmp1 - tmp3;
// }
// vec decomp_homo(const mat& Cki, const mat& W0){
//   vec s, tmp1;
//   mat U, V, WC12;
//   svd(U, s, V, Cki);
//   WC12 = W0 * (U * diagmat(sqrt(s)));
//   tmp1 = sum(WC12 % WC12, 1);
//   return tmp1;
// }  
// 
// mat update_homo_W0(const mat& X, const mat& R, const cube& Ez, const cube& Ci, const double& sigma20, const vec& N){
//   int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
//   mat tmpMat(n,q, fill::zeros), Ezzt(q,q, fill::zeros), tmpMat2;
//   // vec N = sum(R.t(), 1);
//   for(k=0; k<K; ++k){
//     tmpMat2= repmat(R.col(k), 1, q) % Ez.slice(k);
//     tmpMat += tmpMat2;
//     Ezzt+= tmpMat2.t() * Ez.slice(k) + N(k) * sigma20 * Ci.slice(k);
//   }
//   return X.t() * tmpMat * Ezzt.i();
// }
// 
// // update Sigma0
// cube update_homo_Sigma0(const mat& R, const cube& Ez, const cube& Ci, const mat& Mu, const double& sigma20, const vec&  N){
//   int k, K= R.n_cols, q= Mu.n_cols, n= R.n_rows;
//   cube Sigma0(q,q,K);
//   for(k = 0; k<K; ++k){
//     Sigma0.slice(k) = trans(Ez.slice(k) - repmat(Mu.row(k), n, 1)) * diagmat(R.col(k)) * (Ez.slice(k) - repmat(Mu.row(k), n, 1)) + N(k)* sigma20 * Ci.slice(k);
//     Sigma0.slice(k) = Sigma0.slice(k) / N(k);
//   }
//   return(Sigma0);
// }
// 
// // update sigma20
// double update_sigma20(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci, const double& sigma20){
//   int k, K= R.n_cols;
//   mat tmpXk;
//   double term1=0, term2=0;
//   for(k=0; k <K; ++k){
//     tmpXk = (X - Ez.slice(k) * W.t() );
//     term1 += accu(sum(tmpXk % tmpXk, 1) % R.col(k));
//     term2  += accu(decomp_homo(Ci.slice(k), W))* accu(R.col(k));
//   }
//   return((term1+ sigma20* term2)/ X.n_elem);
// }
// 
// // Assume E(X) = 0.
// // [[Rcpp::export]]
// Rcpp:: List EMmPCpp(const arma::mat& X, const arma::vec& Pi_int, const arma::mat& Mu_int, const arma::mat&
//                       W_int, const arma::cube& Sigma_int, const double& sigma2_int,
//                       const int& maxIter, const double& epsLogLik, const int& output){
//   // basic info
//   int p = X.n_cols;
//   int n = X.n_rows;
//   int K = Pi_int.size();
//   int q = Mu_int.n_cols;
//   // obtain intial vlaues
//   vec Pi0(Pi_int);
//   mat Mu0(Mu_int), W0(W_int);
//   cube Sigma0(Sigma_int);
//   double sigma20(sigma2_int);
//   
//   // temperary objects
//   List resList;
//   // If p is sufficient large, loglik can not be computed.
//   // But this can be solved by some programming tricks.
//   vec loglik(maxIter);
//   loglik(0) = INT_MIN;
//   vec maxA1(n,fill::zeros);
//   vec loglik_more_vec =maxA1;
//   
//   // Define the variables that will be used in algorithm
//   // variables usded in updating Pi0
//   double  dSk;
//   cube Cki_ara(q, q, K, fill::zeros);
//   mat A1(n,K, fill::zeros), Ck(q,q, fill::zeros);
//   mat R(A1); 
//   vec mSk(n);
//   int k, iter;
//   
//   // variables usded in updating Mu0
//   cube Ez(n,q, K,  fill::zeros);
//   
//   
//   // variables usded in calculating Q value.
//   // arma::mat L(n,K, fill::zeros);
//   // arma::vec vecL;
//   
//   // begin algorithm
//   for(iter = 1; iter < maxIter; iter++){
//     // cache some objects
//     
//     // compute loglikelihood
//     for(k=0; k<K; k++){
//     
//       Ck = W0.t() * W0 + sigma20* inv(Sigma0.slice(k));
//       Cki_ara.slice(k) = Ck.i();
//       multi_det_Sk_homo(X, W0, Ck, Mu0, 
//                       k, dSk, mSk);
//       
//       A1.col(k) = log(Pi0(k)) + 0.5*log(dSk) -1/(2* sigma20) * mSk;
//       Ez.slice(k) = (X*W0 + sigma20* repmat(Mu0.row(k)*inv(Sigma0.slice(k)), n, 1)) * Ck.i();
//      
//     } 
//     
//     // loglik(iter) =  sum(log(sum(exp(A1),1))) - n* p/2.0 * log(2* M_PI* sigma20); 
//     maxA1 = max(A1, 1);
//     A1 = (A1 - repmat(maxA1, 1, K));
//     
//     loglik_more_vec = sum(exp(A1),1);
//     loglik(iter) = sum(log(loglik_more_vec) + maxA1) - n* p /2.0 * log(2* M_PI* sigma20); 
//     R = exp(A1) / repmat(loglik_more_vec, 1, K);
//     // A.replace(datum::nan, 0);
//     
//     // compute Q(theta^t; theta^t)
//    
//     // update Pi0
//     vec N = arma::sum(R.t(),1);
//     Pi0 = N / accu(N);
//     // cout<<"Pi0="<<trans(Pi0)<<endl;
//     
//     // update Mu0
//     for(k=0; k<K; ++k){
//       Mu0.row(k) = trans(R.col(k)) * Ez.slice(k) / N(k);
//     }
//     
//     // update Sigma20
//     Sigma0 = update_homo_Sigma0(R, Ez, Cki_ara, Mu0, sigma20, N);
//     
//     // update W0
//     W0 = update_homo_W0(X,  R, Ez, Cki_ara, sigma20, N);
//    
//     // update sigma2
//     sigma20 = update_sigma20(R, X, W0, Ez, Cki_ara, sigma20);
//     
//     
//     // calculate loglikelihood
//     if(loglik(iter)  - loglik(iter-1)   < -1e-7){
//       // error_flag = 1;
//       perror("The likelihood failed to increase!");
//     }
//    
//     // // calculate Q function value
//     
//        
//     // output algorithm info.
//     if(output){
//       cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
//       // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;  
//     }
//     if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
//      
//   }
//   mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
//   for(k=0; k<K; k++){
//     
//     Ezz +=  Ez.slice(k) % repmat(R.col(k), 1, q);
//   }
//   
//   // output return value
//   resList["cluster"] = index_max(R,1) + 1;
//   resList["Ez"] =Ezz;
//   resList["Pi"] = Pi0;
//   resList["Mu"] = Mu0;
//   resList["Sigma2"] = Sigma0;
//   resList["W"] = W0;
//   resList["sigma2"] = sigma20;
//   
//   resList["loglik"] = loglik(iter-1);
//   resList["loglik_seq"] = loglik.subvec(0, iter-1);
//  
//   resList["dLogLik"] = loglik(iter-1)  - loglik(iter-2);
//   // resList["aic"] = -2.0* loglik(iter-1) + (p*(q+1) + K*(1+q+q*(q-1)/2.0))* 2;
//   // resList["bic"] = -2.0* loglik(iter-1) + (p*(q+1) + K*(1+q+q*(q-1)/2.0))* log(n) ;
//   return(resList);
// } 