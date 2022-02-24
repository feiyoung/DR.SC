// This script considers the case with the hetergenous error with/without spatial information.
// lightly edited version of question, working fine for me
// #define ARMA_64BIT_WORD 1

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <time.h>
//#include <iostream>
#include "utilSimulDRcluster.h"
#include "simulDRcluster.h"
#include "mt_paral_job.h"
#include "mt_paral_job2.h"

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;






/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 
//' Data X has been centered.

//[[Rcpp::export]]
arma::mat getPairDist(const arma::mat x)	{
  int N = x.n_rows;
  arma::mat D(N, N);
  for (int j = 0; j < N; ++j)
  {
    for (int i = j; i < N; ++i)
    {
      D(i, j)	= norm(x.row(i) - x.row(j), 2);
      D(j, i) = D(i, j);
    }
  }
  
  return D;
}



sp_mat get_spNbs(ivec y, const sp_mat& Adj) {   // ivec是索引型向量
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i. 
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin(); //构造一个sp_mat的常数迭代器,常数迭代器只可读，不可写，实现对矩阵按照列对每个非零元素进行访问。
  //arma::sp_mat::const_iterator end   = Adj.end();
  
  // Calculate number of nonzero points
  //int n = std::distance(start, end);
  int n = Adj.n_nonzero; // 计算Adj的所有非零元的个数
  //cout << "n=" << n << endl;
  //cout << "n=" << Adj.n_nonzero << endl;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
  
  
  arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
  for(int i = 0; i < n; ++i)
  {
    //temp(0) = it.row();
    //temp(1) = it.col();
    spNbs(it.row(), it.col()) = y(it.col()); // it只自加非零元个数次，得到每个i对应的邻居的状态
    ++it; // increment
  }
  
  return spNbs.t(); // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
}

// [[Rcpp::export]]
arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  // Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
  int n = y.n_rows;
  arma::sp_mat spNbs_t = get_spNbs(y, Adj); // transform spNbs to iterate by column.
  arma::mat Uy(n, K);
  double n_sameS;
  int i, k, nn;
  for (k = 0; k < K; k++)
  {
    for (i = 0; i < n; i++)
    {
      arma::sp_mat col(spNbs_t.col(i)); // the class label of neighbors of i-th sample.
      n_sameS = 0;
      
      nn = col.n_nonzero; // the number of neighbors of i-th sample
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += ((*j) == (k+1));
        
      }
      Uy(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
      
      
    }
  }
  
  arma::mat C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
  Uy = -log(C_mat); // normalized Uy, this is the energy of y.
  return Uy;
  
}
// Suppose we know true y, can we estimate beta correctly by maximizing pseudo loglikelihood  
// //[[Rcpp::export]] 
// double obj_beta2(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
//   
//   mat Uy = calYenergy2D_sp(y, Adj, K, alpha, beta);
//   arma::mat C_mat = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE to get the likelihood
//   double loglike = 0;
//   int n = y.n_elem;
//   int i = 0;
//   for(; i < n; ++i)
//   {
//     
//     loglike += log(C_mat(i, y(i)-1));
//   }
//   return loglike;
// }

//[[Rcpp::export]]  
double obj_beta(const arma::ivec& y, const arma::mat& R, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  
  mat Uy = calYenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
  //arma::mat C_mat = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE to get the likelihood
  //return accu(R % log(C_mat)); 
  return -accu(R % Uy);
}

void runICM_sp (const arma::mat& X,  arma::ivec& y, const arma::mat& W0, const arma::vec& Lam_vec0, const arma::mat& Mu0,
                      const arma::cube& Sigma0,const arma::sp_mat& Adj, const arma::vec& alpha, const arma::vec& beta_grid,
                      double& beta, int maxIter_ICM, mat& R, cube& Ez, cube& Cki_ara, double& loglik)	{
  // Target: estimate Y, evaluate R, Ez, Ck inverse, and update beta by using grid search.
  
  // basic info.
  int n = X.n_rows, K = Mu0.n_rows, q= Mu0.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  // cube Cki_ara(q, q, K, fill::zeros), Ez(n,q,K, fill::zeros);
  double  logdSk;
  vec mSk(n);
  mat WtLW = W0.t() * (repmat(1.0/ Lam_vec0, 1, q) %  W0); //O(p^2 q) cache object 
  mat XLW = X * (repmat(1.0/ Lam_vec0, 1, q) % W0);
  // evaluate energy of x, Ux
  arma::mat Ux(n, K), Ck;

  for (k = 0; k < K; k++)	{

    Ck = WtLW +  inv_sympd(Sigma0.slice(k));

    Cki_ara.slice(k) = inv_sympd(Ck);

    multi_det_SkCpp2(X, Lam_vec0,W0, Ck, Mu0.row(k), Sigma0.slice(k), // Use SVD to speed up.
                    logdSk, mSk);

    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) = -0.5*logdSk  + 0.5 * mSk; // calculate energy by column.

    Ez.slice(k) = (XLW + repmat(Mu0.row(k)* inv_sympd(Sigma0.slice(k)), n, 1)) * Ck.i();
  }
 
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------	
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, alpha, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    
    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  loglik = sum(log(loglik_more_vec) + maxA1); //  - n* p /2.0 * log(2* M_PI); 
  R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, alpha, beta_grid(k));
  }
  beta = beta_grid(index_max(objBetaVec));
  
  
  // List output = List::create(
  //   Rcpp::Named("y") = y,
  //   Rcpp::Named("R") = R,
  //   Rcpp::Named("Ez") = Ez,
  //   Rcpp::Named("Cki_ara") = Cki_ara,
  //   Rcpp::Named("loglik") = loglik,
  //   Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  // return output; 
  
}  



Objdrsc drsc(const arma::mat& X, const arma::sp_mat& Adj, arma::ivec& y, mat& Mu0,
                    cube& Sigma0, arma::mat& W0, arma::vec& Lam_vec0,
                    arma::vec& alpha, double& beta0, arma::vec& beta_grid,int& maxIter_ICM, 
                    int& maxIter, double& epsLogLik, bool& verbose, bool& homo,
                    bool& diagSigmak){
    
    int K = Mu0.n_rows;
    int q = Mu0.n_cols;
    int n = X.n_rows;
    
    
  // If p is sufficient large, loglik can not be computed.
  // But this can be solved by some programming tricks.
  vec loglik(maxIter), N(K, fill::zeros);
  loglik(0) = INT_MIN; // -INFINITY
  
  // Define the variables that will be used in algorithm
  // variables usded in eavluating R.
  cube Cki_ara(q, q, K, fill::zeros), Ez(n,q,K, fill::zeros);
  mat  R(n,K, fill::zeros);
  int k, iter;
  double loglikVal =0.0;
    
    // begin algorithm
  for(iter = 1; iter < maxIter; iter++){
    
    
    // estimate Y, update beta, and cache some objects
    // cout<<"Start computing ICM algorithm"<<endl;
    runICM_sp(X, y, W0, Lam_vec0, Mu0, Sigma0, Adj,alpha, beta_grid, beta0, maxIter_ICM,
                        R, Ez, Cki_ara, loglikVal);
    
    loglik(iter) = loglikVal;
    
    
    // cout<< R << Ez << Cki_ara << beta0 <<endl;
    // compute Q(theta^t; theta^t)
    
    // double Q0 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    
    
    // double Q1 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Pi="<< Q1-Q0<<endl;
    // update Mu0
    N = arma::sum(R.t(), 1);
    for(k=0; k<K; ++k){
        
      Mu0.row(k) = trans(R.col(k)) *  Ez.slice(k) / N(k);
    }
    
    // double Q2 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    
    
    // update Sigma0
    Sigma0 = update_Sigma0(R, Ez, Cki_ara, Mu0, N, diagSigmak);
    
    // double Q3 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Sigma0="<< Q3-Q2<<endl;
    
    // update W
    W0 = update_W0(X, R, Ez, Cki_ara, N);
    
    // double Q4 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    
    // update  Lambda
    //cout<<"start updating Lambda_vec"<<endl;
    Lam_vec0 = update_Lam(R, X, W0, Ez, Cki_ara, homo);
    //cout<<"Finish updating Lambda_vec"<<endl;
    // Lam_vec0 = update_Lam2(R, X, W0, Ez, Cki_ara);
    // double Q5 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Lambda="<< Q5-Q4<<endl;
    
    // calculate loglikelihood
      /*
    if(loglik(iter)  - loglik(iter-1)   < -1e-7){
      // perror("The likelihood failed to increase!");
      break;
    }
    */
    
    // output algorithm info.
    if(verbose){
      // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
      // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
      Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
              iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
    }
    if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
    // if(abs(Q  - tmp_Q) < epsLogLik) break;
    
  }
  
  // output return value
  
  
  // mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  // cube Ez = ICM_fit["Ez"]; // Ez is also the problem.
  // Evalute the posterior expectation of z.
  
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    Ezz +=  Ez.slice(k) % repmat(R.col(k), 1, q);
  }
  
    
    Objdrsc output;
    
    output.y = y;
    output.Ezz = Ezz;
    output.beta0 = beta0;
    output.Mu0 = Mu0;
    output.Sigma0 = Sigma0;
    output.W0 = W0;
    output.Lam_vec0 = Lam_vec0;
    output.loglik = loglik(iter-1);
    output.loglik_seq = loglik.subvec(0, iter-1);
    
    return(output);
 }







// Assume E(X) = 0, consider the spatial information.
// [[Rcpp::export]]
Rcpp:: List icmem_heterCpp(const arma::mat& X,const arma::sp_mat& Adj, const arma::imat& y_int, Rcpp::List& Mu_intList, const arma::mat&
  W_int, Rcpp::List& Sigma_intList, arma::vec& Lam_vec_int,
  Rcpp::List& alphaList, const arma::vec& beta_int, const arma::vec& beta_grid,const int& maxIter_ICM, 
  const int& maxIter, const double& epsLogLik, const int& verbose,const bool& homo = false, const bool& diagSigmak = false, const int maxK = 10, const int minK = 2, const int coreNum = 1){
  

  // basic information
  int n = X.n_rows, q= W_int.n_cols;
  
    
  // Initialize the  iterative parameters    
    
    int lengthK = maxK - minK + 1;
    field<mat> Mu0(lengthK);
    field<cube> Sigma0(lengthK);
    imat y0 = y_int;
    mat W0 = W_int;
    vec beta0 = beta_int;
    field<vec> alpha0(lengthK);
    vec Lam_vec0 = Lam_vec_int;
    
    for (int i = 0; i < lengthK; i++){
        mat Mu0_tmp = Mu_intList[i];
        Mu0(i) = Mu0_tmp;
        cube Sigma0_tmp = Sigma_intList[i];
        Sigma0(i) = Sigma0_tmp;
        vec alpha0_tmp = alphaList[i];
        alpha0(i) = alpha0_tmp;     
    }
  
    
    mat out_param = zeros<mat>(lengthK, 6);
    
    //set parallel structure object
    par_DRSC parObj(X, Adj, y0,  Mu0, Sigma0, W0, Lam_vec0, 
                        alpha0, beta0, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose,
                        homo, diagSigmak, maxK, minK, out_param);

	const int n_thread = coreNum;
	std::vector<std::thread> threads(n_thread);
    
	for (int i_thread = 0; i_thread < n_thread; i_thread++){
		threads[i_thread] = std::thread(&par_DRSC::update_by_thread_drsc, &parObj, i_thread);
	}
	for (int i = 0; i < n_thread; i++){
		threads[i].join();
	}

    List ObjdrscRcpp(maxK-minK+1);
    
    out_param = parObj.out_param;
    
    for (int k = 0; k<lengthK; k++){
        ObjdrscRcpp[k] = List::create(
        Rcpp::Named("cluster") = parObj.output[k].y,
        Rcpp::Named("hZ") = parObj.output[k].Ezz,
        Rcpp::Named("beta") = parObj.output[k].beta0,
        Rcpp::Named("Mu") = parObj.output[k].Mu0,
        Rcpp::Named("Sigma") = parObj.output[k].Sigma0,
        Rcpp::Named("W") = parObj.output[k].W0,
        Rcpp::Named("Lam_vec") = parObj.output[k].Lam_vec0,
        Rcpp::Named("loglik") = parObj.output[k].loglik,
        Rcpp::Named("loglik_seq") = parObj.output[k].loglik_seq);
    }
    
  List resList = List::create(
    Rcpp::Named("Objdrsc") = ObjdrscRcpp, 
    Rcpp::Named("out_param") = out_param);
    
  return(resList);
}









// without considering the spatial information
void runEstep(const arma::mat& X,  const arma::mat& W0, const arma::vec& Lam_vec0, const arma::mat& Mu0,
                const arma::cube& Sigma0, const arma::vec& Pi0,  mat& R, cube& Ez, cube& Cki_ara, double& loglik)	{
  // Target: estimate Y, evaluate R, Ez, Ck inverse, and update beta by using grid search.
  
  // basic info.
  int n = X.n_rows, p=X.n_cols, K = Mu0.n_rows, q= Mu0.n_cols;
  int k;

  // two cached objects used for parameters update.
  // cube Cki_ara(q, q, K, fill::zeros), Ez(n,q,K, fill::zeros);
  double  logdSk;
  vec mSk(n);
  vec maxA1(n,fill::zeros);
  vec loglik_more_vec =maxA1;

  mat WtLW = W0.t() * (repmat(1.0/ Lam_vec0, 1, q) %  W0); //O(p^2 q) cache object 
 
  mat XLW = X * (repmat(1.0/ Lam_vec0, 1, q) % W0);

  // evaluate energy of x, Ux
  arma::mat A1(n, K), Ck;

  for (k = 0; k < K; k++)	{
    Ck = WtLW +  inv_sympd(Sigma0.slice(k));
 
    Cki_ara.slice(k) = Ck.i();

    multi_det_SkCpp2(X, Lam_vec0,W0, Ck, Mu0.row(k), Sigma0.slice(k), // Use SVD to speed up.
                     logdSk, mSk);
  
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    A1.col(k) = log(Pi0(k)) + 0.5*logdSk  -0.5 * mSk;

    Ez.slice(k) = (XLW + repmat(Mu0.row(k)* inv_sympd(Sigma0.slice(k)), n, 1)) * Ck.i();
  }
  
  // calculate R and pseudo observed loglikelihood
  maxA1 = max(A1, 1);
  A1 = (A1 - repmat(maxA1, 1, K));
  loglik_more_vec = sum(exp(A1),1);
  loglik = sum(log(loglik_more_vec) + maxA1) - n* p /2.0 * log(2* M_PI); 
  R = exp(A1) / repmat(loglik_more_vec, 1, K);
 
}  








Objdrsc_nonspa drsc_nonspa(const arma::mat& X, const arma::vec& Pi_int, const arma::mat& Mu_int, const arma::mat&
  W_int, const arma::cube& Sigma_int, const  arma::vec& Lam_vec_int,
  const int& maxIter, const double& epsLogLik, const bool& verbose, 
  const bool& homo, const bool& diagSigmak){
    
    // basic info
  int  n = X.n_rows, K = Mu_int.n_rows, q = Mu_int.n_cols; // p = X.n_cols,
  int iter, k;

  // Initialize the  iterative parameters
  mat Mu0(Mu_int), W0(W_int);
  vec Pi0(Pi_int), Lam_vec0(Lam_vec_int);
  cube Sigma0(Sigma_int);
  
  // If p is sufficient large, loglik can not be computed.
  // But this can be solved by some programming tricks.
  vec loglik(maxIter);
  loglik(0) = INT_MIN;

  // Define the variables that will be used in algorithm
  // vec mSk(n); double  logdSk; //mat A1(n,K, fill::zeros), Ck(q,q, fill::zeros);
  cube Cki_ara(q, q, K, fill::zeros), Ez(n,q,K, fill::zeros);
  mat  R(n, K, fill::zeros);
  double loglikVal =0.0;
  

  // begin EM algorithm
  for(iter = 1; iter < maxIter; iter++){

    
    // cache some objects
    
    // compute loglikelihood
    // for(k=0; k<K; k++){
    //   
    //   Ck = W0.t() * (repmat(1.0/ Lam_vec0, 1, q) %  W0) +  inv_sympd(Sigma0.slice(k));
    //   Cki_ara.slice(k) = Ck.i();
    //   multi_det_SkCpp2(X, Lam_vec0,W0, Ck, Mu0.row(k), Sigma0.slice(k),
    //                   logdSk, mSk);
    //   // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    //   A1.col(k) = log(Pi0(k)) + 0.5*logdSk  -0.5 * mSk;
    //   
    //   Ez.slice(k) = (X * (repmat(1.0/ Lam_vec0, 1, q) % W0) + 
    //     repmat(Mu0.row(k)*inv_sympd(Sigma0.slice(k)), n, 1)) * Ck.i();
    //   
    // } 
    // maxA1 = max(A1, 1);
    // A1 = (A1 - repmat(maxA1, 1, K));
    // loglik_more_vec = sum(exp(A1),1);
    // loglik(iter) = sum(log(loglik_more_vec) + maxA1) - n* p /2.0 * log(2* M_PI); 
    // R = exp(A1) / repmat(loglik_more_vec, 1, K);
    // cout<<"Finish R computing in CPP::"<<endl;
    runEstep(X,  W0, Lam_vec0,  Mu0, Sigma0, Pi0,  R, Ez, Cki_ara, loglikVal);
    loglik(iter) = loglikVal;

    // compute Q(theta^t; theta^t)
    // double Q0 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    

    // update Pi0
    vec N = arma::sum(R.t(), 1);
    Pi0 = N/ accu(N);
    // double Q1 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Pi="<< Q1-Q0<<endl;

    // update Mu0
    for(k=0; k<K; ++k){
      Mu0.row(k) = trans(R.col(k)) * Ez.slice(k) / N(k);
    }
    
    // double Q2 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    
    
    // update Sigma0

    Sigma0 = update_Sigma0(R, Ez, Cki_ara, Mu0, N, diagSigmak);
    
    // double Q3 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Sigma0="<< Q3-Q2<<endl;
    
    // update W

    W0 = update_W0(X, R, Ez, Cki_ara, N);
    
    // double Q4 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    
    
    // update  Lambda
    Lam_vec0 = update_Lam(R, X, W0, Ez, Cki_ara,homo);
    // Lam_vec0 = update_Lam2(R, X, W0, Ez, Cki_ara);
    // double Q5 = Q_fun(X, R, Ez,Cki_ara, W0, Mu0, Sigma0, Pi0, Lam_vec0);
    // cout<<"dQ_Lambda="<< Q5-Q4<<endl;
    
    // calculate loglikelihood
    if(loglik(iter)  - loglik(iter-1)   < -1e-7){
      //perror("The likelihood failed to increase!");
      break;
    }
    
    
    // output algorithm info.
    if(verbose){
      // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
      // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
      Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
              iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
    }
    if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
    // if(abs(Q  - tmp_Q) < epsLogLik) break;
    
  }
  
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    Ezz +=  Ez.slice(k) % repmat(R.col(k), 1, q);
  }
  
    
  Objdrsc_nonspa output;
    
    output.y = conv_to<ivec>::from(index_max(R,1) + 1);
    output.Ezz = Ezz;
    output.Pi0 = Pi0;
    output.Mu0 = Mu0;
    output.Sigma0 = Sigma0;
    output.W0 = W0;
    output.Lam_vec0 = Lam_vec0;
    output.loglik = loglik(iter-1);
    output.loglik_seq = loglik.subvec(0, iter-1);
    
    return(output);
    
}







// [[Rcpp::export]]
Rcpp:: List EMmPCpp_heter(const arma::mat& X,  Rcpp::List& Pi_int, Rcpp::List& Mu_int, const arma::mat&
  W_int, Rcpp::List& Sigma_int, arma::vec& Lam_vec_int,
  const int& maxIter, const double& epsLogLik, const bool& verbose, 
  const bool& homo = false, const bool& diagSigmak = false, const int& maxK = 10, const int& minK = 5, const int& coreNum = 1){
    
    
    
      // basic information
  int n = X.n_rows, q= W_int.n_cols;

    
  // Initialize the iterative parameters   
    int lengthK = maxK - minK + 1;
    field<vec> Pi0(lengthK);
    field<mat> Mu0(lengthK);
    field<cube> Sigma0(lengthK);
    mat W0 = W_int;
    vec Lam_vec0 = Lam_vec_int;

    for (int i = 0; i < lengthK; i++){
        vec Pi0_tmp = Pi_int(i);
        Pi0(i) = Pi0_tmp;
        mat Mu0_tmp = Mu_int(i);
        Mu0(i) = Mu0_tmp;
        cube Sigma0_tmp = Sigma_int(i);
        Sigma0(i) = Sigma0_tmp;
    }

    mat out_param = zeros<mat>(lengthK, 6);
    //set parallel structure object
    par_DRSC_nonspa parObj(X, Mu0, Sigma0, W0, Lam_vec0, 
                        Pi0, maxIter, epsLogLik, verbose,
                        homo, diagSigmak, maxK, minK, out_param);

	const int n_thread = coreNum;
	std::vector<std::thread> threads(n_thread);
    
	for (int i_thread = 0; i_thread < n_thread; i_thread++){
		threads[i_thread] = std::thread(&par_DRSC_nonspa::update_by_thread_drsc, &parObj, i_thread);
	}
	for (int i = 0; i < n_thread; i++){
		threads[i].join();
	}

    
    out_param = parObj.out_param;
        
    List Objdrsc_nonspaRcpp(lengthK);
    
    
    for (int i = 0; i < lengthK; i++){
              // output return value
        Objdrsc_nonspaRcpp[i] = List::create(
            Rcpp::Named("cluster") = parObj.output[i].y,
            Rcpp::Named("hZ") = parObj.output[i].Ezz,
            Rcpp::Named("Pi") = parObj.output[i].Pi0,
            Rcpp::Named("Mu") = parObj.output[i].Mu0,
            Rcpp::Named("Sigma2") = parObj.output[i].Sigma0,
            Rcpp::Named("W") = parObj.output[i].W0,
            Rcpp::Named("Lam_vec") = parObj.output[i].Lam_vec0,
            Rcpp::Named("loglik") = parObj.output[i].loglik,
            Rcpp::Named("loglik_seq") = parObj.output[i].loglik_seq);
    }
    
    
  List resList = List::create(
    Rcpp::Named("Objdrsc_nonspa") = Objdrsc_nonspaRcpp, 
    Rcpp::Named("out_param") = out_param);
    
  return(resList);
} 
