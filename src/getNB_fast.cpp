#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;
using namespace std;


//' @title
//' getneighborhood_fast
//' @description
//' an efficient function to find the neighborhood based on the matrix of position and a pre-defined cutoff
//'
//' @param x is a n-by-2 matrix of position.
//' @param radius is a threashold of Euclidean distance to decide whether a spot is an neighborhood of another spot. For example, if the Euclidean distance between spot A and B is less than cutoff, then A is taken as the neighbourhood of B. 
//' @return A sparse matrix containing the neighbourhood
//'
//' @export
// [[Rcpp::export]]
arma::sp_umat getneighborhood_fast(const arma::mat x, double radius)	{
  int N = x.n_rows;
  arma::sp_umat D(N, N);
  double dis;
  uvec idx, idx2;
  for (int j = 0; j < N-1; ++j)
  {    
    idx = find(abs(x(j,0) - x.col(0))<radius); 
    idx2 = find(idx>j);
    int p = idx2.n_elem;
    for (int i = 0; i < p; ++i)
    {
      dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
      if (dis < radius){
        D(idx(idx2(i)),j) = 1;
        D(j,idx(idx2(i))) = 1;
      }
    }
  }
  return D;
}


//' Calculate column-wise or row-wise mean
//' @param sp_data A sparse matrix
//' @param rowMeans A boolean value, whether to calculate row-wise mean
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_means_Rcpp(arma::sp_mat sp_data, bool rowMeans = false) {
  
  arma::sp_mat norm_col_sums;
  arma::mat tmp_mat;
  
  if (rowMeans) {
    norm_col_sums = arma::mean(sp_data, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}
  else {
    
    norm_col_sums = arma::mean(sp_data, 0);
    
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0).t());
  }
  
  return tmp_mat;
}


//' Calculate column-wise or row-wise sum
//' @param sp_data A sparse matrix
//' @param rowSums A boolean value, whether to calculate row-wise sum
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_sums_Rcpp(arma::sp_mat sp_data, bool rowSums = false) {
  arma::mat tmp_mat;
  
  arma::sp_mat norm_col_sums;
  
  if (rowSums) {
    norm_col_sums = arma::sum(sp_data, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}
  else {
    norm_col_sums = arma::sum(sp_data, 0);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0).t());
  }
  
  return tmp_mat;
}

