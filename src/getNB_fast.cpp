#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


//' @title
//' getneighborhood_fast
//' @description
//' an efficient function to find the neighborhood based on the matrix of position and a pre-defined cutoff
//'
//' @param x is a n-by-2 matrix of position.
//' @param cutoff is a threashold of Euclidean distance to decide whether a spot is an neighborhood of another spot. For example, if the Euclidean distance between spot A and B is less than cutoff, then A is taken as the neighbourhood of B. 
//' @return A sparse matrix containing the neighbourhood
//'
//' @export
// [[Rcpp::export]]
arma::sp_umat getneighborhood_fast(const arma::mat x, double cutoff)	{
  int N = x.n_rows;
  arma::sp_umat D(N, N);
  double dis;
  uvec idx, idx2;
  for (int j = 0; j < N-1; ++j)
  {    
    idx = find(abs(x(j,0) - x.col(0))<cutoff); 
    idx2 = find(idx>j);
    int p = idx2.n_elem;
    for (int i = 0; i < p; ++i)
    {
      dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
      if (dis < cutoff){
        D(idx(idx2(i)),j) = 1;
        D(j,idx(idx2(i))) = 1;
      }
    }
  }
  return D;
}
