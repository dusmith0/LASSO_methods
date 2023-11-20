#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  if(a > lambda){
    return(a - lambda);
  }else if(a < -lambda){
    return(a + lambda);
  }else{
    return(0);
  }
}


// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  int n = Ytilde.size();
  double obj = sum(pow(Ytilde - Xtilde * beta,2))/(2 * n) + lambda * sum(arma::abs(beta));
  return obj;
}


// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here

  // Do not forget to include a method for assigning beta_start
  int m = beta_start.size();
  int n = Ytilde.size();
  arma::colvec beta = beta_start;
  arma::colvec beta_update = beta_start;
  double eps_check = 100;
  arma::colvec r = Ytilde - Xtilde * beta_start;
 
  int j = 1;
  while (abs(eps_check) > eps){
  
    for(int j = 1; j <= m; ++j){
      arma::colvec mat = (beta(j-1) + Xtilde.col(j-1).t() * r / n); 
      double a = mat[0];
      beta_update(j-1) = soft_c(a,lambda);
    
      r = r + Xtilde.col(j-1) * (beta(j-1) - beta_update(j-1));
    }
    
      eps_check = lasso_c(Xtilde,Ytilde,beta,lambda) - lasso_c(Xtilde,Ytilde,beta_update,lambda);
      beta = beta_update;
  }
  return beta;
}


// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  int length = lambda_seq.size();
  arma::mat beta_mat;
  beta_mat = beta_mat.zeros(Xtilde.n_cols,lambda_seq.size());
  arma::colvec beta_start;
  beta_start = beta_start.zeros(Xtilde.n_cols);
  for(int k = 1; k <= length; ++k){
    arma::colvec out = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[k-1], beta_start, eps = 0.001);
    arma::colvec beta_start = out;
    beta_mat.col(k-1) = out;
  }
  return beta_mat;
}





