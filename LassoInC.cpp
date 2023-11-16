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
  int n = Ytilde.size();
  arma::colvec beta; //= beta_start;
  arma::colvec beta_update;// = beta_start;
  double eps_check = -.4456;
  arma::colvec r = Ytilde - Xtilde * beta_start;
  
  while (fabs(eps_check) > eps){
    
    for(int j = 1; j <= n; ++j){
      beta_update(j) = soft_c((beta(j) + Xtilde.col(j-1).t()*r / n),lambda);
      r = r + Xtilde.col(j-1) * (beta(j) - beta_update(j));
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
}





