
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
# Test one: Super simple test to see if my function even works
soft_c(3,2)
soft_c(1,2)
soft_c(-3,2)


# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
# Test one using LASSO Example 2 Data.
n = 50
beta = c(1, 0.5)
beta0 = 2 
p = length(beta)
sigma = 0.4 

library(mnormt)
set.seed(983645) 

Sigma = matrix(0.7, p, p) + diag(rep(1-0.7, p)) 
X = rmnorm(n, mean = rep(0, p), varcov = Sigma) 

Y = beta0 + X %*% beta + sigma * rnorm(n)
lambda <- .4 
new <- standardizeXY(X,Y)

lasso_c(new$Xtilde,new$Ytilde,beta,lambda)
lasso(new$Xtilde,new$Ytilde,beta,lambda)

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
