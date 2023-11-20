
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")
# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
  # Test one: Super simple test to see if my function even works
  expect_equal(soft_c(3,2),1)
  expect_equal(soft_c(1,2),0)
  expect_equal(soft_c(-3,2),-1)

  # Test two using the LASSO objective
  n <- length(Ytilde)
  a <- sum((Ytilde - Xtilde %*% beta) ^ 2)/(2 * n) + lambda * sum(abs(beta))
  expect_equal(soft_c(a,.02),soft(a,.02))
  expect_equal(soft_c(a,.2),soft(a,.2))
  expect_equal(soft_c(a,.1),soft(a,.1))

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
  # Test one using LASSO Example 2 Data from the class notes.
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

  expect_equal(lasso_c(new$Xtilde,new$Ytilde,beta,lambda),lasso(new$Xtilde,new$Ytilde,beta,lambda))
     # No value was found, the functons are equal.

  # Test 2 using the RiboflavinData:
  library(hdi)
  data(riboflavin) 
  dim(riboflavin$x)
  class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]
  X = as.matrix(riboflavin$x)
  Y = riboflavin$y
  
  new <- standardizeXY(X,Y)
  lambda <- .04
  beta <- rep(.2,ncol(new$Xtilde))
  expect_equal(lasso_c(new$Xtilde,new$Ytilde,beta,lambda),lasso(new$Xtilde,new$Ytilde,beta,lambda))
    #Again no value was found so the functions are equal.
  
# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################
  ## Test one being done on the Example 2 data set from the notes. 
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
  beta_start <- rep(0,ncol(X))
  
  expected <- fitLASSOstandardized(new$Xtilde, new$Ytilde, lambda, beta_start = NULL, eps = 0.001)
  expect_equal(sum(fitLASSOstandardized_c(new$Xtilde, new$Ytilde, lambda, beta_start, eps = 0.001)),
               sum(expected$beta))
  
  ##Test 2 on the Riboflavin Data
  library(hdi)
  data(riboflavin) 
  dim(riboflavin$x)
  class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]
  X = as.matrix(riboflavin$x)
  Y = riboflavin$y
  
  new <- standardizeXY(X,Y)
  lambda <- .04
  beta_start <- rep(0,ncol(X))
  
  expected <- fitLASSOstandardized(new$Xtilde, new$Ytilde, lambda, beta_start = NULL, eps = 0.001)
  expect_equal(sum(fitLASSOstandardized_c(new$Xtilde, new$Ytilde, lambda, beta_start, eps = 0.001)),
               sum(expected$beta))
  
  # Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
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
  
  lambda_max <- max(crossprod(new$Xtilde,new$Ytilde)/nrow(new$Xtilde))
  lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = 30))
  
  expected <- fitLASSOstandardized_seq(new$Xtilde,new$Ytilde,lambda_seq = NULL, n_lambda = 30, eps = 0.001)
  expect_equal(sum(fitLASSOstandardized_seq_c(new$Xtilde,new$Ytilde,lambda_seq,eps = .001)),(sum(expected$beta_mat)))
  #This test is slightly off. The betas are exaclty the same unitl the 7th interation.
  #Then for some reason the dicimals are off at around the 5th decimal place. It seems to become more off as the iterations increase. 
  
  
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
