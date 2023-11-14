#### This script is to store my tests for the LassoFunctions
library(microbenchmark)
##Super simple set of data, note the Y-intercept should not be present. 
X <- matrix(c(rep(1,5),1:5),nrow=5)
beta <- c(2,3) 
sigma <- 3
Y <- X%*%beta

standardizeXY(X,Y)


## Building a new data set for testing.
set.seed(120)
X <- c(sample(seq(1:25),25,replace=TRUE))
X <- matrix(X, ncol=5)
solve(X%*%t(X))

beta <- c(1,2,3,4,5)
Y <- X %*% beta + .4*rnorm(5)

new <- standardizeXY(X,Y)
n <- nrow(X)
sqrt(t(new$Xtilde) %*% new$Xtilde / n)

new$Xtilde

lasso(new$Xtilde,new$Ytilde,beta,.4)

##Using the exmaple data to see if I am receiving the same output so far.
#++++ Start
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

new <- standardizeXY(X,Y)
n <- nrow(X)
(t(new$Xtilde) %*% new$Xtilde / n)
edit((t(new$Xtilde) %*% new$Xtilde / n))
#++++ Ends
#++++ For ease of inputing into functions
Xtilde <- new$Xtilde
Ytilde <- new$Ytilde

##Example 2 centered calculations
meanY = mean(Y)
Ytilde = Y - mean(Y)
meansX = colMeans(X)

Xcentered = X - matrix(meansX, nrow(X), ncol(X), byrow = T)

normsX = colSums(Xcentered^2)/n
Xtilde = Xcentered %*% diag(1/sqrt(normsX))

lasso(Xtilde,Ytilde,beta,.4)

##Testing to see if my functions actually run.
lambda <- .5
fitLASSOstandardized(new$Xtilde,new$Ytilde,lambda,beta_start = NULL, eps = .001)
fitLASSOstandardized_seq(Xtilde,Ytilde,lambda_seq=NULL,n_lambda = 60,eps=.001)


##checking to ensure that I actually end with zeroed betas at a large lambda.
fitLASSOstandardized(new$Xtilde,new$Ytilde,2,beta_start = NULL, eps = .001)

microbenchmark( #3.5 microseconds
  lasso(new$Xtilde,new$Ytilde,beta,lambda=.4),times=200L
)

microbenchmark( #3.5 microseconds  This one contained the as.numeric(crossprod trick. It is still slower.)
  lasso2(new$Xtilde,new$Ytilde,beta,lambda=.4),times=200L
)

lasso2 <- function(Xtilde, Ytilde, beta, lambda){
  n = length(Ytilde)
  as.numeric(crossprod(Ytilde - Xtilde %*% beta))/(2 * n) + lambda * sum(abs(beta))
}

##Checking to see if fitLASSO runs.
microbenchmark(
fitLASSO(X,Y)
)

##Using Example 2 with my LASSO Update
beta_start = NULL
eps = 0.001
lambda <- 1.44

if(is.null(beta_start)){
  beta_start <- rep(0,ncol(X))
}else if(ncol(Xtilde) != length(beta_start)){
  stop(paste("Warning: Please input a beta_start length that matched the number or columns in X"))
  
}

beta <- beta_update <- beta_start
eps_check <- 100
r <- new$Ytilde - new$Xtilde %*% beta_start

microbenchmark( ##13.7627 milliseconds
while(abs(eps_check) > eps){
  
  for(j in 1:ncol(new$Xtilde)){
    beta_update[j] <- soft((beta[j] + (t(X[,j]) %*% r) / n),lambda)
    r <- r + X[,j]*(beta[j] - beta_update[j])
  }
  eps_check <- lasso(new$Xtilde,new$Ytilde,beta,lambda) - (fmin <- lasso(new$Xtilde,new$Ytilde,beta_update,lambda))
  beta <- beta_update
},
times = 200L
)

microbenchmark( ##12.4169 milliseconds
  while(abs(eps_check) > eps){
    
    for(j in 1:ncol(new$Xtilde)){
      beta_update[j] <- soft((beta[j] + crossprod(X[j],r) / n),lambda)
      r <- r + X[,j]*(beta[j] - beta_update[j])
    }
    eps_check <- lasso(new$Xtilde,new$Ytilde,beta,lambda) - (fmin <- lasso(new$Xtilde,new$Ytilde,beta_update,lambda))
    beta <- beta_update
  },
  times = 2000L
)

microbenchmark( ##12.5669 milliseconds
  while(abs(eps_check) > eps){
    
    for(j in 1:ncol(new$Xtilde)){
      beta_update[j] <- soft3((beta[j] + crossprod(X[j],r) / n),lambda)
      r <- r + X[,j]*(beta[j] - beta_update[j])
    }
    eps_check <- lasso(new$Xtilde,new$Ytilde,beta,lambda) - (fmin <- lasso(new$Xtilde,new$Ytilde,beta_update,lambda))
    beta <- beta_update
  },
  times = 2000L
)


##Tyring to see what the Riboflavin data does.
new <- standardizeXY(X,Y)
nrow(new$Ytilde)


##Checking speeds on my standardizeXY
microbenchmark(#seems to be taking around 30 milliseconds.
  standardizeXY(X,Y),
  times = 20L
)
