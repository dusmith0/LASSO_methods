#### This script is to store my tests for the LassoFunctions

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
sqrt(t(new$Xtilde) %*% new$Xtilde / n)
#++++ Ends

##Example 2 centered calculations
meanY = mean(Y)
Ytilde = Y - mean(Y)
meansX = colMeans(X)

Xcentered = X - matrix(meansX, nrow(X), ncol(X), byrow = T)

normsX = colSums(Xcentered^2)/n
Xtilde = Xcentered %*% diag(1/sqrt(normsX))

lasso(Xtilde,Ytilde,beta,.4)


##Using Example 2 with my LASSO Update
beta_start = NULL
eps = 0.001
lambda <- 3

if(is.null(beta_start)){
  beta_start <- rep(0,ncol(X))
}else if(ncol(Xtilde) != length(beta_start)){
  stop(paste("Warning: Please input a beta_start length that matched the number or columns in X"))
  
}

beta <- beta_update <- beta_start
eps_check <- 100
while(eps_check > eps){
  
  r <- new$Ytilde - new$Xtilde %*% beta_start
  for(j in 1:ncol(new$Xtilde)){
    beta_update[j] <- soft((beta[j] + t(X[,j]) %*% r / n),lambda)
    r <- r + X[,j]*(beta[j] - beta_update[j])
  }
  eps_check <- lasso(new$Xtilde,new$Ytilde,beta,lambda) - (fmin <- lasso(new$Xtilde,new$Ytilde,beta_update,lambda))
  beta <- beta_update
  #fmin <- lasso(beta)
}
