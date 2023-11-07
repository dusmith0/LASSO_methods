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

##Example 2 centered calculations
meanY = mean(Y)
Ytilde = Y - mean(Y)
meansX = colMeans(X)

Xcentered = X - matrix(meansX, nrow(X), ncol(X), byrow = T)

normsX = colSums(Xcentered^2)/n
Xtilde = Xcentered %*% diag(1/sqrt(normsX))

lasso(Xtilde,Ytilde,beta,.4)
