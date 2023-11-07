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
