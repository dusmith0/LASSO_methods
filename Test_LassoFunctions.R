#### This script is to store my tests for the LassoFunctions

##Super simple set of data, note the Y-intercept should not be present. 
X <- matrix(c(rep(1,5),1:5),nrow=5)
beta <- c(2,3) 
sigma <- 3
Y <- X%*%beta

standardizeXY(X,Y)
