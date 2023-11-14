# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  # [ToDo] Center and scale X
  n <- nrow(X)
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, nrow(X), ncol(X), byrow = TRUE)
  weights <- apply(Xcentered,2,function(Xcentered) sqrt(crossprod((Xcentered),Xcentered) / n))
  #Xtilde = scale(X)* sqrt(n/(n-1))
  
  normsX <- colSums(Xcentered ^ 2)/n
  Xtilde <- Xcentered %*% diag(1/sqrt(normsX))

  #Xcentered <- scale(X,scale = FALSE)
  #weights <- apply(Xcentered,2,function(Xcentered) sqrt(crossprod((Xcentered),Xcentered) / n))
  #Xtilde <- scale(Xcentered, center = FALSE, scale = weights)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  ifelse(a > lambda,return(a - lambda),ifelse(a < -lambda, return(a + lambda), return(0)))
}

soft2 <- function(a,lambda){ #I would like to test which is faster.
  if(a > lambda){
    return(a - lambda)
  }else if(a < -lambda){
    return(a + lambda)
  }else{
    return(0)
  }
}

soft3 <- function(a,lambda){
  sign(a) * max(abs(a) - lambda,0)
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n = length(Ytilde)
  sum((Ytilde - Xtilde %*% beta) ^ 2)/(2 * n) + lambda * sum(abs(beta))
}

lasso2 <- function(Xtilde, Ytilde, beta, lambda){
  n = length(Ytilde)
  sum((Ytilde - Xtilde %*% beta) ^ 2)/(2 * n) + lambda * sum(abs(beta))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    stop(paste("It seems that your standardized and centered X and Y do not have the equivalent amount of rows."))
  }
  #[ToDo]  Check that lambda is non-negative
  if(lambda < 0){
    stop(paste("Warning: Please input a postive or zero penalty value."))
  }
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if(is.null(beta_start)){
    beta_start <- rep(0,ncol(X))
  }else if(ncol(Xtilde) != length(beta_start)){
    stop(paste("Warning: Please input a beta_start length that matched the number or columns in X"))
    
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  beta <- beta_update <- beta_start
  eps_check <- 100
  r <- Ytilde - Xtilde %*% beta_start
  
  while(abs(eps_check) > eps){
    
    for(j in 1:ncol(Xtilde)){
      beta_update[j] <- soft((beta[j] + crossprod(Xtilde[,j],r) / n),lambda)
      r <- r + X[,j] * (beta[j] - beta_update[j])
    }

    eps_check <- lasso(Xtilde,Ytilde,beta,lambda) - (fmin <- lasso(Xtilde,Ytilde,beta_update,lambda))
    beta <- beta_update
  }
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    stop(paste("It seems that your standardized and centered X and Y do not have the equivalent amount of rows."))
  }
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if(!is.null(lambda_seq)){
    lambda_seq <- sort(lambda_seq,decreacing = TRUE)
    lambda_seq <- lambda_seq[-which(lambda_seq <= 0)]
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if(is.null(lambda_seq)){
    lambda_max <- max(crossprod(Xtilde,Ytilde)/nrow(Xtilde))
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  beta_mat <- matrix(0,nrow = ncol(Xtilde), ncol = length(lambda_seq))
  fmin_vec <- rep(0,n_lambda)
  
  # Use warm starts strategy discussed in class for setting the starting values.
  beta_start <- NULL
  for(i in seq_along(lambda_seq)){
    new <- fitLASSOstandardized(Xtilde,Ytilde, lambda_seq[i], beta_start = beta_start,eps = eps)
    beta_start <- new$beta
    beta_mat[,i] <- new$beta
    fmin_vec[i] <- new$fmin
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  new <- standardizeXY(X,Y)
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  seq <- fitLASSOstandardized_seq(new$Xtilde,new$Ytilde,lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- seq$lambda_seq
  beta_mat <- seq$beta_mat
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  beta_original <- diag(1/sqrt(new$normsX)) %*% beta_mat
  beta_intercept <- mean(Y) - colSums(colMeans(X) * beta_mat)
  beta0_vec <- rbind(beta_intercept,beta_mat)
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  fitLASSO(X = X, Y = Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if(is.null(fold_ids)){
    fold_ids <- sample(1:n) %% k + 1
  }
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  errors <- matrix(NA,nrow = nrow(X),ncol = n_lambda)
  for(fold in 1:k){
    Xtrain <- X[fold_ids != fold, ]
    Ytrain <- Y[fold_ids != fold]
    
    Xtest <- X[fold_ids == fold,]
    Ytest <- Y[fold_ids == fold]
    
    out <- fitLASSO(X = Xtrain, Y = Ytrain, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
    errors[fold_ids == fold, fold] <- (Ytest - out$beta0_vec[1,] - crossprod(Xtest, out$beta0_vec[-1,])) ^ 2
  }
  cvm <- colMeans(errors)
  
  # [ToDo] Find lambda_min
  lambda_min <- min(cvm)

  # [ToDo] Find lambda_1SE
  
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

