# This function computes the Least Squares parameters
# with a penalty on the L1-norm of the parameters
#
# Method used:
# The Shooting method of [Fu, 1998]
#
# Modifications:
# We precompute the Hessian diagonals, since they do not 
# change between iterations

LassoShooting2 <- function(X,y, lambda, control=list(maxIter=10000, optTol=10^(-5), zeroThreshold=10^(-4))) {
  #browser()
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  # Start from the LS solution for beta
  beta <- solve(XX+diag(as.vector(lambda))%*%diag(1,p))%*%Xy
  # Start the log
  w_old <- beta
  k <-1
  wp <- beta
  m <- 0
  XX2 <- XX*2
  Xy2 <- Xy*2
  
  while (m < control$maxIter) {
  beta_old <- beta    
  for (j in 1:p) {
    # Compute the Shoot and Update the variable
    S0 <- sum(XX2[j,]*beta) - XX2[j,j]*beta[j] - Xy2[j]
    if (S0 > lambda[j]) beta[j] <- (lambda[j] - S0)/XX2[j,j]
    if (S0 < -1*lambda[j]) beta[j] <- (-1*lambda[j] - S0)/XX2[j,j]
    if (abs(S0) <= lambda[j]) beta[j] <- 0
    }
  m <- m+1
  # Update the log
  w_old <- beta
  k <- k+1
  wp <- cbind(wp,beta)
  # Check termination
  if (sum(abs(beta-beta_old)) < control$optTol) {break}
  }
  
  w <- beta
  return(list(w=w, wp=wp, m=m))
}

summary_lasso <- function(lasso, beta=NULL, y, X, threshold=10^(-5)) {
  MSE_obs <- mean((y-X%*%lasso$w)^2)
  MAE_obs <- mean(abs(y-X%*%lasso$w))
  beta_zero <- sum(abs(lasso$w)<threshold)
  MSE_beta <- mean((beta-lasso$w)^2)
  MAE_beta <- mean(abs(beta-lasso$w))
  return(list(MSE_obs= MSE_obs, MAE_obs= MAE_obs, beta_zero=beta_zero, 
              MSE_beta= MSE_beta, MAE_beta= MAE_beta))
}

# input: object of type Lassoshooting 2, X,y
# output lm object

post_lasso <- function(lasso, X, y, intercept="FALSE"){
  ind <- (abs(lasso$w) > 0)
  X0 <- X[,ind]
  if (intercept==FALSE) {
  post_lasso <- lm(y ~ -1+ X0)
  } else {
    post_lasso <- lm(y ~ X0) 
  }
  return(post_lasso)
}


