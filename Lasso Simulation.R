# THIS CODE FOLLOWS FROM COMPUTATIONAL STATISTICS LECTURE SUMMER 2021 PS5 Q2 #

# Useful packages
library(codetools)
library(MASS) # Generate multivariate normal distributions
library(mvtnorm)
library(ggplot2)
library(glmnet)


# Data generating process (DGP)

# Number of observations
n = 1000

# Number of covariates
p = 50

# FUNCTION: Variance-covariance Matrix 
varcov_matrix = function(covariates, sd){
  varcov_matrix = diag(sd^2,covariates,covariates)
  return(varcov_matrix)
}


# True means 
true_means = rep(0,p)

# True coefficients
true_coeff = runif(p,0.1,0.5) # uniform distribution length p from 0.1 to 0.5

# Generate variance-covariance matrix
sd_diag = 1
varcov = varcov_matrix(covariates = p, sd = sd_diag )

# Generate observations from multivariate istribution
X = mvrnorm(n=n,mu=true_means,Sigma=varcov)

# Error terms
error = rnorm(n=n,mean=0,sd=1)

# Outcome variable
Y = X %*% true_coeff + error

# Grid of penalties
grid = 10^seq(2,-2,length.out=100)

###################################################
############## MONTE CARLO SIMULATION #############
##################################################
## SIMULATION 1)

set.seed(123)

# Containers
best_sim_mse_lasso = rep(NA,100) #Store lambda that minimizes MSE
MSE_lasso_sim = rep(NA,100) 
MSE_OLS_sim = rep(NA,100)

# Sample size 
n = 100

for(s in 2:150){ # this interval can be predefined (avoid hard coding)
  #Define DGP
  #First generate training sample
  cov_sim = varcov_matrix(covariates=s, sd=1)
  true_coeff_sim = seq(from=0.1, to=0.5, length.out=s)
  true_means_sim = rep(0,s)
  X_sim = mvrnorm(n=n, mu=true_means_sim, Sigma=cov_sim)
  error = rnorm(n=n,sd=1)
  Y_sim = X_sim %*% true_coeff_sim + error
  
  
  #Now generate test sample
  X_sim_test = mvrnorm(n=n,mu=true_means_sim, Sigma=cov_sim)
  error_test = rnorm(n=n, sd=1)
  Y_sim_test = X_sim_test %*% true_coeff_sim + error_test
  
  #Fit Lasso with 10 fold cross validation
  lasso_sim = glmnet(x=X_sim, y=Y_sim,family="gaussian", lambda=grid, alpha=1)
  lasso_cv_sim = cv.glmnet(X_sim, Y_sim, alpha=1, lambda=grid)
  #Select lambda that minimizes MSE
  best_sim_mse_lasso = lasso_cv_sim$lambda.min
  #Use predict function to compute predicted Y hat values based on fitted model
  lasso_preds_sim = predict(lasso_sim, s=best_sim_mse_lasso, newx=X_sim_test)
  mse_lasso_sim[s] = mean((lasso_preds_sim - Y_sim_test)^2)
  
  #Compute OLS baseline for comparison
  OLS_cv_sim = glmnet(X_sim, Y_sim, alpha=0, lambda=0, family="gaussian")
  OLS_preds_sim = predict(OLS_cv_sim,newx=X_sim_test)
  MSE_OLS_sim[s]=mean((OLS_preds_sim - Y_sim_test)^2)
  
}

# Plot results
index = seq(1,150,1)
dfs = data.frame(index,MSE_lasso_sim,MSE_OLS_sim)

p = ggplot(dfs,aes(index)) +
  geom_point(aes(y=MSE_lasso_sim, colour="Lasso"))
+ geom_point(aes(y=MSE_OLS_sim, colour="OLS"))
+ xlab("Number of Coefficients")
+ ylab("MSE")
p + coord_cartesian(xlim=c(0,50),ylim=c(0,5))
p


###################################################
############## MONTE CARLO SIMULATION #############
##################################################
## SIMULATION 2)

grid_sparse = seq(0,100,2)
d=100
N=1000
#Containers
best_sim_MSE_lasso_sparse = rep(NA,100)
MSE_lasso_sim_sparse = rep(NA,100)
MSE_OLS_sim_sparse = rep(NA,100)

for(i in grid_sparse){
  #Define DGP
  cov_sim_sparse = varcov_matrix(covariates=d, sd=1)
  length = i-d
  true_coeff_sim1 = seq(from=0.1, to=0.5, length.out=length)
  true_coeff_sparse = append(true_coeff_sim1, rep(0,i)) #Appending zeros gradually increases sparsity
  true_means_sim_sparse = rep(0,d)
  
  #Generate training data
  X_sim_sparse = mvrnorm(n=N, mu=true_means_sim_sparse, Sigma=cov_sim_sparse)
  error_sparse = rnorm(n=N,sd=1)
  Y_sim_sparse = X_sim_sparse %*% true_coeff_sparse + error_sparse
  
  #Generate test data
  X_sim_sparse_test = mvrnorm(n=N, mu=true_means_sim_sparse, Sigma=cov_sim_sparse)
  error_sparse_test = rnorm(n=N, sd=1)
  Y_sim_sparse_test = X_sim_sparse_test %*% true_coeff_sparse + error_sparse

  #Fit lasso model
  #lasso.sim.sparse<-glmnet(x=X.sim.sparse,y=Y.sim.sparse,family="gaussian", lambda=grid,alpha=1)
  lasso_cv_sim_sparse = cv.glmnet(X_sim_sparse, Y_sim_sparse, alpha = 1, lambda = grid)
  #Select lambda that minimizes mse using 10 fold cross validation
  best_sim_MSE_lasso_sparse[i]= lasso_cv_sim_sparse$lambda.min
  ##predict Y hat and use to calculate residuals and mse
  lasso_preds_sim_sparse = predict(lasso_cv_sim_sparse, s = "lambda.min", newx = X_sim_sparse_test)
  MSE_lasso_sim_sparse[i]= mean((lasso_preds_sim_sparse - Y_sim_sparse_test)^2)
  
  ##fit ridge model
  #ridge.sim.sparse<-glmnet(x=X.sim.sparse,y=Y.sim.sparse,family="gaussian", lambda=grid,alpha=0)
  #ridge.cv.sim.sparse <- cv.glmnet(X.sim.sparse, Y.sim.sparse, alpha = 0, lambda = grid)
  ##find lambda that minimizes mse using 10 fold cross validation
  #best.sim.mse.ridge.sparse[i]= ridge.cv.sim.sparse$lambda.min
  if (i>=98){
    #plot(ridge.cv.sim.sparse)
    plot(lasso_cv_sim_sparse)
  }
  ##predict Y hat and use to calculate residuals and test mse
  #ridge.preds.sim.sparse = predict(ridge.cv.sim.sparse, s ="lambda.min", newx = X.sim.sparse.test)
  #print(best.sim.mse.ridge.sparse)
  #mse.ridge.sim.sparse[i]= mean((ridge.preds.sim.sparse - Y.sim.sparse.test)^2)
  
  OLS_cv_sim_sparse<-glmnet(X_sim_sparse, Y_sim_sparse, alpha = 0, lambda = 0,family="gaussian")
  OLS_preds_sim_sparse<-predict(OLS_cv_sim_sparse, newx=X_sim_sparse_test)
  MSE_OLS_sim_sparse[i]=mean((OLS_preds_sim_sparse - Y_sim_sparse_test)^2)
}


















