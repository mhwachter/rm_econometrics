rm(list=ls())

# This Simulation tries to show how Lasso outperforms OLS under a sparse
# Setup using MSE as performance indicator. 

# Load relevant packages
library(codetools) # ??
library(MASS) # Generate multivariate normal distributions
library(mvtnorm) # Multivariate DGP (using variance-covariance matrix)
library(ggplot2) # Graphs
library(dplyr)
library(glmnet) # Lasso
library(Rcpp) # Lasso
library(writexl)
library(corrplot)



# FUNCTION: Variance-covariance Matrix 
varcov_matrix = function(covariates, sd){
  matrix_vcov = diag(sd^2,covariates,covariates)
  return(matrix_vcov)
}

# FUNCTION: Toeplitz variance-covariance Matrix
varcov_matrix_toeplitz = function(rho, covariates,sigma_grid){
  lim_sup = covariates-1
  matrix_vcov = toeplitz(rho^(0:lim_sup))
  diag(matrix_vcov) = runif(n=covariates,min=sigma_grid[1],max=sigma_grid[2])
  return(matrix_vcov)
}


length_grid = 100
grid = 10^seq(-2, 2, length.out=length_grid)
#grid_sparse = seq(0,100,2)
grid_sparse = seq(10,100,10)

#d=100 # Number of covariates
N=1000
reps = 500
#Containers
best_sim_MSE_lasso_sparse = matrix(NA, nrow=reps, ncol=length(grid_sparse))
MSE_lasso_sim_sparse = matrix(NA, nrow=reps, ncol=length(grid_sparse))
#MSE_lasso_sim_sparse = rep(NA,100)
MSE_OLS_sim_sparse = matrix(NA, nrow=reps, ncol=length(grid_sparse))
#MSE_OLS_sim_sparse = rep(NA,100)


set.seed(123)

for(i in grid_sparse){ # grid sparse is a sequence from 0 to 100 by 2
  
  length = length_grid - i
  c = i/10
  
  #cat("(Sparse Regressors,Iteration)",i,r)
  #cat("\n")
  
  
  for(r in 1:reps){
    
    cat("(Sparse Regressors,Iteration)",i,r)
    cat("\n")
    #Define DGP
    cov_sim_sparse = varcov_matrix(covariates=length_grid, sd=1)
    #length = d-i
    true_coeff_sim1 = seq(from=0.1, to=0.5, length.out=length)
    true_coeff_sparse = append(true_coeff_sim1, rep(0,i)) #Appending zeros gradually increases sparsity
    true_means_sim_sparse = rep(0,length_grid)
    
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
    best_sim_MSE_lasso_sparse[r]= lasso_cv_sim_sparse$lambda.min
    
    ##predict Y hat and use to calculate residuals and mse
    lasso_preds_sim_sparse = predict(lasso_cv_sim_sparse, s = "lambda.min", newx = X_sim_sparse_test)
    MSE_lasso_sim_sparse[r,c]= mean((lasso_preds_sim_sparse - Y_sim_sparse_test)^2)
    

    #Compute OLS baseline for comparison
    OLS_cv_sim_sparse<-glmnet(X_sim_sparse, Y_sim_sparse, alpha = 0, lambda = 0,family="gaussian")
    OLS_preds_sim_sparse<-predict(OLS_cv_sim_sparse, newx=X_sim_sparse_test)
    MSE_OLS_sim_sparse[r,c]=mean((OLS_preds_sim_sparse - Y_sim_sparse_test)^2)
  }
}  


# Extract the results
# Sample size n=1000
mean_MSE_lasso_sim_sparse = colMeans(MSE_lasso_sim_sparse)
mean_MSE_OLS_sim_sparse = colMeans(MSE_OLS_sim_sparse)

# Plot the results for sample size N=1000
Index<-seq(10,100,10)
dfsim2<-data.frame(Index, mean_MSE_lasso_sim_sparse, mean_MSE_OLS_sim_sparse)
p2<-ggplot(dfsim2, aes(Index)) + 
  geom_point(aes(y = mean_MSE_lasso_sim_sparse, colour = "Lasso")) +
  geom_point(aes(y = mean_MSE_OLS_sim_sparse, colour = "OLS")) +
  xlab("Number of Coefficients with true 0")+ylab("MSE")
p2

################################################################################
###################### WRITE REULTS TO EXCEL FILE ##############################
################################################################################
lasso_results_df = as.data.frame(MSE_lasso_sim_sparse)
OLS_results_df = as.data.frame(MSE_OLS_sim_sparse)
best_sim_MSE_lasso_sparse_df = as.data.frame(best_sim_MSE_lasso_sparse)
write_xlsx(lasso_results_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim21_increasing_sparsity_no_multicollinearity_lasso_results.xlsx")
write_xlsx(OLS_results_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim21_increasing_sparsity_no_multicollinearity_OLS_results.xlsx")
write_xlsx(best_sim_MSE_lasso_sparse_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim21_increasing_sparsity_no_multicollinearity_best_sim_mse_lasso_results.xlsx")