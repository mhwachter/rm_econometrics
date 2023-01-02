rm(list=ls())

# This Simulation tries to show how Lasso outperforms OLS under a sparse
# Setup using MSE as performance indicator. 

# Load relevant packages
library(codetools) # ??
library(MASS) # Generate multivariate normal distributions
library(mvtnorm) # Multivariate DGP (using variance-covariance matrix)
library(ggplot2) # Graphs
library(glmnet) # Lasso
library(Rcpp) # Lasso
library(writexl)


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

# Predefined parameters 
# Sample size 
sample_sizes = 100

# Simulation repetitions
rep = 500

# Number of regressors
#reg_quan = c(50,100)
reg_quan = seq(2,100,1)

# Grid of penalties 
length_grid = 100
grid = 10^seq(-2, 2, length.out=length_grid)

# Containers
MSE_lasso_sim = matrix(NA, nrow=rep, ncol=length(reg_quan))
MSE_OLS_sim = matrix(NA, nrow=rep, ncol=length(reg_quan))
best_sim_mse_lasso = matrix(NA, nrow=rep, ncol=length(reg_quan))
#MSE_lasso_sim = list(matrix(NA, nrow=rep, ncol=length(sample_sizes)),matrix(NA, nrow=rep, ncol=length(sample_sizes)))
#MSE_OLS_sim = list(matrix(NA, nrow=rep, ncol=length(sample_sizes)),matrix(NA, nrow=rep, ncol=length(sample_sizes)))
################################################################################
####################### MONTE CARLO SIMULATION #################################
################################################################################
#container = rep(NA,length(sample_sizes))
set.seed(123)



for (r in 1:length(reg_quan)){
  regressors = reg_quan[r]
  cat("Number Regressors:",regressors)
  cat("\n")
  
  for (i in 1:rep){
    #cat("n:", sample)
    #cat("\n")
    #cat("r:", regressors)
    #cat("\n")
    #cat("i:",i)
    #cat("\n")
    sample = sample_sizes
    #Define DGP
    #First generate training sample
    cov_sim = varcov_matrix(covariates=regressors, sd=1)
    true_coeff_sim = seq(from=0.1, to=0.5, length.out=regressors)
    true_means_sim = rep(0,regressors)
    X_sim = mvrnorm(n=sample, mu=true_means_sim, Sigma=cov_sim)
    error = rnorm(n=sample,sd=1)
    Y_sim = X_sim %*% true_coeff_sim + error
    
    #Now generate test sample
    X_sim_test = mvrnorm(n=sample,mu=true_means_sim, Sigma=cov_sim)
    error_test = rnorm(n=sample, sd=1)
    Y_sim_test = X_sim_test %*% true_coeff_sim + error_test
    
    #Fit Lasso with 10 fold cross validation
    #lasso_sim = glmnet(x=X_sim, y=Y_sim,family="gaussian", lambda=grid, alpha=1)
    lasso_cv_sim = cv.glmnet(X_sim, Y_sim, type.measure="mse", alpha=1, lambda=grid, family ="gaussian")
    
    #Select lambda that minimizes MSE
    best_sim_mse_lasso[i,r] = lasso_cv_sim$lambda.min
    
    #Use predict function to compute predicted Y hat values based on fitted model
    lasso_preds_sim = predict(lasso_cv_sim, s=lasso_cv_sim$lambda.min , newx=X_sim_test)
    MSE_lasso_sim[i,r] = mean((lasso_preds_sim - Y_sim_test)^2)
    
    #Compute OLS baseline for comparison
    OLS_cv_sim = glmnet(X_sim, Y_sim, type.measure="mse", alpha=0, lambda=0, family="gaussian")
    OLS_preds_sim = predict(OLS_cv_sim,newx=X_sim_test)
    MSE_OLS_sim[i,r]=mean((OLS_preds_sim - Y_sim_test)^2)
    
  }
  
}


# Extract the results
# Sample size n=100
mean_MSE_lasso_sim_100 = colMeans(MSE_lasso_sim)
mean_MSE_OLS_sim_100 = colMeans(MSE_OLS_sim)

#PLOT THE RESULTS FOR SAMPLE SIZE n=100
index = seq(1,length(reg_quan),1)
dfs = data.frame(index,mean_MSE_lasso_sim_100,mean_MSE_OLS_sim_100)

p = ggplot(dfs,aes(index)) +
  geom_point(aes(y=mean_MSE_lasso_sim_100, colour="Lasso"))+ geom_point(aes(y=mean_MSE_OLS_sim_100, colour="OLS"))+ xlab("Number of Coefficients")+ ylab("MSE")
p + coord_cartesian(xlim=c(0,50),ylim=c(0,5))
p


################################################################################
###################### WRITE REULTS TO EXCEL FILE ##############################
################################################################################
lasso_results_df = as.data.frame(MSE_lasso_sim)
OLS_results_df = as.data.frame(MSE_OLS_sim)
best_sim_mse_lasso_df = as.data.frame(best_sim_mse_lasso)
write_xlsx(lasso_results_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim11_increasing_regressors_no_multicollinearity_lasso_results.xlsx")
write_xlsx(OLS_results_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim11_increasing_regressors_no_multicollinearity_OLS_results.xlsx")
write_xlsx(best_sim_mse_lasso_df, "C:\\Users\\57300\\Documents\\Cristian Gutierrez\\Uni - Bonn\\5. Semester\\RM ECONOMETRICS\\PROJECT\\R Simulations\\Excel results\\sim11_increasing_regressors_no_multicollinearity_best_sim_mse_lasso_results.xlsx")

