# Clean Workspace ---------------------------------------------------------

rm(list = ls())


# Packages ----------------------------------------------------------------

library(hdm)
library(glmnet)
library(sandwich)
library(tidyverse)


# Working Directory -------------------------------------------------------

setwd("/Users/marcel/Desktop/RM_Econometrics_My_Part/Code")


# Random Seed -------------------------------------------------------------

set.seed(90922)


# Import Functions --------------------------------------------------------

source("functions.R")

# Simulation Function -------------------------------------------------------------

sim <- function(n, p, k, R2, alpha0, nreps) {
  res <- data.frame(matrix(NA, nreps, 12))
  colnames(res) <- c("d_ols", "se_ols", "d_pds", "se_pds",
                     "d_pds_min","se_pds_min", "d_pds_1se", "se_pds_1se",
                     "bias_ols", "bias_pds", "bias_pds_min", "bias_pds_1se")
  
  for (i in 1:nreps) {
    sim_data <- gen_design(n, p, k, R2, R2, alpha0)
    
    y <- sim_data$y
    X <- sim_data$X
    d <- sim_data$d
    
    ols <- lm(y ~ d + X)
    pds <- rlassoEffect(X, y, d, method = "double selection")
    pds_min <- pdl_cv(y,d,X,cv="min")
    pds_1se <- pdl_cv(y,d,X,cv="1se")
    res[i, 1] <- summary(ols)$coefficients["d", 1]
    res[i, 2] <- summary(ols)$coefficients["d", 2]
    res[i, 3] <- summary(pds)$coefficients["d1", 1]
    res[i, 4] <- summary(pds)$coefficients["d1", 2]
    res[i, 5] <- pds_min$alpha_hat
    res[i, 6] <- pds_min$se_hat
    res[i, 7] <- pds_1se$alpha_hat
    res[i, 8] <- pds_1se$se_hat
    res[i, 9] <- res[i, 1] - alpha0
    res[i, 10] <- res[i, 3] - alpha0
    res[i, 11] <- res[i, 5] - alpha0
    res[i, 12] <- res[i, 7] - alpha0
    
  }
  return(res)
}

# Simulation Function IV --------------------------------------------------

sim_iv <- function(n, p, k, alpha0, nreps) {
  res <- data.frame(matrix(NA, nreps, 12))
  colnames(res) <- c("d_ols", "se_ols", "d_pds", "se_pds",
                     "d_pds_min","se_pds_min", "d_pds_1se", "se_pds_1se",
                     "bias_ols", "bias_pds", "bias_pds_min", "bias_pds_1se")
  
  for (i in 1:nreps) {
    sim_data <- gen_design_iv(n, p, k, alpha0)
    
    y <- sim_data$y
    X <- sim_data$X
    d <- sim_data$d
    
    ols <- lm(y ~ d + X)
    pds <- rlassoEffect(X, y, d, method = "double selection")
    pds_min <- pdl_cv(y,d,X,cv="min")
    pds_1se <- pdl_cv(y,d,X,cv="1se")
    res[i, 1] <- summary(ols)$coefficients["d", 1]
    res[i, 2] <- summary(ols)$coefficients["d", 2]
    res[i, 3] <- summary(pds)$coefficients["d1", 1]
    res[i, 4] <- summary(pds)$coefficients["d1", 2]
    res[i, 5] <- pds_min$alpha_hat
    res[i, 6] <- pds_min$se_hat
    res[i, 7] <- pds_1se$alpha_hat
    res[i, 8] <- pds_1se$se_hat
    res[i, 9] <- res[i, 1] - alpha0
    res[i, 10] <- res[i, 3] - alpha0
    res[i, 11] <- res[i, 5] - alpha0
    res[i, 12] <- res[i, 7] - alpha0
    
  }
  return(res)
}

# Simulations -------------------------------------------------------------

res1 <- sim(n = 500, p = 200, k = 5, R2 = 0.5, alpha0 = 1, nreps = 100)
res2 <- sim(n = 500, p = 300, k = 5, R2 = 0.5, alpha0 = 1, nreps = 100)
res3 <- sim(n = 500, p = 400, k = 5, R2 = 0.5, alpha0 = 1, nreps = 100)
res4 <- sim(n = 500, p = 500, k = 5, R2 = 0.5, alpha0 = 1, nreps = 100)

res1_iv <- sim_iv(n = 500, p = 200, k = 5, alpha0 = 1, nreps = 100)
res2_iv <- sim_iv(n = 500, p = 300, k = 5, alpha0 = 1, nreps = 100)
res3_iv <- sim_iv(n = 500, p = 400, k = 5, alpha0 = 1, nreps = 100)
res4_iv <- sim_iv(n = 500, p = 500, k = 5, alpha0 = 1, nreps = 100)


# Plots -------------------------------------------------------------------

d_res1 <- data.frame(OLS = res1$d_ols, PDS = res1$d_pds, PDS_Min = res1$d_pds_min, PDS_1SE = res1$d_pds_1se)
d_res1_long <- stack(d_res1)
ggplot(data = d_res1_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res2 <- data.frame(OLS = res2$d_ols, PDS = res2$d_pds, PDS_Min = res2$d_pds_min, PDS_1SE = res2$d_pds_1se)
d_res2_long <- stack(d_res2)
ggplot(data = d_res2_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res3 <- data.frame(OLS = res3$d_ols, PDS = res3$d_pds, PDS_Min = res3$d_pds_min, PDS_1SE = res3$d_pds_1se)
d_res3_long <- stack(d_res3)
ggplot(data = d_res3_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res4 <- data.frame(OLS = res4$d_ols, PDS = res4$d_pds, PDS_Min = res4$d_pds_min, PDS_1SE = res4$d_pds_1se)
d_res4_long <- stack(d_res4)
ggplot(data = d_res4_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res1_iv <- data.frame(OLS = res1_iv$d_ols, PDS = res1_iv$d_pds, PDS_Min = res1_iv$d_pds_min, PDS_1SE = res1_iv$d_pds_1se)
d_res1_iv_long <- stack(d_res1_iv)
ggplot(data = d_res1_iv_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res2_iv <- data.frame(OLS = res2_iv$d_ols, PDS = res2_iv$d_pds, PDS_Min = res2_iv$d_pds_min, PDS_1SE = res2_iv$d_pds_1se)
d_res2_iv_long <- stack(d_res2_iv)
ggplot(data = d_res2_iv_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res3_iv <- data.frame(OLS = res3_iv$d_ols, PDS = res3_iv$d_pds, PDS_Min = res3_iv$d_pds_min, PDS_1SE = res3_iv$d_pds_1se)
d_res3_iv_long <- stack(d_res3_iv)
ggplot(data = d_res3_iv_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))

d_res4_iv <- data.frame(OLS = res4_iv$d_ols, PDS = res4_iv$d_pds, PDS_Min = res4_iv$d_pds_min, PDS_1SE = res4_iv$d_pds_1se)
d_res4_iv_long <- stack(d_res4_iv)
ggplot(data = d_res4_iv_long, mapping = aes(x = values, color = ind, fill = ind)) + 
  geom_density(alpha = 0.2) + geom_vline(xintercept = 1, linetype = "dashed",
                                         color = "red") + theme_bw() + scale_x_continuous(limits = c(0.75, 1.25))
