# General Design ----------------------------------------------------------

gen_design <- function(n, p, k, RsqY, RsqD, alpha0) {
  tbeta <- tgamma <- rep(0,p)
  tbeta[1:k] <- tgamma[1:k] <- 1
  
  cD <- sqrt(RsqD/(k-RsqD*k))
  
  if (alpha0==0){
    cY <- sqrt(RsqY/(k-RsqY*k))
  }
  if (alpha0!=0){
    a   <- k*(RsqY-1)
    b   <- 2*alpha0*cD*k*(RsqY-1)
    c   <- alpha0^2*cD^2*k*(RsqY-1)+(alpha0^2+1)*RsqY
    cY  <- max(Re(polyroot(c(c,b,a))))
  }
  
  gamma <-  cD*tgamma
  beta  <-  cY*tbeta
  
  X <- matrix(rnorm(n*p),n,p)
  d <- X%*%gamma + rnorm(n)
  y <- alpha0*d + X%*%beta + rnorm(n)
  
  return(list(X=X,d=d,y=y))
}

# General Design IV -------------------------------------------------------

gen_design_iv <- function(n, p, k, alpha0) {
  tbeta <- tgamma <- rep(0,p)
  tbeta[1:k] <- tgamma[1:k] <- 1
  
  gamma <-  tgamma
  beta  <-  tbeta
  
  X <- matrix(rnorm(n*p),n,p)
  d <- X%*%gamma + rnorm(n)
  y <- alpha0*d + rnorm(n)
  
  return(list(X=X,d=d,y=y))
}

# PDS with CV -------------------------------------------------------------

pdl_cv <- function(y,d,X,cv){
  
  cvfit1  <- cv.glmnet(X,d,nfolds=5)
  cvfit2  <- cv.glmnet(X,y,nfolds=5)
  
  if (cv=="min"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.min")[-1])
    pi      <- as.matrix(coef(cvfit2, s = "lambda.min")[-1])
  }
  if (cv=="1se"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.1se")[-1])
    pi      <- as.matrix(coef(cvfit2, s = "lambda.1se")[-1])
  }
  
  ind <- (abs(gamma) + abs(pi) > 0)
  
  if (sum(ind)==0) {
    obj <- lm(y ~ d)
  }
  if (sum(ind)>0){
    Xsel  <- X[,ind]
    obj   <- lm(y ~ d + Xsel)
  } 
  
  vcov      <- vcovHC(obj, type = "HC1")
  robust_se <- sqrt(diag(vcov))
  
  alpha_hat <- obj$coefficients[2]
  se_hat    <- robust_se[2]
  
  return(list(alpha_hat=alpha_hat,se_hat=se_hat,sel_index=ind))
  
}