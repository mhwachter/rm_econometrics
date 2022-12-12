# tsls: computes 2sls estimate of coefficients b and
# variance covariance matrix VC assuming homoskedasticity for outcome 
# variable y where d are endogenous variables, in structural equation,
# x are exogensous variables in structural equation and z are 
# instruments.  x should include the constant term.
# 
# If the optional argument VG is supplied, the variance covariance matrix
# assuming that sqrt(n) times the unidentified regression coefficient on z 
# is a normal random variable with mean 0 and variance VG is computed as VC2.


tsls <- function(y,d,x,z,VG=NULL) {
  #browser()
  n <- length(y)
  a1 <- dim(d)[2]
  a2 <- dim(x)[2]
  if (is.null(x)) {a2 <- 0}
  if (is.vector(x)) {a2 <- 1}
  if (is.vector(d)) {a1 <- 1}
  k <- a1 + a2
  X <- cbind(d,x)
  Z <- cbind(z,x)
  
  Mxz <- t(X)%*%Z
  Mzz <- solve(t(Z)%*%Z)
  M <- solve(Mxz%*%Mzz%*%t(Mxz))
  
  b <- M%*%Mxz%*%Mzz%*%(t(Z)%*%y)
  e <- y - X%*%b
  VC1 <- (t(e)%*%e/(n-k))*M
  
 if(!is.null(VG)) {VC2 <- VC1 + M%*%Mxz%*%VG%*%t(Mxz)%*%M} else {VC2 <- NULL}
 return(list(b=b, VC1=VC1, VC2=VC2))
}

tsls2a <- function(y,d_pred,d,x,z,VG=NULL) {
  #browser()  
  n <- length(y)
  k <- dim(d_pred)[2] + dim(x)[2]
  if (is.null(x)) {k <- dim(d_pred)[2]}
  X <- cbind(d_pred,x)
  Z <- cbind(z,x)
  
  Mxz <- t(X)%*%Z
  Mzz <- solve(t(Z)%*%Z)
  M <- solve(Mxz%*%Mzz%*%t(Mxz))
  
  #b <- M%*%Mxz%*%Mzz%*%(t(Z)%*%y)
  b <- solve(t(d_pred)%*%d_pred)%*%t(d_pred)%*%y
  e <- y - X%*%b
  VC1 <- (t(e)%*%e/(n-k))*M
  
  if(!is.null(VG)) {VC2 <- VC1 + M%*%Mxz%*%VG%*%t(Mxz)%*%M} else {VC2 <- NULL}
  return(list(b=b, VC1=VC1, VC2=VC2))
} 

tsls2b <- function(X,Y,Z, Xhat) {
  #browser()
  n <- length(Y)
  K <- dim(X)[2]
  P <- Z%*%solve(t(Z)%*%Z)%*%t(Z)
  nom <- t(X)%*%P%*%Y
  den <- t(X)%*%Xhat
  b <- nom/den
  #c <- solve(t(Xhat)%*%Xhat)%*%t(Xhat)%*%Y
  s2 <- sum((Y-X*as.numeric(b))^2)/(n-K)
  VC1 <- s2*solve(t(X)%*%Xhat)
  return(list(b=b,VC1=VC1))
}
