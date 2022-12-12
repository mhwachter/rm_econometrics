hetero_se <- function(x,e,XpXinv) {
  k <- dim(x)[2]
  n <- dim(x)[1]
  
  V <- t((x*((e^2)%*%t(rep(1,k)))))%*%x
  
  vhetero <- n/(n-k)*XpXinv%*%V%*%t(XpXinv)
  se <- sqrt(diag(vhetero))
  
  return(list(se=se, vhetero=vhetero, V=V))
} 