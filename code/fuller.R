# [b VC VCb VCc] = fuller(y,d,x,z,c) computes FULLER estimate of coefficients b and
# variance covariance matrix using usual asymptotic approximation (VC) and
# Bekker asymptotic approximation (VCb) assuming homoskedasticity for outcome 
# variable y where where d are endogenous variables, in structural equation,
# x are exogensous variables in structural equation and z are 
# instruments.  x should include the constant term.  c is a user specified
# parameter (c = 1 higher order unbiased; c >= 4 higher order admissible
# under quadratic loss)

fuller <- function(y,d,x=NULL,z,c) {
 #browser()
 n <- length(y)
 a1 <- dim(d)[2]
 a2 <- dim(x)[2]
 if (is.null(x)) {a2 <- 0}
 if (is.vector(x)) {a2 <- 1}
 if (is.vector(d)) {a1 <- 1}
 k <- a1 + a2
 
 #k <- dim(x)[2] + dim(d)[2]
 #if (!is.matrix(z)) {z <- as.matrix(z)}
 #if (is.null(x)) {k <- dim(d)[2]}
 if (!is.matrix(z)) {z <- as.matrix(z)}
 
 if(!is.null(x)) {
   Mxinv <- solve(t(x)%*%x)
   My <- y - x%*%Mxinv%*%(t(x)%*%y)
   Md <- d -x%*%Mxinv%*%(t(x)%*%d)
   Mz <- solve(t(z)%*%z-(t(z)%*%x)%*%Mxinv%*%(t(x)%*%z))
 } else {
   Mxinv <- NULL
   My <- y
   Md <- d
   Mz <- solve(t(z)%*%z)
 }
 
 M <- solve(t(cbind(My, Md))%*%cbind(My, Md))%*%((t(cbind(My, Md))%*%z)%*%Mz%*%(t(z)%*%cbind(My, Md)))
 alpha <- min(eigen(M)$values)
 alpha1 <-  (alpha - (1 - alpha)*c/(n - k - dim(z)[2]))/(1 - (1 - alpha)*c/(n - k - dim(z)[2]))    

 X     <- cbind(d,x)
 Z     <- cbind(z,x)
 Mxy   <- t(X)%*%y
 Mzy   <- t(Z)%*%y
 Mzx   <- t(Z)%*%X
 Mxx   <- t(X)%*%X
 Mzz   <- solve(t(Z)%*%Z)
 Mxzzx <- t(Mzx)%*%Mzz%*%Mzx

 H     <- Mxzzx - alpha1*Mxx
 b     <- solve(H)%*%(t(Mzx)%*%Mzz%*%Mzy - alpha1*Mxy) #error??
 e     <- y - X%*%b
                            
 J     <- Mxzzx - alpha1*(t(X)%*%e)%*%(t(e)%*%X)%*%solve((t(e)%*%e))
 S     <- (1 - alpha1)*J - alpha1*H
            
 VC    <- (t(e)%*%e/(n - k))*solve(H)
 VCb   <- (t(e)%*%e/(n - k))%*%solve(H)%*%S%*%solve(H)

 s2 <- (t(e)%*%e)/(n-k)
 atilde <- (t(e)%*%Z)%*%Mzz%*%(t(Z)%*%e)%*%solve((t(e)%*%e))
 Ups <- Z%*%Mzz%*%(t(Z)%*%X)
 Xhat <- X - e%*%(t(e)%*%X)%*%solve((t(e)%*%e))
 Vhat <- Xhat - Z%*%Mzz%*%(t(Z)%*%Xhat)
 tau <- k/n
 kappa <- 0
 A1 <- 0
 A2 <- 0
 B1 <- 0
 for (i in 1:n) {
  pii <- Z[i,]%*%Mzz%*%Z[i,] # dimension??
  kappa <- kappa + pii
  A1 <- A1+(pii-tau)%*%Ups[i,]
  A2 <- A2+(e[i]^2)*Vhat[i,]/n
  B1 = B1+((e[i]^2-s2)*(t(Vhat[i,])%*%Vhat[i,]))
  }
 kappa <- kappa/k
 B <- k*(kappa-tau)*B1/(n*(1-2*tau+kappa*tau))
 A <- t(A1)%*%A2
 SB <- s2*(((1-atilde)^2)%*%((t(Xhat)%*%Z)%*%Mzz%*%(t(Z)%*%Xhat))+(atilde^2)%*%((t(Xhat)%*%Xhat)-(t(Xhat)%*%Z)%*%Mzz%*%(t(Z)%*%Xhat)))
 VCc = solve(Mxzzx - atilde%*%Mxx)%*%(SB+A+t(A)+B)%*%solve(Mxzzx - atilde%*%Mxx)

 return(list(b=b,VC=VC, VCb=VCb, VCc=VCc))
}