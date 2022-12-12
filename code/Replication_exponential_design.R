rm(list=ls())
# setwd("") # setting working directory
# libraries
library(mvtnorm)
library(lassoshooting)
library(xlsx)
library(glmnet)
library(MASS)
source('LassoShooting2.R')
source('fuller.R')
source('tsls.R')
set.seed(12345)
# parameter setting
nRep <- 500
n <- c(100,250)
nn <- length(n)
p <- 100
pnz <- 5
indp <- 0:(p-1)
Fstat <- c(30,180)
nF <- length(Fstat)
s2e <- 1
Cev <- .6
s2z <- 1
szz <- .5
pi1  <- .7
alpha <- 1
K <- 15

# Common design elements
SZ <- s2z*toeplitz((szz)^indp)
cSZ <- chol(SZ)
cFS <- pi1^indp
scale <- matrix(0, nrow=nF, ncol=nn)
s2v <- matrix(0, nrow=nF, ncol=nn)
for (ii in 1:nn) {
  for (jj in 1:nF) {
  scale[jj,ii] <- sqrt(Fstat[jj]/((Fstat[jj]+n[ii])*t(cFS)%*%SZ%*%cFS))
  s2v[jj,ii] <- 1-(scale[jj,ii]^2)*t(cFS)%*%SZ%*%cFS
  }
}

sev <- Cev*sqrt(s2e)*sqrt(s2v)


# Initialize matrices for conventional estimators
b2sls <- array(0, dim=c(nRep,nn,nF))
bfull <- array(0, dim=c(nRep,nn,nF))
s2sls <- array(0, dim=c(nRep,nn,nF))
sfull <- array(0, dim=c(nRep,nn,nF))
# Initialize matrices for IV with highest corr Z
b2sls1 <- array(0, dim=c(nRep,nn,nF))
s2sls1 <- array(0, dim=c(nRep,nn,nF))
# Initialize matrices for IV-LASSO estimators
FS <- array(0, dim=c(nRep,nn,nF))
blassoC <- array(0, dim=c(nRep,nn,nF))
blassoCF <- array(0, dim=c(nRep,nn,nF))
slassoC <- array(0, dim=c(nRep,nn,nF))
slassoCF <- array(0, dim=c(nRep,nn,nF))
FSV <- array(0, dim=c(nRep,nn,nF))
blassoCV <- array(0, dim=c(nRep,nn,nF))
blassoCFV <- array(0, dim=c(nRep,nn,nF))
slassoCV <- array(0, dim=c(nRep,nn,nF))
slassoCFV <- array(0, dim=c(nRep,nn,nF))
FSX <- array(0, dim=c(nRep,nn,nF))
blassoCX <- array(0, dim=c(nRep,nn,nF))
blassoCFX <- array(0, dim=c(nRep,nn,nF))
slassoCX <- array(0, dim=c(nRep,nn,nF))
slassoCFX <- array(0, dim=c(nRep,nn,nF))
FSn <- array(0, dim=c(nRep,nn,nF))
blassoCn <- array(0, dim=c(nRep,nn,nF))
blassoCFn <- array(0, dim=c(nRep,nn,nF))
slassoCn <- array(0, dim=c(nRep,nn,nF))
slassoCFn <- array(0, dim=c(nRep,nn,nF))
# Initialising for sup score test
aTest <- seq(from= .5, to=1.5, by=0.01);
na <- length(aTest);
supScore <- array(0, dim=c(nRep,na,nn,nF))
supScore05 <- array(0, dim=c(nRep,na,nn,nF))
lambdaSS05 <- qnorm(1-.05/(2*p))
# Initialize matrix for CV Ridge penalty value
LambdaRidgeA  <- array(0, dim=c(nRep,nn,nF))
LambdaRidgeB <-  array(0, dim=c(nRep,nn,nF))
# Initialize matrices for IV with highest corr Z including ridge
Rb2sls1A <- array(0, dim=c(nRep,nn,nF))
Rs2sls1A <- array(0, dim=c(nRep,nn,nF))
Rb2sls1B <- array(0, dim=c(nRep,nn,nF))
Rs2sls1B <- array(0, dim=c(nRep,nn,nF))
# Initialize matrices for IV-LASSO estimators including ridge
RFSA <- array(0, dim=c(nRep,nn,nF))
RblassoCA <- array(0, dim=c(nRep,nn,nF))
RblassoCFA <- array(0, dim=c(nRep,nn,nF))
RslassoCA <- array(0, dim=c(nRep,nn,nF))
RslassoCFA <- array(0, dim=c(nRep,nn,nF))
RFSB <-  array(0, dim=c(nRep,nn,nF))
RblassoCB <- array(0, dim=c(nRep,nn,nF))
RblassoCFB <-  array(0, dim=c(nRep,nn,nF))
RslassoCB <-  array(0, dim=c(nRep,nn,nF))
RslassoCFB <-  array(0, dim=c(nRep,nn,nF))
RFSC <-  array(0, dim=c(nRep,nn,nF))
RblassoCC <-  array(0, dim=c(nRep,nn,nF))
RblassoCFC <-  array(0, dim=c(nRep,nn,nF))
RslassoCC <-  array(0, dim=c(nRep,nn,nF))
RslassoCFC <-  array(0, dim=c(nRep,nn,nF))
# Initialize matrices for sup-score test inference including ridge
RaTestA <- seq(.5,1.5, by=0.01)
RnaA <- length(RaTestA)
RsupScoreA <-  array(0, dim=c(nRep,RnaA,nn,nF))
RsupScore05A  <-  array(0, dim=c(nRep,RnaA,nn,nF))
RaTestB <- seq(0.5,1.5, by=0.01)
RnaB <- length(RaTestB)
RsupScoreB <-  array(0, dim=c(nRep,RnaA, nn,nF))
RsupScore05B <-  array(0, dim=c(nRep,RnaA, nn,nF))
RlambdaSS05 <- qnorm(1-.05/(2*(p+1)))


for (ii in 1:nn) {
  for (jj in 1:nF){
    nUse <- n[ii]
    SU <- matrix(c(s2e,sev[jj,ii],sev[jj,ii],s2v[jj,ii]), ncol=2)
    cSU <- chol(SU)
    lambda0C <- 2.2*sqrt(2*log(2*p*log(nUse)/.1))
    Rlambda0C <- 2.2*sqrt(2*log(2*(p+1)*log(nUse/2)/.1))
    #  Calculate all the estimators, etc. inside this loop
    for (kk in 1:nRep){
      if (floor((kk-1)/10)==(kk-1)/10) {print(kk)}
      # Data Generating
      zorig <- matrix(rnorm(nUse*p), nrow=nUse,ncol=p)%*%cSZ
      U <- matrix(rnorm(nUse*2), ncol=2)%*%cSU
      xorig <- scale[jj,ii]*zorig%*%cFS+U[,2]
      yorig <- alpha*xorig+U[,1]
      Z <- zorig - rep(1,nUse)*mean(zorig)
      X <- xorig - mean(xorig)
      Y <- yorig - mean(yorig)
      
      #2SLS and Fuller
      
      if (nUse <=p) {
        pUse <- sample(p,nUse-2, replace=FALSE, prob=rep(1/p,p))
        ZUSE <- Z[,pUse]
      } else {
        ZUSE <- Z
      }
      twosls <- tsls(Y,X,x=NULL,ZUSE)
      #[btemp1,VCtemp1] = tsls(y,x,[],ZUSE);
      ful <- fuller(Y,X,x=NULL,ZUSE,1)
      #[btemp2,~,~,VCtemp2] = fuller(y,x,[],ZUSE,1);
      
      b2sls[kk,ii,jj] = twosls$b
      s2sls[kk,ii,jj] = sqrt(twosls$VC1)
      bfull[kk,ii,jj] = ful$b
      sfull[kk,ii,jj] = sqrt(ful$VCc)
      
      lambda0C = 2.2*sqrt(2*log(2*p*log(nUse)/.1))
      Rlambda0C = 2.2*sqrt(2*log(2*(p+1)*log(nUse/2)/.1))
      
      # Highest correlation 2SLS
      Z1 <- Z[,which.max(cor(Z,X))]
      
      twosls <- tsls(Y,X,x=NULL,Z1)
      
      b2sls1[kk,ii,jj] <- twosls$b
      s2sls1[kk,ii,jj] <- sqrt(twosls$VC1)
      
      # LASSO estimators
      e0 <- X
      Ups0 <- sqrt(t(t(e0^2)%*%(Z^2)))
      lambda <- lambda0C*Ups0 # choice paper    
      coefTemp <- LassoShooting2(Z,X, lambda)
      ind0 <- (abs(coefTemp$w) > 0)
      Z0 <- as.matrix(Z[,ind0])
      Z1 <- Z0
      ind1 <- ind0
             for (mm in 1:K) {
               if (dim(Z0)[2]==0) {break}
               e1 <- X-Z0%*%coef(lm(X ~ -1 + Z0))
               Ups1 <- sqrt(t(t(e1^2)%*%(Z^2)))
               coefTemp <- LassoShooting2(Z,X,lambda0C*Ups1)$w
               #print(Ups1)
               ind1 <- (abs(coefTemp) > 0)
               Z0 <- as.matrix(Z[,ind1])
             }
            if (dim(Z0)[2]!=0) {Z1 <- as.matrix(Z[,ind1])}
      
      if (dim(Z0)[2]==0) {
        blassoC[kk,ii,jj] <- NaN
        slassoC[kk,ii,jj] <- NaN
        blassoCF[kk,ii,jj] <- NaN
        slassoCF[kk,ii,jj] <- NaN
        FS[kk,ii,jj] <- 0 } else {
          bfs <- lm(X ~ -1 + Z1)$coef
          efs <- X - Z1%*%as.matrix(bfs)
          Vfs <- solve(t(Z1)%*%Z1)%*%(t(Z1*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%Z1)%*%solve(t(Z1)%*%Z1)
          FS[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
          twosls <- tsls(Y,X,x=NULL,Z1)
          ful <- fuller(Y,X,x=NULL,Z1,1)
          blassoC[kk,ii,jj] <-  twosls$b
          slassoC[kk,ii,jj] <- sqrt(twosls$VC1)
          blassoCF[kk,ii,jj] <- ful$b
          slassoCF[kk,ii,jj] <-  sqrt(ful$VCc)      
        }
      
      # X dependent
      R <- 500 # number of simulations
      sim <- vector("numeric", length=R)
      for (l in 1:R) {
        g <- matrix(rep(rnorm(nUse), each=p), ncol=p, byrow=TRUE)
        sim[l] <- nUse*max(2*colMeans(Z*g))
      }
      sigma0 <- sqrt(var(X))
      lambda0 <- rep(1.1*quantile(sim, probs=1-0.05)*sigma0,p)
      coefTemp <- LassoShooting2(Z,X, lambda0)
      lambda <- lambda0 # !!
      
#       sigma1 <- sqrt(mean((X - Z%*%coefTemp$w)^2))
#       sim <- vector("numeric", length=R)
#       for (l in 1:R) {
#         g <- matrix(rep(rnorm(nUse), each=p), ncol=p, byrow=TRUE)
#         sim[l] <- nUse*max(2*colMeans(Z*g))
#       }
#       lambda <- rep(1.1*quantile(sim, probs=1-0.05)*sqrt(sigma1),p)
         
      coefTemp <- LassoShooting2(Z,X, lambda)
      ind0 <- (abs(coefTemp$w) > 0)
      Z0 <- as.matrix(Z[,ind0])
      Z1 <- Z0
      ind1 <- ind0
      if (dim(Z0)[2]==0) {
        blassoCn[kk,ii,jj] <- NaN
        slassoCn[kk,ii,jj] <- NaN
        blassoCFn[kk,ii,jj] <- NaN
        slassoCFn[kk,ii,jj] <- NaN
        FSn[kk,ii,jj] <- 0 } else {
          bfs <- lm(X ~ -1 + Z1)$coef
          efs <- X - Z1%*%as.matrix(bfs)
          Vfs <- solve(t(Z1)%*%Z1)%*%(t(Z1*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%Z1)%*%solve(t(Z1)%*%Z1)
          FSn[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
          twosls <- tsls(Y,X,x=NULL,Z1)
          ful <- fuller(Y,X,x=NULL,Z1,1)
          blassoCn[kk,ii,jj] <-  twosls$b
          slassoCn[kk,ii,jj] <- sqrt(twosls$VC1)
          blassoCFn[kk,ii,jj] <- ful$b
          slassoCFn[kk,ii,jj] <-  sqrt(ful$VCc)      
        }
      
      # CV choice
      cv <- cv.glmnet(x=Z,y=X, alpha=1, lambda=NULL)
      lambda <- rep(cv$lambda.min*nUse*2,p)      
      coefTemp <- LassoShooting2(Z,X, lambda)
      ind0 <- (abs(coefTemp$w) > 0)
      Z0 <- as.matrix(Z[,ind0])
      Z1 <- Z0
      ind1 <- ind0
      if (dim(Z0)[2]==0) {
        blassoCV[kk,ii,jj] <- NaN
        slassoCV[kk,ii,jj] <- NaN
        blassoCFV[kk,ii,jj] <- NaN
        slassoCFV[kk,ii,jj] <- NaN
        FSV[kk,ii,jj] <- 0 } else {
          bfs <- lm(X ~ -1 + Z1)$coef
          efs <- X - Z1%*%as.matrix(bfs)
          Vfs <- solve(t(Z1)%*%Z1)%*%(t(Z1*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%Z1)%*%solve(t(Z1)%*%Z1)
          FSV[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
          twosls <- tsls(Y,X,x=NULL,Z1)
          ful <- fuller(Y,X,x=NULL,Z1,1)
          blassoCV[kk,ii,jj] <-  twosls$b
          slassoCV[kk,ii,jj] <- sqrt(twosls$VC1)
          blassoCFV[kk,ii,jj] <- ful$b
          slassoCFV[kk,ii,jj] <-  sqrt(ful$VCc)      
        } 
      
      # X independent choice
      lambda <- rep(2*1.1*sqrt(nUse)*qnorm(1-0.1/(2*p))*sqrt(var(X)),p)   
      coefTemp <- LassoShooting2(Z,X, lambda)
      ind0 <- (abs(coefTemp$w) > 0)
      Z0 <- as.matrix(Z[,ind0])
      ind1 <- ind0
      Z1 <- Z0
      if (dim(Z0)[2]==0) {
        blassoCX[kk,ii,jj] <- NaN
        slassoCX[kk,ii,jj] <- NaN
        blassoCFX[kk,ii,jj] <- NaN
        slassoCFX[kk,ii,jj] <- NaN
        FSX[kk,ii,jj] <- 0 } else {
          bfs <- lm(X ~ -1 + Z1)$coef
          efs <- X - Z1%*%as.matrix(bfs)
          Vfs <- solve(t(Z1)%*%Z1)%*%(t(Z1*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%Z1)%*%solve(t(Z1)%*%Z1)
          FSX[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
          twosls <- tsls(Y,X,x=NULL,Z1)
          ful <- fuller(Y,X,x=NULL,Z1,1)
          blassoCX[kk,ii,jj] <-  twosls$b
          slassoCX[kk,ii,jj] <- sqrt(twosls$VC1)
          blassoCFX[kk,ii,jj] <- ful$b
          slassoCFX[kk,ii,jj] <-  sqrt(ful$VCc)      
        }    
      
      # Sup-Score Tests
      for (mm in 1:na) {
      aEval <- aTest[mm]
      eTemp <- Y-aEval*X
      ScoreVec <- t(eTemp)%*%Z
      ScoreStd <- sqrt(t(eTemp^2)%*%(Z^2))
      ScaledScore <- ScoreVec/(1.1*ScoreStd)
      supScore05[kk,mm,ii,jj] <- max(abs(ScaledScore)) < lambdaSS05
      supScore[kk,mm,ii,jj] <- max(abs(ScaledScore))
      }
      
      
      
      # Calculate estimators putting in ridge fit with coefficients 
      # estimated from 1/2 of sample
      
      #Split sample
      index <- sample(1:nUse,nUse/2, replace=F)
      UseRidgeA <- index
      nRidgeA <- length(UseRidgeA)
      UseRidgeB <- setdiff(1:nUse,index)
      nRidgeB <- length(UseRidgeB)
           
      yA <- yorig[UseRidgeA]
      xA <- xorig[UseRidgeA]
      ZA <- zorig[UseRidgeA,]
      yB <- yorig[UseRidgeB]
      xB <- xorig[UseRidgeB]
      ZB <- zorig[UseRidgeB,]
      
      
      cv.ridgeA <- cv.glmnet(x=ZA,y=xA, family="gaussian", alpha=0, nfolds=nUse )
      LambdaRidgeA[kk,ii,jj] <- cv.ridgeA$lambda.min 
      
      cv.ridgeB <- cv.glmnet(x=ZB,y=xB, family="gaussian", alpha=0, nfolds=nUse )
      LambdaRidgeA[kk,ii,jj] <- cv.ridgeB$lambda.min
      
      #out <- glmnet(x=X,y=Y,alpha =1, lambda = bestlam, intercept=FALSE, standardize=FALSE)
      #lasso.coef <- as.vector(predict(out, type ="coefficients", s=bestlam))[-1]
      
      RidgeFitB <- predict(cv.ridgeA, newx= ZB, type ="response", s= LambdaRidgeA[kk,ii,jj])
      RidgeFitA <- predict(cv.ridgeB, newx= ZA, type ="response", s= LambdaRidgeB[kk,ii,jj])
      
      # to avoid multicollinearity
      if (sum(RidgeFitB==RidgeFitB[1])==length(RidgeFitB)) {RidgeFitB <- rnorm(length(RidgeFitB))}
      
      ZB <- cbind(ZB, RidgeFitB)
      ZA <- cbind(ZA, RidgeFitA)
       
      ZA <- ZA - rep(1,nRidgeA)%*%t(as.matrix(colMeans(ZA)))
      xA <- xA - mean(xA)
      yA <- yA - mean(yA)
    
      ZB <- ZB - rep(1,nRidgeB)%*%t(as.matrix(colMeans(ZB)))
      xB <- xB - mean(xB)
      yB <- yB - mean(yB)
      
      # Highest correlation 2SLS
      ZL1A <- ZA[,which.max(cor(ZA,xA))]
      twosls <- tsls(yA,xA,x=NULL,ZL1A)
      Rb2sls1A[kk,ii,jj] <- twosls$b
      Rs2sls1A[kk,ii,jj] <- sqrt(twosls$VC1)
      
      ZL1B <- ZB[,which.max(cor(ZB,xB))]
      twosls <- tsls(yB,xB,x=NULL,ZL1B)
      Rb2sls1B[kk,ii,jj] <- twosls$b
      Rs2sls1B[kk,ii,jj] <- sqrt(twosls$VC1)
   
                  
      # LASSO estimators
      # A Sample
      e0 <- xA
      Ups0 <- sqrt(t(t(e0^2)%*%(ZA^2)))
      lambda <- lambda0C*Ups0
      coefTemp <- LassoShooting2(ZA,xA,lambda)
      ind0 <- (abs(coefTemp$w) > 0)
      ZL0 <-as.matrix(ZA[,ind0])
      ZL1 <- ZL0
      ind1 <- ind0
      for (mm in 1:15) {
        if (dim(ZL0)[2]==0) {break}
        e1 <- xA-ZL0%*%coef(lm(xA ~ -1 + ZL0))
        Ups1 <- sqrt(t(t(e1^2)%*%(ZA^2)))
        coefTemp <- LassoShooting2(ZA,xA,Rlambda0C*Ups1)$w
        ind1 <- (abs(coefTemp) > 0)
        ZL0 <- as.matrix(ZA[,ind1])
    }
    if (dim(ZL0)[2]!=0) {ZL1A <- as.matrix(ZA[,ind1])}
                              
       if (dim(ZL0)[2]==0) {
          RblassoCA[kk,ii,jj] <- NaN
          RslassoCA[kk,ii,jj] <- NaN
          RblassoCFA[kk,ii,jj] <- NaN
          RslassoCFA[kk,ii,jj] <- NaN
          RFSA[kk,ii,jj] <- 0
          WA <- 0} else {
            bfs <- lm(xA ~ -1 + ZL1A)$coef
            efs <- xA - ZL1A%*%as.matrix(bfs)
            Vfs <- solve(t(ZL1A)%*%ZL1A)%*%(t(ZL1A*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%ZL1A)%*%solve(t(ZL1A)%*%ZL1A)
            RFSA[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
            twosls <- tsls(yA,xA,x=NULL,ZL1A)
            ful <- fuller(yA,xA,x=NULL,ZL1A,1)
            RblassoCA[kk,ii,jj] <- twosls$b
            RslassoCA[kk,ii,jj] <- sqrt(twosls$VC1)
            RblassoCFA[kk,ii,jj] <- ful$b
            RslassoCFA[kk,ii,jj] <-  sqrt(ful$VCc)
            WA <- RFSA[kk,ii,jj]*sum(ind1)
          }
          
                         
       # B Sample
       e0 <- xB
       Ups0 <- sqrt(t(t(e0^2)%*%(ZB^2)))
       lambda <- lambda0C*Ups0
       coefTemp <- LassoShooting2(ZB,xB,lambda) # here problem!
       ind0 <- (abs(coefTemp$w) > 0)
       ZL0 <-as.matrix(ZB[,ind0])
       ZL1 <- ZL0
       ind1 <- ind0
       for (mm in 1:15) {
            if (dim(ZL0)[2]==0) {break}
            e1 <- xB-ZL0%*%coef(lm(xB ~ -1 + ZL0))
            Ups1 <- sqrt(t(t(e1^2)%*%(ZB^2)))
            coefTemp <- LassoShooting2(ZB,xB,Rlambda0C*Ups1)$w
            ind1 <- (abs(coefTemp) > 0)
            ZL0 <- as.matrix(ZB[,ind1])
            }
       if (dim(ZL0)[2]!=0) {ZL1B <- as.matrix(ZB[,ind1])}
                                              
       if (dim(ZL0)[2]==0) {
             RblassoCB[kk,ii,jj] <- NaN
             RslassoCB[kk,ii,jj] <- NaN
             RblassoCFB[kk,ii,jj] <- NaN
             RslassoCFB[kk,ii,jj] <- NaN
             RFSB[kk,ii,jj] <- 0
             WB <- 0} else {
             bfs <- lm(xB ~ -1 + ZL1B)$coef
             efs <- xB - ZL1B%*%as.matrix(bfs)
             Vfs <- solve(t(ZL1B)%*%ZL1B)%*%(t(ZL1B*(efs^2%*%matrix(rep(1,sum(ind1)), nrow=1)))%*%ZL1B)%*%solve(t(ZL1B)%*%ZL1B)
             RFSB[kk,ii,jj] = t(bfs)%*%solve(Vfs)%*%bfs/sum(ind1)
             twosls <- tsls(yB,xB,x=NULL,ZL1B)
             ful <- fuller(yB,xB,x=NULL,ZL1B,1)
             RblassoCB[kk,ii,jj] <- twosls$b
             RslassoCB[kk,ii,jj] <- sqrt(twosls$VC1)
             RblassoCFB[kk,ii,jj] <- ful$b
             RslassoCFB[kk,ii,jj] <-  sqrt(ful$VCc)
             WB <- RFSB[kk,ii,jj]*sum(ind1);
             }
                                              

              l <- dim(ZL1A)[2]
              k <- dim(ZL1B)[2]
              if (is.vector(ZL1A)) {l <- 1}
              if (is.vector(ZL1B)) {k <- 1}
            # Combine estimators
             if (l==0 && k==0) {
              RblassoCC[kk,ii,jj] <- NaN
              RslassoCC[kk,ii,jj] <- NaN
              RblassoCFC[kk,ii,jj] <- NaN
              RslassoCFC[kk,ii,jj] <- NaN
              RFSC[kk,ii,jj] <-  0 }
              if (l==0 && k!=0) {
              RblassoCC[kk,ii,jj] <- RblassoCB[kk,ii,jj]
              RslassoCC[kk,ii,jj] <- RslassoCB[kk,ii,jj]
              RblassoCFC[kk,ii,jj] <- RblassoCFB[kk,ii,jj]
              RslassoCFC[kk,ii,jj] <- RslassoCFB[kk,ii,jj]
              RFSC[kk,ii,jj] <- RFSB[kk,ii,jj]}
             if (l!=0 && k==0) {            
              RblassoCC[kk,ii,jj] <- RblassoCA[kk,ii,jj]
              RslassoCC[kk,ii,jj] <- RslassoCA[kk,ii,jj]
              RblassoCFC[kk,ii,jj] <- RblassoCFA[kk,ii,jj]
              RslassoCFC[kk,ii,jj] <- RslassoCFA[kk,ii,jj]
              RFSC[kk,ii,jj] <- RFSA[kk,ii,jj]}
            
             if (l!=0 && k!=0) {              
              weightA <- RslassoCB[kk,ii,jj]^2/(RslassoCA[kk,ii,jj]^2 + RslassoCB[kk,ii,jj]^2)
              weightB <- RslassoCA[kk,ii,jj]^2/(RslassoCA[kk,ii,jj]^2 + RslassoCB[kk,ii,jj]^2)
              RblassoCC[kk,ii,jj] <- weightA*RblassoCA[kk,ii,jj] + weightB*RblassoCB[kk,ii,jj]
              RslassoCC[kk,ii,jj] <- sqrt((weightA*RslassoCA[kk,ii,jj])^2 + (weightB*RslassoCB[kk,ii,jj])^2)
              weightA <- RslassoCFB[kk,ii,jj]^2/(RslassoCFA[kk,ii,jj]^2 + RslassoCFB[kk,ii,jj]^2)
              weightB <- RslassoCFA[kk,ii,jj]^2/(RslassoCFA[kk,ii,jj]^2 + RslassoCFB[kk,ii,jj]^2)
              RblassoCFC[kk,ii,jj] <- weightA*RblassoCFA[kk,ii,jj] + weightB*RblassoCFB[kk,ii,jj]
              RslassoCFC[kk,ii,jj] <- sqrt((weightA*RslassoCFA[kk,ii,jj])^2 + (weightB*RslassoCFB[kk,ii,jj])^2)
              RFSC[kk,ii,jj] <- (WA+WB)/(l+k)
             }
                                   
                  
             # Sup-Score Tests
             for (mm in 1:RnaA) {
                aEval <- RaTestA[mm]
                eTemp <- yA-aEval*xA
                ScoreVec <- t(eTemp)%*%ZA
                ScoreStd <- sqrt(t(eTemp^2)%*%(ZA^2))
                ScaledScore <- ScoreVec/(1.1*ScoreStd)
                RsupScore05A[kk,mm,ii,jj] <- max(abs(ScaledScore)) < RlambdaSS05 #>?
                RsupScoreA[kk,mm,ii,jj] <- max(abs(ScaledScore))
                eTemp <- yB-aEval*xB
                ScoreVec <- t(eTemp)%*%ZB
                ScoreStd <- sqrt(t(eTemp^2)%*%(ZB^2))
                ScaledScore <- ScoreVec/(1.1*ScoreStd)
                RsupScore05B[kk,mm,ii,jj] <- max(abs(ScaledScore)) < RlambdaSS05 #>?
                RsupScoreB[kk,mm,ii,jj] <- max(abs(ScaledScore))
             }
      
      
      }
    }
  }
  
# Preparation results

Results <- matrix(NA, ncol=4, nrow=12)
colnames(Results) <- c("N(0)", "Bias", "MAD", "rp(0.05)")
rownames(Results) <- c("2SLS(100)", "FULL(100)", "Post-LASSO (paper)", "Post-LASSO-F (paper)",  
                       "Post-LASSO (X-dep.)", "Post-LASSO-F (X-dep.)",
                       "Post-LASSO (X-indep.)", "Post-LASSO-F (X-indep.)", "Post-Lasso (Ridge)",  "Post-Lasso-F (Ridge)", "Post-LASSO (CV)", "Post-LASSO-F (CV)")
# n=100, mu=30 [1,1]
#sfull[is.nan(sfull)] <- 10000
Results[1,1] <- sum(is.nan(b2sls[,1,1]))
Results[1,2] <- mean(b2sls[,1,1])-1
Results[1,3] <- mean(abs(b2sls[,1,1]-1))
Results[1,4] <- (500 - sum(abs((b2sls[,1,1]-1)/s2sls[,1,1]) < qnorm(0.95)))/500

Results[2,1] <- sum(is.nan(bfull[,1,1]))
Results[2,2] <- mean(bfull[,1,1])-1
Results[2,3] <- mean(abs(bfull[,1,1]-1))
Results[2,4] <- (500 - sum(abs((bfull[,1,1]-1)/sfull[,1,1]) < qnorm(0.95)))/500

Results[3,1] <- sum(is.nan(blassoC[,1,1]))
blassoC_aux <- blassoC
blassoC_aux[is.nan(blassoC_aux)] <- b2sls1[is.nan(blassoC_aux)]
Results[3,2] <- mean(blassoC_aux[,1,1])-1
Results[3,3] <- mean(abs(blassoC_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoC[,1,1]))
test1 <- sum(abs((blassoC[element_chosen,1,1]-1)/slassoC[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[3,4] <- (500 - test1 - test2)/500

Results[4,1] <- sum(is.nan(blassoCF[,1,1]))
blassoCF_aux <- blassoCF
blassoCF_aux[is.nan(blassoCF_aux)] <- b2sls1[is.nan(blassoCF_aux)]
Results[4,2] <- mean(blassoCF_aux[,1,1])-1
Results[4,3] <- mean(abs(blassoCF_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCF[,1,1]))
test1 <- sum(abs((blassoCF[element_chosen,1,1]-1)/slassoCF[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[4,4] <- (500 - test1 - test2)/500


Results[5,1] <- sum(is.nan(blassoCn[,1,1]))
blassoCn_aux <- blassoCn
blassoCn_aux[is.nan(blassoCn_aux)] <- b2sls1[is.nan(blassoCn_aux)]
Results[5,2] <- mean(blassoCn_aux[,1,1])-1
Results[5,3] <- mean(abs(blassoCn_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCn[,1,1]))
test1 <- sum(abs((blassoCn[element_chosen,1,1]-1)/slassoCn[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[5,4] <- (500 - test1 - test2)/500

Results[6,1] <- sum(is.nan(blassoCFn[,1,1]))
blassoCFn_aux <- blassoCFn
blassoCFn_aux[is.nan(blassoCFn_aux)] <- b2sls1[is.nan(blassoCFn_aux)]
Results[6,2] <- mean(blassoCFn_aux[,1,1])-1
Results[6,3] <- mean(abs(blassoCFn_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCFn[,1,1]))
test1 <- sum(abs((blassoCFn[element_chosen,1,1]-1)/slassoCFn[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[6,4] <- (500 - test1 - test2)/500


Results[7,1] <- sum(is.nan(blassoCX[,1,1]))
blassoCX_aux <- blassoCX
blassoCX_aux[is.nan(blassoCX_aux)] <- b2sls1[is.nan(blassoCX_aux)]
Results[7,2] <- mean(blassoCX_aux[,1,1])-1
Results[7,3] <- mean(abs(blassoCX_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCX[,1,1]))
test1 <- sum(abs((blassoCX[element_chosen,1,1]-1)/slassoCX[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[7,4] <- (500 - test1 - test2)/500

Results[8,1] <- sum(is.nan(blassoCFX[,1,1]))
blassoCFX_aux <- blassoCFX
blassoCFX_aux[is.nan(blassoCFX_aux)] <- b2sls1[is.nan(blassoCFX_aux)]
Results[8,2] <- mean(blassoCFX_aux[,1,1])-1
Results[8,3] <- mean(abs(blassoCFX_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCFX[,1,1]))
test1 <- sum(abs((blassoCFX[element_chosen,1,1]-1)/slassoCFX[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[8,4] <- (500 - test1 - test2)/500
      
Results[9,1] <- sum(is.nan(RblassoCC[,1,1]))
RblassoCC_aux <- RblassoCC
RblassoCC_aux[is.nan(RblassoCC_aux)] <- b2sls1[is.nan(RblassoCC_aux)]
Results[9,2] <- mean(RblassoCC_aux[,1,1])-1
Results[9,3] <- mean(abs(RblassoCC_aux[,1,1]-1))
element_chosen <- which(!is.nan(RblassoCC[,1,1]))
test1 <- sum(abs((RblassoCC[element_chosen,1,1]-1)/RslassoCC[element_chosen,1,1]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,1,1])+sum(RsupScore05B[-element_chosen,51,1,1]))
if (sum(element_chosen)==0) {
  test1 <- 0
  test2 <-  1/2*(sum(RsupScore05A[,51,1,1])+sum(RsupScore05B[,51,1,1]))
}
Results[9,4] <- (500 - test1 - test2)/500
      
Results[10,1] <- sum(is.nan(RblassoCFC[,1,1]))
RblassoCFC_aux <- RblassoCFC
RblassoCFC_aux[is.nan(RblassoCFC_aux)] <- b2sls1[is.nan(RblassoCFC_aux)]
Results[10,2] <- mean(RblassoCFC_aux[,1,1])-1
Results[10,3] <- mean(abs(RblassoCFC_aux[,1,1]-1))
element_chosen <- which(!is.nan(RblassoCFC[,1,1]))
test1 <- sum(abs((RblassoCFC[element_chosen,1,1]-1)/RslassoCFC[element_chosen,1,1]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,1,1])+sum(RsupScore05B[-element_chosen,51,1,1]))
if (sum(element_chosen)==0) {
  test1 <- 0
  test2 <-  1/2*(sum(RsupScore05A[,51,1,1])+sum(RsupScore05B[,51,1,1]))
}
Results[10,4] <- (500 - test1 - test2)/500

Results[11,1] <- sum(is.nan(blassoCV[,1,1]))
blassoCV_aux <- blassoCV
blassoCV_aux[is.nan(blassoCV_aux)] <- b2sls1[is.nan(blassoCV_aux)]
Results[11,2] <- mean(blassoCV_aux[,1,1])-1
Results[11,3] <- mean(abs(blassoCV_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCV[,1,1]))
test1 <- sum(abs((blassoCV[element_chosen,1,1]-1)/slassoCV[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[11,4] <- (500 - test1 - test2)/500

Results[12,1] <- sum(is.nan(blassoCFV[,1,1]))
blassoCFV_aux <- blassoCFV
blassoCFV_aux[is.nan(blassoCFV_aux)] <- b2sls1[is.nan(blassoCFV_aux)]
Results[12,2] <- mean(blassoCFV_aux[,1,1])-1
Results[12,3] <- mean(abs(blassoCFV_aux[,1,1]-1))
element_chosen <- which(!is.nan(blassoCFV[,1,1]))
test1 <- sum(abs((blassoCFV[element_chosen,1,1]-1)/slassoCFV[element_chosen,1,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,1])
Results[12,4] <- (500 - test1 - test2)/500


Results_100_30 <- Results

# n=250, mu=30 [2,1]
sfull[is.nan(sfull)] <- 10000
Results[1,1] <- sum(is.nan(b2sls[,2,1]))
Results[1,2] <- mean(b2sls[,2,1])-1
Results[1,3] <- mean(abs(b2sls[,2,1]-1))
Results[1,4] <- (500 - sum(abs((b2sls[,2,1]-1)/s2sls[,2,1]) < qnorm(0.95)))/500

Results[2,1] <- sum(is.nan(bfull[,2,1]))
Results[2,2] <- mean(bfull[,2,1])-1
Results[2,3] <- mean(abs(bfull[,2,1]-1))
Results[2,4] <- (500 - sum(abs((bfull[,2,1]-1)/sfull[,2,1]) < qnorm(0.95)))/500

Results[3,1] <- sum(is.nan(blassoC[,2,1]))
blassoC_aux <- blassoC
blassoC_aux[is.nan(blassoC_aux)] <- b2sls1[is.nan(blassoC_aux)]
Results[3,2] <- mean(blassoC_aux[,2,1])-1
Results[3,3] <- mean(abs(blassoC_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoC[,2,1]))
test1 <- sum(abs((blassoC[element_chosen,2,1]-1)/slassoC[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[3,4] <- (500 - test1 - test2)/500

Results[4,1] <- sum(is.nan(blassoCF[,2,1]))
blassoCF_aux <- blassoCF
blassoCF_aux[is.nan(blassoCF_aux)] <- b2sls1[is.nan(blassoCF_aux)]
Results[4,2] <- mean(blassoCF_aux[,2,1])-1
Results[4,3] <- mean(abs(blassoCF_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCF[,2,1]))
test1 <- sum(abs((blassoCF[element_chosen,2,1]-1)/slassoCF[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[4,4] <- (500 - test1 - test2)/500


Results[5,1] <- sum(is.nan(blassoCn[,2,1]))
blassoCn_aux <- blassoCn
blassoCn_aux[is.nan(blassoCn_aux)] <- b2sls1[is.nan(blassoCn_aux)]
Results[5,2] <- mean(blassoCn_aux[,2,1])-1
Results[5,3] <- mean(abs(blassoCn_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCn[,2,1]))
test1 <- sum(abs((blassoCn[element_chosen,2,1]-1)/slassoCn[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[5,4] <- (500 - test1 - test2)/500

Results[6,1] <- sum(is.nan(blassoCFn[,2,1]))
blassoCFn_aux <- blassoCFn
blassoCFn_aux[is.nan(blassoCFn_aux)] <- b2sls1[is.nan(blassoCFn_aux)]
Results[6,2] <- mean(blassoCFn_aux[,2,1])-1
Results[6,3] <- mean(abs(blassoCFn_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCFn[,2,1]))
test1 <- sum(abs((blassoCFn[element_chosen,2,1]-1)/slassoCFn[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[6,4] <- (500 - test1 - test2)/500

Results[7,1] <- sum(is.nan(blassoCX[,2,1]))
blassoCX_aux <- blassoCX
blassoCX_aux[is.nan(blassoCX_aux)] <- b2sls1[is.nan(blassoCX_aux)]
Results[7,2] <- mean(blassoCX_aux[,2,1])-1
Results[7,3] <- mean(abs(blassoCX_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCX[,2,1]))
test1 <- sum(abs((blassoCX[element_chosen,2,1]-1)/slassoCX[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[7,4] <- (500 - test1 - test2)/500

Results[8,1] <- sum(is.nan(blassoCFX[,2,1]))
blassoCFX_aux <- blassoCFX
blassoCFX_aux[is.nan(blassoCFX_aux)] <- b2sls1[is.nan(blassoCFX_aux)]
Results[8,2] <- mean(blassoCFX_aux[,2,1])-1
Results[8,3] <- mean(abs(blassoCFX_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCFX[,2,1]))
test1 <- sum(abs((blassoCFX[element_chosen,2,1]-1)/slassoCFX[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[8,4] <- (500 - test1 - test2)/500

Results[9,1] <- sum(is.nan(RblassoCC[,2,1]))
RblassoCC_aux <- RblassoCC
RblassoCC_aux[is.nan(RblassoCC_aux)] <- b2sls1[is.nan(RblassoCC_aux)]
Results[9,2] <- mean(RblassoCC_aux[,2,1])-1
Results[9,3] <- mean(abs(RblassoCC_aux[,2,1]-1))
element_chosen <- which(!is.nan(RblassoCC[,2,1]))
test1 <- sum(abs((RblassoCC[element_chosen,2,1]-1)/RslassoCC[element_chosen,2,1]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,2,1])+sum(RsupScore05B[-element_chosen,51,2,1]))
Results[9,4] <- (500 - test1 - test2)/500
      
Results[10,1] <- sum(is.nan(RblassoCFC[,2,1]))
RblassoCFC_aux <- RblassoCFC
RblassoCFC_aux[is.nan(RblassoCFC_aux)] <- b2sls1[is.nan(RblassoCFC_aux)]
Results[10,2] <- mean(RblassoCFC_aux[,2,1])-1
Results[10,3] <- mean(abs(RblassoCFC_aux[,2,1]-1))
element_chosen <- which(!is.nan(RblassoCFC[,2,1]))
test1 <- sum(abs((RblassoCFC[element_chosen,2,1]-1)/RslassoCFC[element_chosen,2,1]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,2,1])+sum(RsupScore05B[-element_chosen,51,2,1]))
Results[10,4] <- (500 - test1 - test2)/500

Results[11,1] <- sum(is.nan(blassoCV[,2,1]))
blassoCV_aux <- blassoCV
blassoCV_aux[is.nan(blassoCV_aux)] <- b2sls1[is.nan(blassoCV_aux)]
Results[11,2] <- mean(blassoCV_aux[,2,1])-1
Results[11,3] <- mean(abs(blassoCV_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCV[,2,1]))
test1 <- sum(abs((blassoCV[element_chosen,2,1]-1)/slassoCV[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[11,4] <- (500 - test1 - test2)/500

Results[12,1] <- sum(is.nan(blassoCFV[,2,1]))
blassoCFV_aux <- blassoCFV
blassoCFV_aux[is.nan(blassoCFV_aux)] <- b2sls1[is.nan(blassoCFV_aux)]
Results[12,2] <- mean(blassoCFV_aux[,2,1])-1
Results[12,3] <- mean(abs(blassoCFV_aux[,2,1]-1))
element_chosen <- which(!is.nan(blassoCFV[,2,1]))
test1 <- sum(abs((blassoCFV[element_chosen,2,1]-1)/slassoCFV[element_chosen,2,1]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,1])
Results[12,4] <- (500 - test1 - test2)/500

Results_250_30 <- Results


# n=100, mu=180 [1,2]
sfull[is.nan(sfull)] <- 10000
Results[1,1] <- sum(is.nan(b2sls[,1,2]))
Results[1,2] <- mean(b2sls[,1,2])-1
Results[1,3] <- mean(abs(b2sls[,1,2]-1))
Results[1,4] <- (500 - sum(abs((b2sls[,1,2]-1)/s2sls[,1,2]) < qnorm(0.95)))/500

Results[2,1] <- sum(is.nan(bfull[,1,2]))
Results[2,2] <- mean(bfull[,1,2])-1
Results[2,3] <- mean(abs(bfull[,1,2]-1))
Results[2,4] <- (500 - sum(abs((bfull[,1,2]-1)/sfull[,1,2]) < qnorm(0.95)))/500

Results[3,1] <- sum(is.nan(blassoC[,1,2]))
blassoC_aux <- blassoC
blassoC_aux[is.nan(blassoC_aux)] <- b2sls1[is.nan(blassoC_aux)]
Results[3,2] <- mean(blassoC_aux[,1,2])-1
Results[3,3] <- mean(abs(blassoC_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoC[,1,2]))
test1 <- sum(abs((blassoC[element_chosen,1,2]-1)/slassoC[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[3,4] <- (500 - test1 - test2)/500

Results[4,1] <- sum(is.nan(blassoCF[,1,2]))
blassoCF_aux <- blassoCF
blassoCF_aux[is.nan(blassoCF_aux)] <- b2sls1[is.nan(blassoCF_aux)]
Results[4,2] <- mean(blassoCF_aux[,1,2])-1
Results[4,3] <- mean(abs(blassoCF_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCF[,1,2]))
test1 <- sum(abs((blassoCF[element_chosen,1,2]-1)/slassoCF[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[4,4] <- (500 - test1 - test2)/500


Results[5,1] <- sum(is.nan(blassoCn[,1,2]))
blassoCn_aux <- blassoCn
blassoCn_aux[is.nan(blassoCn_aux)] <- b2sls1[is.nan(blassoCn_aux)]
Results[5,2] <- mean(blassoCn_aux[,1,2])-1
Results[5,3] <- mean(abs(blassoCn_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCn[,1,2]))
test1 <- sum(abs((blassoCn[element_chosen,1,2]-1)/slassoCn[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[5,4] <- (500 - test1 - test2)/500

Results[6,1] <- sum(is.nan(blassoCFn[,1,2]))
blassoCFn_aux <- blassoCFn
blassoCFn_aux[is.nan(blassoCFn_aux)] <- b2sls1[is.nan(blassoCFn_aux)]
Results[6,2] <- mean(blassoCFn_aux[,1,2])-1
Results[6,3] <- mean(abs(blassoCFn_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCFn[,1,2]))
test1 <- sum(abs((blassoCFn[element_chosen,1,2]-1)/slassoCFn[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[6,4] <- (500 - test1 - test2)/500


Results[7,1] <- sum(is.nan(blassoCX[,1,2]))
blassoCX_aux <- blassoCX
blassoCX_aux[is.nan(blassoCX_aux)] <- b2sls1[is.nan(blassoCX_aux)]
Results[7,2] <- mean(blassoCX_aux[,1,2])-1
Results[7,3] <- mean(abs(blassoCX_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCX[,1,2]))
test1 <- sum(abs((blassoCX[element_chosen,1,2]-1)/slassoCX[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[7,4] <- (500 - test1 - test2)/500

Results[8,1] <- sum(is.nan(blassoCFX[,1,2]))
blassoCFX_aux <- blassoCFX
blassoCFX_aux[is.nan(blassoCFX_aux)] <- b2sls1[is.nan(blassoCFX_aux)]
Results[8,2] <- mean(blassoCFX_aux[,1,2])-1
Results[8,3] <- mean(abs(blassoCFX_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCFX[,1,2]))
test1 <- sum(abs((blassoCFX[element_chosen,1,2]-1)/slassoCFX[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[8,4] <- (500 - test1 - test2)/500

Results[9,1] <- sum(is.nan(RblassoCC[,1,2]))
RblassoCC_aux <- RblassoCC
RblassoCC_aux[is.nan(RblassoCC_aux)] <- b2sls1[is.nan(RblassoCC_aux)]
Results[9,2] <- mean(RblassoCC_aux[,1,2])-1
Results[9,3] <- mean(abs(RblassoCC_aux[,1,2]-1))
element_chosen <- which(!is.nan(RblassoCC[,1,2]))
test1 <- sum(abs((RblassoCC[element_chosen,1,2]-1)/RslassoCC[element_chosen,2,1]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,1,2])+sum(RsupScore05B[-element_chosen,51,1,2]))
Results[9,4] <- (500 - test1 - test2)/500
      
Results[10,1] <- sum(is.nan(RblassoCFC[,1,2]))
RblassoCFC_aux <- RblassoCFC
RblassoCFC_aux[is.nan(RblassoCFC_aux)] <- b2sls1[is.nan(RblassoCFC_aux)]
Results[10,2] <- mean(RblassoCFC_aux[,1,2])-1
Results[10,3] <- mean(abs(RblassoCFC_aux[,1,2]-1))
element_chosen <- which(!is.nan(RblassoCFC[,1,2]))
test1 <- sum(abs((RblassoCFC[element_chosen,1,2]-1)/RslassoCFC[element_chosen,2,2]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,1,2])+sum(RsupScore05B[-element_chosen,51,1,2]))
Results[10,4] <- (500 - test1 - test2)/500

Results[11,1] <- sum(is.nan(blassoCV[,1,2]))
blassoCV_aux <- blassoCV
blassoCV_aux[is.nan(blassoCV_aux)] <- b2sls1[is.nan(blassoCV_aux)]
Results[11,2] <- mean(blassoCV_aux[,1,2])-1
Results[11,3] <- mean(abs(blassoCV_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCV[,1,2]))
test1 <- sum(abs((blassoCV[element_chosen,1,2]-1)/slassoCV[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[11,4] <- (500 - test1 - test2)/500

Results[12,1] <- sum(is.nan(blassoCFV[,1,2]))
blassoCFV_aux <- blassoCFV
blassoCFV_aux[is.nan(blassoCFV_aux)] <- b2sls1[is.nan(blassoCFV_aux)]
Results[12,2] <- mean(blassoCFV_aux[,1,2])-1
Results[12,3] <- mean(abs(blassoCFV_aux[,1,2]-1))
element_chosen <- which(!is.nan(blassoCFV[,1,2]))
test1 <- sum(abs((blassoCFV[element_chosen,1,2]-1)/slassoCFV[element_chosen,1,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,1,2])
Results[12,4] <- (500 - test1 - test2)/500


Results_100_180 <- Results


# n=250, mu=180 [2,2]
sfull[is.nan(sfull)] <- 10000
Results[1,1] <- sum(is.nan(b2sls[,2,2]))
Results[1,2] <- mean(b2sls[,2,2])-1
Results[1,3] <- mean(abs(b2sls[,2,2]-1))
Results[1,4] <- (500 - sum(abs((b2sls[,2,2]-1)/s2sls[,2,2]) < qnorm(0.95)))/500

Results[2,1] <- sum(is.nan(bfull[,2,2]))
Results[2,2] <- mean(bfull[,2,2])-1
Results[2,3] <- mean(abs(bfull[,2,2]-1))
Results[2,4] <- (500 - sum(abs((bfull[,2,2]-1)/sfull[,2,2]) < qnorm(0.95)))/500

Results[3,1] <- sum(is.nan(blassoC[,2,2]))
blassoC_aux <- blassoC
blassoC_aux[is.nan(blassoC_aux)] <- b2sls1[is.nan(blassoC_aux)]
Results[3,2] <- mean(blassoC_aux[,2,2])-1
Results[3,3] <- mean(abs(blassoC_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoC[,2,2]))
test1 <- sum(abs((blassoC[element_chosen,2,2]-1)/slassoC[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[3,4] <- (500 - test1 - test2)/500

Results[4,1] <- sum(is.nan(blassoCF[,2,2]))
blassoCF_aux <- blassoCF
blassoCF_aux[is.nan(blassoCF_aux)] <- b2sls1[is.nan(blassoCF_aux)]
Results[4,2] <- mean(blassoCF_aux[,2,2])-1
Results[4,3] <- mean(abs(blassoCF_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCF[,2,2]))
test1 <- sum(abs((blassoCF[element_chosen,2,2]-1)/slassoCF[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[4,4] <- (500 - test1 - test2)/500


Results[5,1] <- sum(is.nan(blassoCn[,2,2]))
blassoCn_aux <- blassoCn
blassoCn_aux[is.nan(blassoCn_aux)] <- b2sls1[is.nan(blassoCn_aux)]
Results[5,2] <- mean(blassoCn_aux[,2,2])-1
Results[5,3] <- mean(abs(blassoCn_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCn[,2,2]))
test1 <- sum(abs((blassoCn[element_chosen,2,2]-1)/slassoCn[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[5,4] <- (500 - test1 - test2)/500

Results[6,1] <- sum(is.nan(blassoCFn[,2,2]))
blassoCFn_aux <- blassoCFn
blassoCFn_aux[is.nan(blassoCFn_aux)] <- b2sls1[is.nan(blassoCFn_aux)]
Results[6,2] <- mean(blassoCFn_aux[,2,2])-1
Results[6,3] <- mean(abs(blassoCFn_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCFn[,2,2]))
test1 <- sum(abs((blassoCFn[element_chosen,2,2]-1)/slassoCFn[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[6,4] <- (500 - test1 - test2)/500


Results[7,1] <- sum(is.nan(blassoCX[,2,2]))
blassoCX_aux <- blassoCX
blassoCX_aux[is.nan(blassoCX_aux)] <- b2sls1[is.nan(blassoCX_aux)]
Results[7,2] <- mean(blassoCX_aux[,2,2])-1
Results[7,3] <- mean(abs(blassoCX_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCX[,2,2]))
test1 <- sum(abs((blassoCX[element_chosen,2,2]-1)/slassoCX[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[7,4] <- (500 - test1 - test2)/500

Results[8,1] <- sum(is.nan(blassoCFX[,2,2]))
blassoCFX_aux <- blassoCFX
blassoCFX_aux[is.nan(blassoCFX_aux)] <- b2sls1[is.nan(blassoCFX_aux)]
Results[8,2] <- mean(blassoCFX_aux[,2,2])-1
Results[8,3] <- mean(abs(blassoCFX_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCFX[,2,2]))
test1 <- sum(abs((blassoCFX[element_chosen,2,2]-1)/slassoCFX[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[8,4] <- (500 - test1 - test2)/500


Results[9,1] <- sum(is.nan(RblassoCC[,2,2]))
RblassoCC_aux <- RblassoCC
RblassoCC_aux[is.nan(RblassoCC_aux)] <- b2sls1[is.nan(RblassoCC_aux)]
Results[9,2] <- mean(RblassoCC_aux[,2,2])-1
Results[9,3] <- mean(abs(RblassoCC_aux[,2,2]-1))
element_chosen <- which(!is.nan(RblassoCC[,2,2]))
test1 <- sum(abs((RblassoCC[element_chosen,2,2]-1)/RslassoCC[element_chosen,2,2]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,2,2])+sum(RsupScore05B[-element_chosen,51,2,2]))
Results[9,4] <- (500 - test1 - test2)/500
      
Results[10,1] <- sum(is.nan(RblassoCFC[,2,2]))
RblassoCFC_aux <- RblassoCFC
RblassoCFC_aux[is.nan(RblassoCFC_aux)] <- b2sls1[is.nan(RblassoCFC_aux)]
Results[10,2] <- mean(RblassoCFC_aux[,2,2])-1
Results[10,3] <- mean(abs(RblassoCFC_aux[,2,2]-1))
element_chosen <- which(!is.nan(RblassoCFC[,2,2]))
test1 <- sum(abs((RblassoCFC[element_chosen,2,2]-1)/RslassoCFC[element_chosen,2,2]) < qnorm(0.95))
test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,2,2])+sum(RsupScore05B[-element_chosen,51,2,2]))
Results[10,4] <- (500 - test1 - test2)/500

Results[11,1] <- sum(is.nan(blassoCV[,2,2]))
blassoCV_aux <- blassoCV
blassoCV_aux[is.nan(blassoCV_aux)] <- b2sls1[is.nan(blassoCV_aux)]
Results[11,2] <- mean(blassoCV_aux[,2,2])-1
Results[11,3] <- mean(abs(blassoCV_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCV[,2,2]))
test1 <- sum(abs((blassoCV[element_chosen,2,2]-1)/slassoCV[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[11,4] <- (500 - test1 - test2)/500

Results[12,1] <- sum(is.nan(blassoCFV[,2,2]))
blassoCFV_aux <- blassoCFV
blassoCFV_aux[is.nan(blassoCFV_aux)] <- b2sls1[is.nan(blassoCFV_aux)]
Results[12,2] <- mean(blassoCFV_aux[,2,2])-1
Results[12,3] <- mean(abs(blassoCFV_aux[,2,2]-1))
element_chosen <- which(!is.nan(blassoCFV[,2,2]))
test1 <- sum(abs((blassoCFV[element_chosen,2,2]-1)/slassoCFV[element_chosen,2,2]) < qnorm(0.95))
test2 <- sum(supScore05[-element_chosen,51,2,2])
Results[12,4] <- (500 - test1 - test2)/500

Results_250_180 <- Results

#############################

# Writing to Excel

#write.xlsx(x = Results, file = "Simulation_Results.xlsx", sheetName = "Results")

wb <- createWorkbook()
sheet  <- createSheet(wb, sheetName="Results")
addDataFrame(Results_100_30, sheet, startRow=2, startColumn=1)
addDataFrame(Results_250_30, sheet, startRow=16, startColumn=1)
addDataFrame(Results_100_180, sheet, startRow=30, startColumn=1)
addDataFrame(Results_250_180, sheet, startRow=44, startColumn=1)
saveWorkbook(wb, file="Simulation_Results_Lasso_expa.xlsx") 

