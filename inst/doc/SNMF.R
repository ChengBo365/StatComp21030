## -----------------------------------------------------------------------------
  cbNMF<-function(V,r,maxiter){
  V<-as.matrix(V)
  h<-nrow(V)
  u<-ncol(V)
  W=matrix(runif(h*r),nrow = h,ncol = r)
  H=matrix(runif(r*u),nrow = r,ncol = u)
  for (i in 1:maxiter) {
    W=W*(V%*%t(H))/(W%*%H%*%t(H))
    W=W/(matrix(rep(1,h),h,1)%*%colSums(W))
    H=H*((t(W)%*%V)/(t(W)%*%W%*%H))
  }
  MV=W%*%H
  list(W=W,H=H,MV=MV)
    }

## -----------------------------------------------------------------------------
library(StatComp21030)
load("~/StatComp21030/data/v.rda")
cbNMF(v,r=5,maxiter=1000)

## -----------------------------------------------------------------------------
cbNMF(v,r=10,maxiter=1000)

## -----------------------------------------------------------------------------
SNMF<-function(V,r,maxiter,kong){
  V<-as.matrix(V)
  h<-nrow(V)
  u<-ncol(V)
  W<-matrix(runif(h*r),nrow = h,ncol = r)
  H<-matrix(runif(r*u),nrow = r,ncol = u)
  for (i in 1:maxiter) {
    W<-W*(V%*%t(H))/(W%*%H%*%t(H))
    W<-W/(matrix(rep(1,h),h,1)%*%colSums(W))
    H<-H*((t(W)%*%V)/(t(W)%*%W%*%H))
    y<-t(W)%*%W%*%H-t(W)%*%V
    for (k in 1:r) {
      y[k,]<-y[k,]/c(t(W[,k])%*%W[,k])
    }
    deta0<-0.5*sqrt((H-y)*(H-y))
    for (k in 1:r) {
      deta0[k,]<-deta0[k,]%*%(t(W[,k])%*%W[,k])
    }
    As<-order(deta0,decreasing = TRUE)
    deta0[As]
    H[deta0<deta0[As[kong]]]=0.00001
  }
  MV<-W%*%H
  list(W=W,H=H,MV=MV)
}

## -----------------------------------------------------------------------------
SNMF(v,r=5,maxiter=1000,kong=50)

## -----------------------------------------------------------------------------
SNMF(v,r=10,maxiter=1000,kong=100)

