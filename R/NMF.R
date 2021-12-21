#' @title Standard non-negative matrix factorization algorithm.
#' @description Apply Non-negative matrix factorization to a picture.
#' @param V the data
#' @param r The number of columns in W, or the number of rows in H.
#' @param maxiter Number of cycles
#' @return Generated H matrix
#' @importFrom stats runif
#' @examples
#' \dontrun{
#' data(V)
#' H1<-cbNMF(V,30,500)
#' }
#' @export
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
  list(W=w,H=H,MV=MV)
}



#' @title Sparse non-negative matrix factorization using BeSS.
#' @description Apply Sparse nonnegative matrix factorization to a picture.
#' @param V the data
#' @param r the number of columns in W, or the number of rows in H.
#' @param maxiter the number of cycles
#' @param kong t he number of non-zero elements of the matrix H
#' @return Generated H matrix
#' @importFrom stats runif
#' @examples
#' \dontrun{
#' data(V)
#' H2<-SNMF(V,30,500,800)
#' }
#' @export
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
  list(W=w,H=H,MV=MV)
}
