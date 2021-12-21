## -----------------------------------------------------------------------------
#Some simple codes are written as follows:
x1<-c(1.2,2.6, 3.3,4.5,6.1,6.9,7.5,8.4,9.8,10.1,11.2,12.8,13.5,14.4,15.4,16.7)
x2<-c(0.98,0.87,0.85,0.80,0.74,0.69,0.65,0.60,0.52,0.49,0.44,0.39,0.28,0.19,0.15,0.08)
x3<-c(12.6,13.1,14,15.1,16.8,17.2,18.3,19.1,20.1,21.2,22.3,23.5,24.6,25.4,26.7,27.8)
y<-c(53,56,57,58.9,59.9,60.4,61.9,62.8,63.9, 65, 66.4,68.2,69.9,70.8,72.2,74)
fun1<-lm(y~I(x1^2)+x1+x2+x3)
fun1
summary(fun1)$coef 

## ----islands------------------------------------------------------------------
#Let's take a look at the data about the island by using a table.
knitr::kable(head(islands))

## ----rock---------------------------------------------------------------------
knitr::kable(head(rock))

## -----------------------------------------------------------------------------
#the graph of fun1 is shown as follows:
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
plot(fun1)

## -----------------------------------------------------------------------------
 x4=matrix(rnorm(100,mean=1.5,sd=2),ncol=2)
 plot(x4)

## -----------------------------------------------------------------------------
#Develop an algorithm to generate random samples from a Rayleigh(σ) distribution
n <- 1000
a<-runif(10,0,3) #取随机数
for (i in seq_along(a)) {
  sigma <- a[[i]]
  u <- runif(n)
  x <- (-2*sigma^2*log(u))^(1/2) # F(x) = 1-exp(-x^2/2/sigma^2)  x>=0
  #作图
  hist(x, prob = TRUE, main = expression(f(x)==x/sigma^2*exp(-x^2/2/sigma^2)))
  y <- seq(0, 8, .01)
  lines(y, y/sigma^2*exp(-y^2/2/sigma^2))
  print(sigma)  #输出结果
}

## -----------------------------------------------------------------------------
# Generate a random sample of size 1000 from a normal location mixture
n <- 1000
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,3,1)
Z1 <- 0.75*X1+0.25*X2
hist(Z1, prob = T)
p<-sample(c(0,1),n,replace = TRUE)
Z2 <- p*X1+(1-p)*X2
hist(Z2, prob = T)
sum(p)/n

## -----------------------------------------------------------------------------
#Write a program to simulate a compound Poisson(λ)–Gamma process (Y has a Gamma distribution)
r<-2;beta<-4
X1<-c()
y1<-c()
t<-10;lambda<-5   #lambda为5;t为10时
for (j in 1:100) {
  N<-rpois(1,t*lambda)  
  y<-rgamma(N,r,beta);
  X1[j]<-(sum(y))
  y1[j]<-(y[1])
}
X1   #模拟产生X1的100个值
mean(X1)
var(X1)
EX<-t*lambda*mean(y1);EX  #理论期望
VarX<-t*lambda*(var(y1)+(mean(y1))^2);VarX  #理论方差

## -----------------------------------------------------------------------------
r<-2;beta<-4
X2<-c()
y1<-c()
t<-10;lambda<-10   #lambda为10;t为10时
for (j in 1:100) {
  N<-rpois(1,t*lambda) 
  y<-rgamma(N,r,beta)
  X2[j]<-(sum(y))
  y1[j]<-(y[1])
}
X2   #模拟产生X2的100个值
mean(X2)
var(X2)
EX<-t*lambda*mean(y1);EX  #理论期望
VarX<-t*lambda*(var(y1)+(mean(y1))^2);VarX  #理论方差

## -----------------------------------------------------------------------------
r<-2;beta<-2
X3<-c()
y1<-c()
t=10;lambda=10  #lambda为10;t为10时
for (j in 1:100) {
  N<-rpois(1,t*lambda) 
  y<-rgamma(N,r,beta)
  X3[j]<-(sum(y))
  y1[j]<-(y[1])
}
X3   #模拟产生X3的100个值
mean(X3)
var(X3)
EX<-t*lambda*mean(y1);EX  #理论期望
VarX<-t*lambda*(var(y1)+(mean(y1))^2);VarX  #理论方差

## -----------------------------------------------------------------------------
#Write a function to compute a Monte Carlo estimate of the Beta(3, 3) cdf
MC.Phi<-function(x){
  m<-1e4
  u<-runif(m,0,x)
  cdf<-mean(x*u^2*(1-u)^2);cdf 
}
#Compare the estimates with the values returned by the pbeta function in R.
a<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for(i in 1:9){
  x<-a[i]
  theta.hat<-MC.Phi(x)/MC.Phi(1)
  theta<-pbeta(x,3,3)
  print(c(theta.hat,theta)) #输出结果
}

## -----------------------------------------------------------------------------
#Implement a function to generate samples from a Rayleigh(σ) distribution,using antithetic variables
MC.P1 <- function(x, R = 10000) {
  a <- runif(R/2)
  v <- 1- a
  u <- c(a, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g<- x[i] *x[i]*u *exp(-(x[i]*u)^2 / 2)
    cdf[i]<-mean(g)
  }
  cdf
}

MC.P2 <- function(x, R = 10000) {
  a <- runif(R/2)
  v <- runif(R/2)
  u<-c(a,v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g<- x[i] *x[i]*u *exp(-(x[i]*u)^2 / 2)
    cdf[i]<-mean(g)
  }
  cdf
}
m <- 1000
MC1<-MC2 <- numeric(m)
x<-3
for (i in 1:m) {
  MC1[i] <- MC.P1(x, R = 1000)
  MC2[i] <- MC.P2(x, R = 1000)
}
V1<-var(MC1);V1
V2<-var(MC2)/2;V2
V2/V1

## -----------------------------------------------------------------------------
#Find two importance functions
x <- seq(0, 1, .01)
w <- 2
g <- x^2*exp(-x^2/2)/sqrt(2*pi)
f1 <- x*exp(x^2/2)/(1-exp(-1/2))
f2 <- x^(1/2)*exp(-x)
gs <- c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),
        expression(f[1](x)==x*e^{x^2/2}/(1-e^{-1/2})),
        expression(f[2](x)==x^{1/2}*e^{-x}))
par(mfrow=c(1,2))
#figure (a)
plot(x, g, type = "l", ylab = "",
     ylim = c(0,1), lwd = w,col=1,main='(A)')
lines(x, f1, lty = 2, lwd = w,col=2)
lines(x, f2, lty = 3, lwd = w,col=3)
legend("topright", legend = gs,
       lty = 1:3, lwd = w, inset = 0.02,col=1:3)
#figure (b)
plot(x, g, type = "l", ylab = "",
     ylim = c(0,2.2), lwd = w, lty = 2,col=2,main='(B)')
lines(x, g/f1, lty = 3, lwd = w,col=3)
lines(x, g/f2, lty = 4, lwd = w,col=4)
legend("topright", legend = gs[-1],
       lty = 2:4, lwd = w, inset = 0.02,col=2:4)

## ----echo=FALSE---------------------------------------------------------------
#Which of your two importance functions should produce the smaller variance in estimating
m <- 10000
me <- sd <- numeric(3)
g <- function(x) {
  x^2*exp(-x^2/2) / sqrt(2*pi) * (x > 1) 
}
x<-rexp(m,1)+1
fg<-g(x)/1
me[1] <- mean(fg)
sd[1] <- sd(fg)

x <- rexp(m,1)+1 #using f1
fg <- g(x) / x/exp(-x^2/2)*(1-exp(-1/2))
me[2] <- mean(fg)
sd[2] <- sd(fg)

x <- rexp(m,1) #using f2
fg <- g(x) / sqrt(x)/exp(-x)
me[3] <- mean(fg)
sd[3] <- sd(fg)

res <- rbind(me=round(me,3), sd=round(sd,3))
colnames(res) <- paste0('f',0:2)
knitr::kable(res, format = "latex",align='c')
print(res)


## -----------------------------------------------------------------------------
#Obtain a Monte Carlo estimate
g1 <- function(x) {
  x^2*exp(-x^2/2) / sqrt(2*pi) * (x > 1) 
}
g2 <- function(x) {
  exp(-1/2/x^2) /x^4/ sqrt(2*pi) * (x>0) *(x<1)
}

m<-1e5
u<-runif(m)
fg <- g2(u)/exp(-u)*(1-exp(-1))  #取f(x)=e^(-x)/(1-e^(-1))
cdf<-mean(fg);cdf


a1<-integrate(g1,1,Inf)
a2<-integrate(g2,0,1)
print(c(a1,a2))

## -----------------------------------------------------------------------------
# a function is written to get confidence interval
wf1<-function(n,alpha){
  x <- rchisq(n,2)
  UCL11 <- mean(x)+sd(x) *qt(alpha/2,n-1)/ sqrt(n)
  UCL12<-mean(x)+ sd(x) *qt(1-alpha/2,n-1)/ sqrt(n)
  c(UCL11,UCL12)
}

## -----------------------------------------------------------------------------
#calculate the probablity of UCL11<2 and UCL12>2
n=10000
data1=matrix(0,nrow = 2,ncol = n)
for(i in 1:n){
  data1[,i]=wf1(20,0.05)
}
length(which(data1[1,]<2&data1[2,]>2))/n
new.t1<-mean(data1[1,])
new.t2<-mean(data1[2,])
new.t1
new.t2

## -----------------------------------------------------------------------------
n=10000
wf1<-function(n,alpha){
  x <- rchisq(n,1)
  UCL11 <- mean(x)+sd(x) *qt(alpha/2,n-1)/ sqrt(n)
  UCL12<-mean(x)+ sd(x) *qt(1-alpha/2,n-1)/ sqrt(n)
  u<-mean(x)
  c(UCL11,UCL12,u)
}
data1=matrix(0,nrow = 3,ncol = n)
for(i in 1:n){
  data1[,i]=wf1(20,0.05)
}
VCL1<-mean(data1[1,])
UCL1<-mean(data1[2,])
U1<-mean(data1[3,])
print(c(VCL1,UCL1,U1))

## -----------------------------------------------------------------------------
n=10000
wf2<-function(n,alpha){
  x <- runif(n,0,2)
  UCL11 <- mean(x)+sd(x) *qt(alpha/2,n-1)/ sqrt(n)
  UCL12<-mean(x)+ sd(x) *qt(1-alpha/2,n-1)/ sqrt(n)
  u<-mean(x)
  c(UCL11,UCL12,u)
}
data2=matrix(0,nrow = 3,ncol = n)
for(i in 1:n){
  data2[,i]=wf2(20,0.05)
}
VCL2<-mean(data2[1,])
UCL2<-mean(data2[2,])
U2<-mean(data2[3,])
print(c(VCL2,UCL2,U2))

## -----------------------------------------------------------------------------
n=10000
wf3<-function(n,alpha){
  x <- rexp(n,1)
  UCL11 <- mean(x)+sd(x) *qt(alpha/2,n-1)/ sqrt(n)
  UCL12<-mean(x)+ sd(x) *qt(1-alpha/2,n-1)/ sqrt(n)
  u<-mean(x)
  c(UCL11,UCL12,u)
}
data3=matrix(0,nrow = 3,ncol = n)
for(i in 1:n){
  data3[,i]=wf3(20,0.05)
}
VCL3<-mean(data3[1,])
UCL3<-mean(data3[2,])
U3<-mean(data3[3,])
print(c(VCL3,UCL3,U3))

## -----------------------------------------------------------------------------
#Skewness chisq.test
n <- c(10, 20, 30, 50, 80, 100) #sample sizes
d<-2    #for ease of calculation,we could choose d=2.
cv <- qchisq(.95,d*(d+1)*(d+2)/6) #crit. values for each n

## -----------------------------------------------------------------------------
#Repeat Examples 6.8 and 6.10 for Mardia’s multivariate skewness test
library(MASS)
#计算b_{1, d}
sk <- function(x,y,n) {    
  S<-matrix(0,n,n)
  xbar <- colMeans(x)
  ybar <- colMeans(y)
  sigma<-cov(x,y)    #计算协方差
  sigma.hat<-ginv(sigma) 
  for(i in 1:n){
    for(j in 1:n){
      S[i,j]<-((x[i,]-xbar)%*%sigma.hat%*%t(t(y[j,]-ybar)))^3
    }
  }
  return(mean(S))
}

## -----------------------------------------------------------------------------
p.reject <- numeric(length(n))  
m <- 1000                   
for (i in 1:length(n)) {
  sktests <- numeric(m)      
  a<-n[i]
  for (j in 1:m) {
    x <- matrix(rnorm(a*2),nrow=a,ncol=2)   #生成X，Y随机阵
    y <- matrix(rnorm(a*2),nrow=a,ncol=2)
    sktests[j] <- as.integer(sk(x,y,a) >= cv ) 
  }
  p.reject[i] <- mean(sktests)         #计算拒绝P值
}
p.reject


## -----------------------------------------------------------------------------
#Power of the skewness test of chi-squared
alpha <- 0.1
n <- 25
m <- 2500
epsilon <- c(seq(0, 0.15, 0.01), seq(0.15, 1, 0.05))
N <- length(epsilon)
pw <- numeric(N)
cv <- qchisq(.95, d*(d+1)*(d+2)/6)
for (j in 1:N) { 
  e <- epsilon[j]
  skt<- numeric(m)
  for (i in 1:m) { 
    sigma <- sample(c(1, 10), replace = TRUE,size = n, prob = c(1-e, e))
    a <- rnorm(n, 0, sigma);b <- rnorm(n, 0, sigma)
    x <- matrix(a,nrow=5,ncol=5)
    y <- matrix(b,nrow=5,ncol=5)
    skt[i] <- as.integer(sk(x,y,5) >= cv)
  }
  pw[j] <- mean(skt)
}
se <- sqrt(pw * (1-pw) / m) #add standard errors
pw;se
#plot power vs epsilon
plot(epsilon, pw, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = 0.1, lty = 3)

## -----------------------------------------------------------------------------
#Use bootstrap to estimate the bias and standard error of θ
set.seed(0)   #生成种子
library(boot)  #加载程辑包
library(bootstrap)
lambda_hat <- eigen(cov(scor))$values #初始lambda值
the.hat <- lambda_hat[1] / sum(lambda_hat)
B <- 2000 # number of bootstrap samples
n <- nrow(scor) # number of rows (data size)
# 进行Bootstrap抽样
func <- function(dat, index){
  x <- dat[index,]
  lambda <- eigen(cov(x))$values
  theta <- lambda[1] / sum(lambda)
  return(theta)
}
bootstrap_result <- boot(
  data = cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),
  statistic = func, R = B)
theta_b <- bootstrap_result$t
bias_boot <- mean(theta_b) - the.hat
# the estimated bias of theta_hat, using bootstrap
se_boot <- sqrt(var(theta_b))
# the estimated standard error (se) of theta_hat, using bootstrap
round(c(original=the.hat,estimate=mean(theta_b),bias_boot=bias_boot,se_boot=se_boot),5)

## -----------------------------------------------------------------------------
#calculate the data covariance matrix
library(bootstrap) #for the score data
cov(scor)
#show a part of data
scor[1:5,]

## -----------------------------------------------------------------------------
#A function that defines a covariance matrix
mcov <- function(x,j) cov(x[j,])

## -----------------------------------------------------------------------------
#Computing data dimension
n<-lengths(scor)[1]
#Calculate the real covariance matrix
real.cov <- mcov(scor,1:n)
#Calculate the real theta value
the.hat<-eigen(real.cov)$values[1]/sum(eigen(real.cov)$values)
#Define variables to store values
the.jack <- numeric(n)
#Loop n times
for(i in 1:n){
  the.cov <- mcov(scor,(1:n)[-i])
  #Get matrix eigenvalues
  eig1<-eigen(the.cov)$values
  #Find the estimate of lambda
  the.jack[i]<-eig1[1]/sum(eig1)
}
#Calculation bias
bias1 <- (n-1)*(mean(the.jack)-the.hat)
#Calculating  standard deviation
se1 <- sqrt((n-1)*mean((the.jack-the.hat)^2))
round(c(original=the.hat,estimate=mean(the.jack),bias1=bias1,
        se1=se1),4)

## ----eval=FALSE---------------------------------------------------------------
#  #Compute 95% percentile and BCa confidence intervals for θ
#  set.seed(0)
#  library(boot)
#  library(bootstrap)
#  lambda_hat <- eigen(cov(scor))$values
#  theta_hat <- lambda_hat[1] / sum(lambda_hat)
#  B <- 2000 # number of bootstrap samples
#  n <- nrow(scor) # number of rows (data size)
#  # Bootstrap
#  func <- function(dat, index){
#    x <- dat[index,]
#    lambda <- eigen(cov(x))$values
#    theta <- lambda[1] / sum(lambda)
#    return(theta)
#  }
#  #对数据进行B=2000次Bootstrap抽样
#  bootstrap_result <- boot(
#    data = cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),
#    statistic = func, R = B)
#  ci<-boot.ci(bootstrap_result,type=c("perc","bca"))
#  ci.perc<-ci$percent[4:5]   #计算得到95%的percentile和BCa的置信区间
#  ci.bca<-ci$bca[4:5]
#  round(c(ci.perc,ci.bca),4)

## -----------------------------------------------------------------------------
#Calculate the skewness with the following function
sk <- function(x,j) {
  #computes the sample skewness coeff.
  xbar <- mean(x[j])
  m3 <- mean((x[j] - xbar)^3)
  m2 <- mean((x[j] - xbar)^2)
  return( m3 / m2^1.5 )
}

## ----eval=FALSE---------------------------------------------------------------
#  #Compare the coverage rates
#  me<-0
#  n<-20
#  m<-1000
#  set.seed(12345)
#  library(boot)
#  nor.norm<-nor.basic<-nor.perc<-matrix(NA,m,2)
#  for (i in 1:m) {
#    data.nor<-rnorm(n,0,2)
#    nor.ske<-boot(data.nor,statistic=sk,R=1000)
#    nor<- boot.ci(nor.ske,type=c("norm","basic","perc"))
#    nor.norm[i,]<-nor$norm[2:3];
#    nor.basic[i,]<-nor$basic[4:5];
#    nor.perc[i,]<-nor$percent[4:5];
#  }
#  #Calculate the coverage probability of a normal distribution
#  cat('norm =',mean(nor.norm[,1]<=me & nor.norm[,2]>=me),
#      'basic =',mean(nor.basic[,1]<=me & nor.basic[,2]>=me),
#      'perc =',mean(nor.perc[,1]<=me & nor.perc[,2]>=me))
#  #Calculate the probability of the left side of the normal distribution
#  cat('norm.left=',mean(nor.norm[,1]>=me ),
#      'basic.left =',mean(nor.basic[,1]>=me ),
#      'perc.left =',mean(nor.perc[,1]>=me))
#  #Calculate the right side probability of a normal distribution
#  cat('norm.right=',mean(nor.norm[,2]<=me ),
#      'basic.right =',mean(nor.basic[,2]<=me ),
#      'perc.right =',mean(nor.perc[,2]<=me))

## -----------------------------------------------------------------------------
pian.real<-mean(replicate(1000,expr = {
  x0<-rchisq(1000,5)
  mean((x0-5)^3)/10^1.5
}))
print(pian.real)

## ----eval=FALSE---------------------------------------------------------------
#  me<-pian.real
#  n<-20
#  m<-1000
#  set.seed(12345)
#  library(boot)
#  chi.norm<-chi.basic<-chi.perc<-matrix(NA,m,2)
#  for (i in 1:m) {
#    data.chisq<-rchisq(n,5)
#    chisq.ske<-boot(data.chisq,statistic=sk,R=1000)
#    chi<- boot.ci(chisq.ske,type=c("norm","basic","perc"))
#    chi.norm[i,]<-chi$norm[2:3];
#    chi.basic[i,]<-chi$basic[4:5];
#    chi.perc[i,]<-chi$percent[4:5];
#  }
#  #Calculate the coverage probability of the chi-square distribution
#  cat('norm =',mean(chi.norm[,1]<=me & chi.norm[,2]>=me),
#      'basic =',mean(chi.basic[,1]<=me & chi.basic[,2]>=me),
#      'perc =',mean(chi.perc[,1]<=me & chi.perc[,2]>=me))
#  #Calculate the probability of the left side of the chi-square distribution
#  cat('norm.left=',mean(chi.norm[,1]>=me ),
#      'basic.left =',mean(chi.basic[,1]>=me ),
#      'perc.left =',mean(chi.perc[,1]>=me))
#  #Calculate the right side probability of the chi-square distribution
#  cat('norm.right=',mean(chi.norm[,2]<=me ),
#      'basic.right =',mean(chi.basic[,2]<=me ),
#      'perc.right =',mean(chi.perc[,2]<=me))

## ----eval=FALSE---------------------------------------------------------------
#  # Implement the bivariate Spearman rank correlation test for independence as a permutation test
#  #自行生成0到10的十个随机数
#  x<-runif(10,0,10)
#  y<-runif(10,0,10)
#  R <- 999 #重复次数
#  z <- c(x, y) #集中样本
#  N <- 1:20
#  reps <- numeric(R)
#  t0 <- cor.test(x, y)$statistic
#  for (i in 1:R) {
#    n <- sample(N, size = 10, replace = FALSE)  #产生样本指标
#    x1 <- z[n]
#    y1 <- z[-n] #z中拿走x1后剩下的给y1
#    reps[i] <- cor.test(x1, y1)$statistic #相关系数检验
#  }
#  p <- mean(abs(c(t0, reps)) >= abs(t0))
#  p_value<-cor.test(x,y)$p.value
#  round(c(p=p,p_value=p_value),3)

## -----------------------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
#NN检验函数
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z)
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 <= n1); i2 <- sum(block2 > n1)
  (i1 + i2) / (k * n)
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e3
#  k<-3
#  set.seed(12345)
#  n1 <- n2 <- 25
#  R<-999
#  n <- n1+n2
#  N = c(n1,n2)
#  p_values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*2,0,1),ncol=2)  #生成均值为0，方差为1的随机阵
#    y <- matrix(rnorm(n1*2,0,2),ncol=2)  #生成均值为0，方差为2的随机阵
#    z <- rbind(x,y)
#    p_values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p_values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p_values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.1
#  pow <- colMeans(p_values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e3
#  k<-3
#  set.seed(12345)
#  n1 <- n2 <- 25
#  R<-999
#  n <- n1+n2
#  N = c(n1,n2)
#  p_values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*2,0,1),ncol=2)  #生成均值为0，方差为1的随机阵
#    y <- matrix(rnorm(n1*2,1,2),ncol=2)  #生成均值为1，方差为2的随机阵
#    z <- rbind(x,y)
#    p_values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p_values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p_values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.1
#  pow <- colMeans(p_values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e3
#  k<-3
#  set.seed(12345)
#  n1 <- n2 <- 25
#  R<-999
#  n <- n1+n2
#  N = c(n1,n2)
#  p_values <- matrix(NA,m,3)
#  for(i in 1:m){
#    #t distribution with 1 df
#    x <- matrix(rt(n1*2,1),ncol=2)  #生成自由度为1的t分布的随机阵
#    y <- matrix(rt(n1*2,1),ncol=2)  #生成自由度为1的t分布的随机阵
#    z <- rbind(x,y)
#    p_values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p_values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p_values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.05
#  pow <- colMeans(p_values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e3
#  k<-3
#  set.seed(12345)
#  n1 <- n2 <- 25
#  R<-999
#  n <- n1+n2
#  N = c(n1,n2)
#  p_values <- matrix(NA,m,3)
#  for(i in 1:m){
#    #bimodel distribution
#    a<-rnorm(n1*2+n2*2);b<-rnorm(n1*2+n2*2,1,2);
#    M1<-0.3*a+0.7*b
#    x <- matrix(M1,nrow=n1,ncol=2)  #mixture of two normal distributions
#    y <- matrix(M1,nrow=n2,ncol=2)  #mixture of two normal distributions
#    z <- rbind(x,y)
#    p_values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p_values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p_values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.05
#  pow <- colMeans(p_values<alpha)
#  pow

## -----------------------------------------------------------------------------
m <- 1e3
k<-3
set.seed(12345)
n1 <- 5;n2 <- 50
R<-999
n <- n1+n2
N = c(n1,n2)
p_values <- matrix(NA,m,3)
for(i in 1:m){
  a<-rnorm(20);b<-rnorm(200,1,2);
  #分别抽取10个和100个样本
  m1 <- sample(a, size = 10, replace = FALSE)   
  m2 <- sample(b, size = 100, replace = FALSE)
  #生成随机阵
  x <- matrix(m1,ncol=2)  
  y <- matrix(m2,ncol=2)  
  z <- rbind(x,y)
  p_values[i,1] <- eqdist.nn(z,N,k)$p.value
  p_values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p_values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
}
alpha <- 0.1
pow <- colMeans(p_values<alpha)
pow

## -----------------------------------------------------------------------------
f <- function(x) {
  return(1 / (pi*(1+x^2)))
}

## -----------------------------------------------------------------------------
#Use the Metropolis-Hastings sampler to generate random variables from a standard Cauchy distribution
set.seed(123)
m <- 10000
sigma <- 4
x <- numeric(m)
x[1] <- rnorm(1)
k <- 0
u <- runif(m)
for (i in 2:m) {
  xt <- x[i-1]
  y <- rnorm(1, 0, abs(xt))
  num <- f(y) * dnorm(xt, 0, abs(y))
  den <- f(xt) * dnorm(y,0, abs(xt))
  r<-num/den
  if(r>1) r<-1
  if (u[i] <= num/den) x[i] <- y else {
    x[i] <- xt
    k <- k+1 #y is rejected
  }
}
print(k)

## -----------------------------------------------------------------------------
#compare the deciles of the generated observations with the deciles of the standard Cauchy distribution
index <- 1001:10000
y0 <- x[index]
y1<-quantile(y0, c(1:9)/10)
y2<-qcauchy(c(1:9)/10,0,1) 
rbind(y1,y2)

## -----------------------------------------------------------------------------
#Use the Gibbs sampler to generate a chain with target joint density f(x, y)
set.seed(123)
m <- 1000
n <- 20   #n取20,a,b取2
a <- 2
b <- 2
x <- matrix(0,m, 2)
x[1,] <- c(0, 0.5)
for (i in 2:m) {
  x2 <- x[i-1, 2]
  x[i,1] <- rbinom(1, n, x2)
  x1 <- x[i,1]
  x[i,2] <- rbeta(1, x1+a, n-x1+b)
}
x <- x[101:m, ]
colMeans(x)
plot(x[,1],x[,2],xlab = "x",ylab = "y")

## -----------------------------------------------------------------------------
set.seed(123)
#生成卡方马氏链函数
f.chain <- function(m,x1){
  a <- 2;b <- 2;n <- 20
  x <- matrix(0, m, 2)
  x[1,] <- x1 
  for (i in 2:m) {
    x2 <- x[i-1, 2]
    #生成二项分布的边际发布
    x[i,1] <- rbinom(1, n, x2)
    x1 <- x[i,1]
    #生成贝塔分布的边际分布
    x[i,2] <- rbeta(1, x1+a, n-x1+b)
  }
  return(x)
}

## -----------------------------------------------------------------------------
#a Write a function to compute
set.seed(1234)
f<-function(a,d,k){
  b<-(-1)^k/2^k/gamma(k+1)
  c<-(sum(a*a))^(2*k+2)/(2*k+1)/(2*k+2)
  e<-gamma((d+1)/2)*gamma(k+3/2)/gamma(k+d/2+1)
  s<-b*c*e
  return(s)
}

## -----------------------------------------------------------------------------
#b Modify the function so that it computes and returns the sum
g<-function(n){
  s=0
  for (k in 0:n){
    s = s + f(a,d,k)
    return(s)
  }
}

## -----------------------------------------------------------------------------
#c Evaluate the sum
a<-c(1,2)
d<-2
n<-1000
s<-g(n)
print(s)

## -----------------------------------------------------------------------------
#First define the function to be called
rm(list=ls())
k=4
ck<-function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
}
wf1 <- function(u) {
  (1+u*u/(k-1))^(-k/2)
}
wf2<- function(u) {
  (1+u*u/k)^(-(k+1)/2)
}
wf3<-function(a){
  2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(wf1,0,ck(k-1,a))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(wf2,0,ck(k,a))$value
}

## -----------------------------------------------------------------------------
#the root value of a under different values of k
kt<-c(seq(4,25,1),100,200)
solution11.5<-solution11.4<-numeric(length(kt))
js=1
for (i in kt) {
  k<-kt[js]
  solution11.5[js]<-uniroot(wf3,c(0.001,sqrt(k)/2+1))$root
  js=js+1
}

## -----------------------------------------------------------------------------
#the root value of a with different values of k
js=1
wf4<-function(a){pt(sqrt(a^2*(k-1)/(k-a^2)),k-1)-pt(sqrt(a^2*k/(k+1-a^2)),k)}
for (i in kt) {
  k<-kt[js]
  solution11.4[js]<-uniroot(wf4,c(0.00001,sqrt(k)-0.0001))$root
  js=js+1
}
data.frame(kt,solution11.4,solution11.5)

## -----------------------------------------------------------------------------
#EM algorithe
set.seed(1234)
N<-100
lambda<-numeric(N)
n<-10;n0<-7
y<-c(0.54,0.48,0.33,0.43,1.0,1.0,0.91,1.0,0.21,0.85)
for(i in 1:N){
  lambda[1]<-0
  lambda[i+1]<-mean(y)+(n-n0)/n*lambda[i]   #迭代过程
  if(abs(lambda[i+1]-lambda[i])<1e-6) break  #设定误差为小于1e-6
  lambda_bar<-lambda[i+1]
}

## -----------------------------------------------------------------------------
#MLE
set.seed(1234)
n<-10;n0<-7;tau<-1
y0<-c(0.54,0.48,0.33,0.43,0.91,0.21,0.85)
lambda_hat<-(tau*(n-n0)+sum(y0))/n0
round(c(lambda_bar,lambda_hat),7)

## -----------------------------------------------------------------------------
# For the first model
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
wf <- numeric(4)
rsq <- function(mod) summary(mod)$r.squared
j = 1:4
wf <- lapply(j,function(j) {rsq(lm(formulas[[j]], data = mtcars))})
wf

## -----------------------------------------------------------------------------
# For the second model
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
wf <- numeric(10)
rsq <- function(mod) summary(mod)$r.squared
j = 1:10
wf <- lapply(j,function(j) {rsq(lm(formulas[[1]], data = bootstraps[[j]]))})
wf


## -----------------------------------------------------------------------------
#a
su<-function(x){
  funs<-c(sd)
  sapply(funs,function(f) f(x,na.rm=TRUE))
}
df<-data.frame(replicate(3,runif(10,1,10)))
round(vapply(df,su,FUN.VALUE = c(sd=0)),3)

## -----------------------------------------------------------------------------
#b
su<-function(x){
  funs<-c(sd)
  sapply(funs,function(f) f(x,na.rm=TRUE))
}
df<-data.frame(x<-replicate(3,runif(10,1,10)),c(list("A","B","C")))
round(vapply(df,su,FUN.VALUE = c(sd=0)),3)

## -----------------------------------------------------------------------------
library(parallel)
formulas <- list(
  l1<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp),l2<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp)),l3<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp + mtcars$wt),l4<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp) +mtcars$wt)
)
cl<-makeCluster(4)
bootstraps2 <- lapply(1:3000, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]})
system.time(sapply(bootstraps2,formulas[[1]]))
system.time(parSapply(cl,bootstraps2,formulas[[1]]))

## -----------------------------------------------------------------------------
#生成R函数
R_Gibbs<-function(n,a,b){
  N <- 10000 #连的长度
  X <- matrix(0, N, 2) #定义一个二元变量
  X[1, ] <- c(1, 0.5) #初始值
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    #生成二项分布的随机数
    X[i, 1] <- rbinom(1,n,x2)
    x1 <- X[i, 1]
    #生成贝塔分布的随机数
    X[i, 2] <- rbeta(1, x1+a, n-x1+b)
  }

  return(X)
}

## -----------------------------------------------------------------------------
#生成Rcpp函数
library(Rcpp)
cppFunction('NumericMatrix C_Gibbs(int n,double a, double b){
            int N=10000;
            int x;
            double y;
            NumericMatrix M(N,2);
            M(0,0)=1;
            M(0,1)=0.5;
            for(int i=1;i<N;i++){
            y=M(i-1,1);
            M(i,0)=rbinom(1,n,y)[0];
            x=M(i,0);
            M(i,1)=rbeta(1,x+a,n-x+b)[0];
            }
            return(M) ;
            }')


## -----------------------------------------------------------------------------
set.seed(0)#设定随机种子
burn<-1000
N<-10000 #循环次数
GibbsR=R_Gibbs(15,3,3) #调用函数
GibbsR<-GibbsR[(burn+1):N,]
GibbsC=C_Gibbs(15,3,3)
GibbsC<-GibbsC[(burn+1):N,]
#分别抽取边际变量值
GibbsR_x<-GibbsR[,1]
GibbsR_y<-GibbsR[,2]

GibbsC_x<-GibbsC[,1]
GibbsC_y<-GibbsC[,2]
#画QQ图
qqplot(GibbsR_x,GibbsC_x)
qqplot(GibbsR_y,GibbsC_y)


## -----------------------------------------------------------------------------
library(microbenchmark)
#比较计算速度
time<-microbenchmark(GibbsR=R_Gibbs(15,3,3),GibbsC=C_Gibbs(15,3,3))
summary(time)[, c(1,3,5,6)]

