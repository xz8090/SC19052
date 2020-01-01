## -----------------------------------------------------------------------------
x<- seq(-10,10,length=20)
y<- x
f<- function(x,y) x^2+y^2
z<- outer(x,y,f)
persp(x,y,z,theta=30,phi=20,expand=0.8,col="lightblue")

## -----------------------------------------------------------------------------
n <- 100
x <- rnorm(n)
y <- 3*x + rnorm(n)
plot(x, y)
out <- lm(y ~ x)
library(pander)
panderOptions("digits", 4)
pander(out)

## -----------------------------------------------------------------------------
rrayl <- function(n,sigma) {
u<- runif(n)
v <- sqrt(-2*(sigma)^2*log(1-u))
return(v)
}
x<-rrayl(10,1)
print(x)

## -----------------------------------------------------------------------------
x <- rrayl(1000,1)
sigma <-1
hist(x, prob = TRUE, main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
lines(density(x),col="red")
y <- seq(0, 4, .04)
z <- y/((sigma)^2)*exp(-y^2/(2*(sigma)^2))
lines(y,z)

## -----------------------------------------------------------------------------
x <- rrayl(1000,4)
sigma <-4
hist(x, prob = TRUE, main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
lines(density(x),col="red")
y <- seq(0, 15, .1)
z <- y/((sigma)^2)*exp(-y^2/(2*(sigma)^2))
lines(y,z)

## -----------------------------------------------------------------------------
x <- rrayl(1000,10)
sigma <-10
hist(x, prob = TRUE, main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
lines(density(x),col="red")
y <- seq(0, 40, .4)
z <- y/((sigma)^2)*exp(-y^2/(2*(sigma)^2))
lines(y,z)

## -----------------------------------------------------------------------------
x <- rrayl(1000,100)
sigma <-100
hist(x, prob = TRUE, main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
lines(density(x),col="red")
y <- seq(0, 400, 4)
z <- y/((sigma)^2)*exp(-y^2/(2*(sigma)^2))
lines(y,z)

## -----------------------------------------------------------------------------
n <- 1000
X <- rnorm(n)
y <- rnorm(n,3,1)
Z1<- 0.75*x +0.25*y
p1 <- sample(c(0,1),n,replace=TRUE)
Z2 <- p1*x+(1-p1)*y
hist(Z1)
hist(Z2)

## -----------------------------------------------------------------------------
n <- 1000
x <- rnorm(n)
y <- rnorm(n,3,1)
W1<- 0.2*x +0.8*y
W2<- 0.4*x +0.6*y
W3<- 0.5*x +0.5*y
W4<- 0.6*x +0.4*y
W5<- 0.8*x +0.2*y
W6<- 0.7*x +0.3*y
hist(W1)
hist(W2)
hist(W3)
hist(W4)
hist(W5)
hist(W6)

## -----------------------------------------------------------------------------
rWishart <- function(n,p,Sigma) {
L <- chol(Sigma) #Choleski factorization of Sigma
A <- matrix(rnorm(p*p), nrow=p, ncol=p) 
B <- A
B[lower.tri(A)]<-0 # B is the mataix of the upper trig elements in A.
C <- A-B           # C is a lower triangular matrix with 0 diagonal elements.
E <- array(1:p)    
for(i in 1:p) {
E[i] <- sqrt(rchisq(1,n-i+1))} # E is a vector with E[i]~chisq(1,n-i+1)
G <-  C+diag(E,nrow=p) 
X <- L%*%G%*%t(G)%*%t(L) # The Bartlett decomposition of a matrix X
return(X)
}

## -----------------------------------------------------------------------------
Y <-rWishart(4,2,matrix(c(1,0.9,0.9,1),nrow=2,ncol=2))
print(Y)
print(cov(Y))

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 1e4 
t <- runif(n, min=0, max=(pi)/3)
estimator <- mean(sin(t)) * (pi)/3
exact.value <- cos(0)-cos((pi)/3)
print(estimator)
print(exact.value)

## -----------------------------------------------------------------------------
MC.Anti <- function(x,R = 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic) v <- runif(R/2) else
v <- 1 - u
u <- c(u, v)
g <- exp(-u)/(1+u^2)
cdf<- mean(g)*x
cdf
}
x <- 1
set.seed(12345)
MC1<-MC.Anti(x,antithetic = FALSE) # generate a simple estimator 
set.seed(12345)
MC2<-MC.Anti(x) # generate an antithetic variable estimator
print(c(MC1,MC2))
set.seed(12345)
m <- 1000
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
MC1[i] <- MC.Anti(1,R=1000,antithetic= FALSE)
MC2[i] <- MC.Anti(1,R=1000)
}
print((var(MC1) - var(MC2))/var(MC1)) #compute the approximate reduction in variance                                         as a percentage of the variance

## -----------------------------------------------------------------------------
set.seed(12345)
M <- 10000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0, N, 2)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- runif(M) #f3, inverse transform method
x <- - log ( (1 - exp(-1))*u)
h <- function(x){g(x) / (exp(-x) / (1 - exp(-1)))}
for (i in 1:N) {
estimates[i, 1] <- mean(g(runif(M)))
for (j in 1:k)
T2[j] <- mean(h(runif(M/k, (j-1)/k, j/k)))
estimates[i, 2] <- mean(T2)
}
apply(estimates,2,mean)
apply(estimates,2,var)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 20 # simple size
alpha <- .05
m <- 10000 
sum <- 0 # count as required
for(i in 1 : m){
x <- rchisq(n,2)
LCL <- mean(x)-sd(x)*qt(1-alpha/2,df=n-1)/sqrt(n) 
# generate the lower bound of the confidence interval
UCL <- mean(x)+sd(x)*qt(1-alpha/2,df=n-1)/sqrt(n) 
# generate the upper bound of the confidence interval
if(LCL<2&UCL>2)
  sum <- sum+1
else
  sum <- sum+0
} 
print(c(mean(LCL),mean(UCL))) # the confidence interval
print(sum) #count the number of intervals that contain mu=2
print(sum/m) # compute the mean to get the confidence level

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 1000 # generate the number of cycles
b <- numeric(m)
for(i in 1:m){
x <- rnorm(1000)
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
b[i] <-sk(x)
}
quantile(b,c(.025,.05,.95,.975)) # compute the quantiles of the skewness
q <- c(.025,.05,.95,.975)
Q <- unname(quantile(b,c(.025,.05,.95,.975))) # just take the result below the percent
f <- function(x){1/sqrt(2*pi*6/m)*exp(-x^2/2)} 
# the density function of normal distribution
Var <- q*(1-q)/m/(f(Q))^2 # according to the formula(2.14)
se  <- sqrt(Var/m)
print(cbind(q,se))

## -----------------------------------------------------------------------------
n <- 1000
cv <- qnorm(c(.025,.05,.95,.975),0,sqrt(6/n))
print(cv)

## -----------------------------------------------------------------------------
set.seed(123)
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
n <- c(10,20,30,50,100,500)
cv <- qnorm(.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))
power <- numeric(length(n))
m<-1000
for (i in 1 : length(n)){
sktests <- numeric(m)
for (j in 1:m){
  alpha <- 100
 x <- rbeta(n[i],alpha,alpha) #generate random numbers from Beta distribution
 sktests[j] <- as.integer(abs(sk(x))>=cv[i]) #test decision is 1 (reject) or 0
}
power[i]<- mean(sktests) #proportion rejicted
}
print(rbind(n,power))

## -----------------------------------------------------------------------------
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
set.seed(123)
alpha <- .05
n <- 30
m <- 2000
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
power <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
a<- sample(c(10, 500), replace = TRUE,
size = n, prob = c(1-e, e))
x <- rbeta(n,a,a)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
}
power[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, power, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="Power of Beta distribution")
abline(h = .1, lty = 3)

## -----------------------------------------------------------------------------
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
set.seed(123)
alpha <- .05
n <- 30
m <- 2000
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
power <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
a<- sample(c(1, 10), replace = TRUE,
size = n, prob = c(1-e, e))
x <- rt(n,a)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
}
power[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, power, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="Power of t distribution")
abline(h = .1, lty = 3)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 100
alpha <- .05
mu0 <- 1
m <- 10000
p <- numeric(m)
for(j in 1:m){
  x <- rchisq(n,mu0)
  ttest <- t.test(x,alternative="greater",mu=mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p<=alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
set.seed(123)
n <- 100
alpha <- .05
mu0 <- 1
m <- 10000
p <- numeric(m)
for(j in 1:m){
  x <- runif(n,0,2*mu0)
  ttest <- t.test(x,alternative="greater",mu=mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p<=alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
set.seed(123)
n <- 100
alpha <- .05
mu0 <- 1
m <- 10000
p <- numeric(m)
for(j in 1:m){
  x <- rexp(n,mu0)
  ttest <- t.test(x,alternative="greater",mu=mu0)
  p[j] <- ttest$p.value
}
p.hat <- mean(p<=alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
library(boot)
data(scor,package="bootstrap")
pairs(~.,data=scor,main="Scatterplot Matrix",col="blue")
cor(scor)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(123)
B <- 200 #number of replicates
n <- nrow(scor) #sample size
R12 <- R34 <- R35 <- R45 <- numeric(B) #storage for replicates
for (b in 1:B) {
#randomly select the indices
i <- sample(1:n, size = n, replace = TRUE)
mec <- scor$mec[i] #i is a vector of indices
vec <- scor$vec[i]
alg <- scor$alg[i]
ana <- scor$ana[i]
sta <- scor$sta[i]
R12[b] <- cor(mec,vec)
R34[b] <- cor(alg,ana)
R35[b] <- cor(alg,sta)
R45[b] <- cor(ana,sta)
}
#output
cat('SE.R12 =',sd(R12),
    'SE.R34 =',sd(R34),
    'SE.R35 =',sd(R35),
    'SE.R45 =',sd(R45))

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
sk <- function(data,indices){
x <- data[indices]
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
data <- rnorm(100)
boot.obj <- boot(data,statistic=sk,R=2000)
print(boot.ci(boot.obj,type=c("norm","basic","perc")))

## -----------------------------------------------------------------------------
mu <- 0 # skewness of normal distribution
m <- 1000
n <- 100
library(boot)
set.seed(123)
boot.sk <- function(x,i){ #function to compute the skewness statistic
  x <- x[i]
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  U <- rnorm(n)
  boot.obj <- boot(data=U,statistic=boot.sk,R=999)
  ci <- boot.ci(boot.obj,type=c("norm","basic","perc"))
ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
ci.perc[i,]<-ci$percent[4:5]
}
Pr.nrom <- c(mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),mean(ci.norm[,1]>=mu),
             mean(ci.norm[,2]<=mu))
Pr.basic <- c(mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),mean(ci.basic[,1]>=mu),
             mean(ci.basic[,2]<=mu))
Pr.perc <- c(mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),mean(ci.perc[,1]>=mu),
             mean(ci.perc[,2]<=mu))
data.frame(Pr.nrom,Pr.basic,Pr.perc,row.names=c("coverage","miss.left","miss.right"))

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
sk <- function(data,indices){
x <- data[indices]
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
data <- rchisq(100,5)
boot.obj <- boot(data,statistic=sk,R=2000)
print(boot.ci(boot.obj,type=c("norm","basic","perc")))

## -----------------------------------------------------------------------------
library(moments)
set.seed(123)
x <- rchisq(100,5)
mu <- skewness(x) # skewness of ??2(5) distribution
library(boot)
m <- 1000
n <- 100
boot.sk <- function(x,i){
  x <- x[i]
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  U <- rchisq(n,5)
  boot.obj <- boot(data=U,statistic=boot.sk,R=999)
  ci <- boot.ci(boot.obj,type=c("norm","basic","perc"))
ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
ci.perc[i,]<-ci$percent[4:5]
}
Pr.nrom <- c(mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),mean(ci.norm[,1]>=mu),
             mean(ci.norm[,2]<=mu))
Pr.basic <- c(mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),mean(ci.basic[,1]>=mu),
             mean(ci.basic[,2]<=mu))
Pr.perc <- c(mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),mean(ci.perc[,1]>=mu),
             mean(ci.perc[,2]<=mu))
data.frame(Pr.nrom,Pr.basic,Pr.perc,row.names=c("coverage","miss.left","miss.right"))

## -----------------------------------------------------------------------------
First.pc <- function(data){
  m <- nrow(data)
  cov.MLE <- cov(data)*(m-1)/m
  evals <- eigen(cov.MLE)$values
  First.pc <- evals[1]/sum(evals[1:5])
  return(First.pc)
}
library(boot)
data(scor,package="bootstrap")
theta.hat <- First.pc(scor)
n <- nrow(scor)
theta.jack <- numeric(n)
for(i in 1:n){
  theta.jack[i] <- First.pc(scor[-i,])
}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- sqrt((n-1) *mean((theta.jack - mean(theta.jack))^2))
print(c(theta.hat,bias,se))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)
L4 <- lm(magnetic ~ chemical + I(chemical^2)+I(chemical^3))
plot(chemical,magnetic, main="Cubic", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] *a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)
cat('A.Rs-L1 =',summary(L1)$adj.r.squared,
    'A.Rs-L2 =',summary(L2)$adj.r.squared,
    'A.Rs-L3 =',summary(L3)$adj.r.squared,
    'A.Rs-L4 =',summary(L4)$adj.r.squared)


## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
J4 <- lm(y ~ x + I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
}
print(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)))
cat('E-J1 =',mean(e1^2),
    'E-J2 =',mean(e2^2),
    'E-J3 =',mean(e3^2),
    'E-J4 =',mean(e4^2))

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000
alphahat <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
count5test(x, y)
}))
print(alphahat)

## -----------------------------------------------------------------------------
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(max(c(outx, outy)))
} #counts the maximum number of extreme points
set.seed(12345)
n1 <- 20
n2 <- 30
x <- rnorm(n1,0,1)
y <- rnorm(n2,0,1)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:50
C <- numeric(R) #storage for replicates
options(warn = -1)
C0 <- maxout(x, y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = 20, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
C[i] <- maxout(x1, y1)
}
p <- mean(c(C0, C) >= C0)
options(warn = 0)
print(p)

## -----------------------------------------------------------------------------
dCov <- function(x, y) {
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m)
return(b + M)
}
A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
dCov
}
ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q1 <- dims[2] + 1
d <- p + dims[2]
x <- z[ , 1:p] #leave x as is
y <- z[ix, q1:d] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
library(boot)
m <- 10
p.value <- numeric(m)
n <- c(seq(10,150,10))
N <- length(n)
power1 <- numeric(N)
for(i in 1:N){
  for(j in 1:m){
X <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
e <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
Y <- (1/4)*X+e
z <- cbind(X,Y)
boot.obj <- boot(data = z, statistic = ndCov2, R = 99,
sim = "permutation", dims = c(2, 2))
tb <- c(boot.obj$t0, boot.obj$t)
p.value[j] <- mean(tb >= boot.obj$t0)
  }
power1[i] <- mean(p.value<.05)  
}
plot(n, power1, type = "b",pch=15,lty=1,col="DarkTurquoise",xlab = bquote(n), ylim = c(0,1),main="Power Comparsion-Model 1")

m <- 10
p.ball <- numeric(m)
n <- c(seq(10,150,10))
N <- length(n)
power2 <- numeric(N)
for(i in 1:N){
for(j in 1:m){
X <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
e <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
Y <- (1/4)*X+e
library(Ball)
p.ball[j] <- bcov.test(X,Y,R=99,seed=i*123)$p.value
}
power2[i] <- mean(p.ball<.05)
}
print(cbind(n,power1,power2))
par(new=TRUE)
plot(n, power2, type = "o",pch=16,lty=2,col="DeepPink",xlab = bquote(n), ylim = c(0,1))
legend(120,0.3,c("distance","ball"),col=c("DarkTurquoise","DeepPink"),text.col=c("DarkTurquoise","DeepPink"),pch=c(15,16),lty=c(1,2))

## -----------------------------------------------------------------------------
dCov <- function(x, y) {
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m)
return(b + M)
}
A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
dCov
}
ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q1 <- dims[2] + 1
d <- p + dims[2]
x <- z[ , 1:p] #leave x as is
y <- z[ix, q1:d] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
library(boot)
m <- 10
p.value <- numeric(m)
n <- c(seq(10,150,10))
N <- length(n)
power1 <- numeric(N)
for(i in 1:N){
  for(j in 1:m){
X <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
e <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
Y <- (1/4)*X*e
z <- cbind(X,Y)
boot.obj <- boot(data = z, statistic = ndCov2, R = 99,
sim = "permutation", dims = c(2, 2))
tb <- c(boot.obj$t0, boot.obj$t)
p.value[j] <- mean(tb >= boot.obj$t0)
  }
power1[i] <- mean(p.value<.05)  
}
plot(n, power1, type = "b",pch=15,lty=1,col="DarkTurquoise",xlab = bquote(n), ylim = c(0,1),main="Power Comparsion-Model 2")

m <- 10
p.ball <- numeric(m)
n <- c(seq(10,150,10))
N <- length(n)
power2 <- numeric(N)
for(i in 1:N){
for(j in 1:m){
X <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
e <- matrix(rnorm(n[i]*2),nrow=n[i],ncol=2)
Y <- (1/4)*X*e
library(Ball)
p.ball[j] <- bcov.test(X,Y,R=99,seed=i*123)$p.value
}
power2[i] <- mean(p.ball<.05)
}
print(cbind(n,power1,power2))
par(new=TRUE)
plot(n, power2, type = "o",pch=16,lty=2,col="DeepPink",xlab = bquote(n), ylim = c(0,1))
legend(120,0.3,c("distance","ball"),col=c("DarkTurquoise","DeepPink"),text.col=c("DarkTurquoise","DeepPink"),pch=c(15,16),lty=c(1,2))

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic)
dslap <- function(x){(1/2)*exp(-abs(x))}
x <- dskewlap(25)
y <- dslap(25)
print(c(x,y))

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic) 
rw.Metropolis <- function( sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (dskewlap(y) / dskewlap(x[i-1]))) #dskewlap() can compute the density of                                                   Lapace distribution 
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
} #generate the sample from standard Laplace distribution with rejected number  
set.seed(123)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
rej.num <- c(rw1$k, rw2$k, rw3$k, rw4$k)
acce.rate <- c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N)
print(cbind(sigma,rej.num,acce.rate))

refline <- qskewlap(c(.025, .975)) #qskewlap() can compute the quantile of                                                   Lapace distribution
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
 plot((rw)[,j], type="l",
xlab=bquote(sigma == .(round(sigma[j],3))),
ylab="X", ylim=range(rw[,j]))
abline(h=refline)
}

## -----------------------------------------------------------------------------
x <- c(seq(10,100,10))
y <- log(exp(x))
z <- exp(log(x))
y==z
isTRUE(all.equal(y,z))

## -----------------------------------------------------------------------------
# creat a function
c.k<-function(k,a){return(sqrt((k*(a^2))/(k+1-a^2)))}
f1<-function(u){(1+(u^2)/(k-1))^(-k/2)}
f2<-function(u){(1+(u^2)/k)^(-(k+1)/2)}
f<-function(a){
(2*gamma(k/2))/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(f1,0,c.k(k-1,a))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(f2,0,c.k(k,a))$value
}  
# compute roots
library(rootSolve)
t<-c(4:25,100)
n<-length(t)
root<-root2<-numeric(n)
for (i in 1:n) {
  k<-t[i]
  root[i]=uniroot(f,c(0.05,sqrt(k)/2+1))$root
}
f2<-function(a){
  pt(sqrt((k-1)*a^2/(k-a^2)),k-1)-pt(sqrt((k*a^2)/(k+1-a^2)),k)
}
f.4<-function(a){
  pt(sqrt((k-1)*a^2/(k-a^2)),k-1)-pt(sqrt((k*a^2)/(k+1-a^2)),k)
}
K1<-c(4:25,100,500,1000)
n<-length(K1)
root.4<-numeric(n)
for (i in 1:n) {
  k<-K1[i]
  root.4[i]<-uniroot(f.4,c(0.5,sqrt(k)/2+1))$root
}
print(cbind(t,root,root.4),) ##the roots of 11.5 and 11.4


## -----------------------------------------------------------------------------
N<-1e3
# max. number of the iteration
n1<-28;n2<-24;n3<-41;n4<-70
L<-c(.5,.4)
# initial estimates 
tol<-.Machine$double.eps^0.5
L.old<-L+1
E<-numeric(N)
for(j in 1:N){
  E[j]<-2*L[1]*n1*log(L[1])/(2-L[1]-2*L[2])+2*L[2]*n2*log(L[2])/(2-L[2]-2*L[1])+2*n3*log(1-L[1]-L[2])+n1*(2-2*L[1]-2*L[2])*log(2*L[1]*(1-L[1]-L[2]))/(2-L[1]-2*L[2])+n2*(2-2*L[1]-2*L[2])*log(2*L[2]*(1-L[1]-L[2]))/(2-L[2]-2*L[1])+n4*log(2*L[1]*L[2])
  model<-function(x){
    F1<-2*L[1]*n1/((2-L[1]-2*L[2])*x[1])-2*n3/(1-x[1]-x[2])+n1*(2-2*L[1]-2*L[2])*(1-2*x[1]-x[2])/((2-L[1]-2*L[2])*x[1]*(1-x[1]-x[2]))-n2*(2-2*L[1]-2*L[2])/((2-L[2]-2*L[1])*(1-x[1]-x[2]))+n4/x[1]
    F2<-2*L[2]*n2/((2-L[2]-2*L[1])*x[2])-2*n3/(1-x[1]-x[2])-n1*(2-2*L[1]-2*L[2])/((2-L[1]-2*L[2])*(1-x[1]-x[2]))+n2*(2-2*L[1]-2*L[2])*(1-2*x[2]-x[1])/((2-L[2]-2*L[1])*x[2]*(1-x[1]-x[2]))+n4/x[2]
    c(F1=F1,F2=F2)
  }
  ss<-multiroot(f=model,star=c(.1,.1))
  L<-ss$root
  # update p and q
  if (sum(abs(L-L.old)/L.old)<tol) break
  L.old<-L
}
print(L.old) #the estimator of p and q
plot(E,type = "l",main="Log-maximum likelihood values in M-steps")

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
out <- vector("list", length(formulas))
for (i in seq_along(formulas)) {
out[[i]] <- lm(formulas[[i]])
}
print(out)

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
lapply(formulas,lm)


## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
Mtc <-mtcars[rows, ]
Mtc$mpg~Mtc$disp
})
out <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)) {
out[[i]] <- lm(bootstraps[[i]])
}
print(out)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
Mtc <-mtcars[rows, ]
lm (Mtc$mpg~Mtc$disp)
})
print(bootstraps )

## -----------------------------------------------------------------------------
bootstraps <- out <-  vector("list",10)
for (i in seq_along(out)) {
 rows <- sample(1:nrow(mtcars), rep = TRUE)
 Mtc <-mtcars[rows, ]
 bootstraps[[i]] <- Mtc$mpg~Mtc$disp
 out[[i]] <- lm(bootstraps[[i]])
}
print(out)

## -----------------------------------------------------------------------------
mod <- list(
lm(mpg ~ disp),
lm(mpg ~ I(1 / disp)),
lm(mpg ~ disp + wt),
lm(mpg ~ I(1 / disp) + wt)
)
rsq <- function(mod) summary(mod)$r.squared
lapply(mod,rsq) # extract R2 for model in Exercise 3

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
Mtc <-mtcars[rows, ]
L <- lm (Mtc$mpg~Mtc$disp)
rsq <- function(mod) summary(mod)$r.squared
bootstraps[i] <- rsq(L)
})
print(bootstraps) # extract R2 for model in Exercise 4

## -----------------------------------------------------------------------------

p.value <- sapply(1:100,function(i){
  p.value <- t.test(rpois(10, 10), rpois(7, 10))$p.value
 })
print(p.value) # extract the p-value from every trial

## -----------------------------------------------------------------------------
sapply(1:10, sqrt)
library(parallel)
cl <- makeCluster(getOption("cl.cores",4))
parSapply(cl,1:10, sqrt)

## -----------------------------------------------------------------------------
fun <- function(x){
  return(x+1)
}
system.time(sapply(1:10000, fun))

library(parallel)
cl <- makeCluster(getOption("cl.cores",4))
system.time(parSapply(cl,1:10000, fun))
stopCluster(cl)

## -----------------------------------------------------------------------------
mcsapply <- function(x, f, ...) {
out <- vector("list", length(x))
for (i in sample(seq_along(x))) {
out[[i]] <- f(x[[i]], ...)
}
out
}

boot_df <- function(x) x[sample(nrow(x), rep = T), ]
rsquared <- function(mod) summary(mod)$r.square
boot_lm <- function(i) {
rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
}
system.time(sapply(1:100, boot_lm))
system.time(mcsapply(1:100, boot_lm))

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- 'C:/Users/Administrator/Desktop/Rcpp/'
# Can create source file in Rstudio
sourceCpp(paste0(dir_cpp,"rwC.cpp"))
set.seed(123)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rwC(sigma[1], x0, N)
rw2 <- rwC(sigma[2], x0, N)
rw3 <- rwC(sigma[3], x0, N)
rw4 <- rwC(sigma[4], x0, N)
rej.num <- c(rw1$k, rw2$k, rw3$k, rw4$k)
acce.rate <- c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N)
print(cbind(sigma,rej.num,acce.rate))

library(GeneralizedHyperbolic) 
refline <- qskewlap(c(.025, .975)) 
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
 plot((rw)[,j], type="l",
xlab=bquote(sigma == .(round(sigma[j],3))),
ylab="X", ylim=range(rw[,j]))
abline(h=refline)
}

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic) 
rw.Metropolis <- function( sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (dskewlap(y) / dskewlap(x[i-1])))
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}
set.seed(123)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rwR1 <- rw.Metropolis(sigma[1], x0, N)
rwR2 <- rw.Metropolis(sigma[2], x0, N)
rwR3 <- rw.Metropolis(sigma[3], x0, N)
rwR4 <- rw.Metropolis(sigma[4], x0, N)
rwR <- cbind(rwR1$x, rwR2$x, rwR3$x, rwR4$x)

library(Rcpp)
dir_cpp <- 'C:/Users/Administrator/Desktop/Rcpp/'
sourceCpp(paste0(dir_cpp,"rwC.cpp"))
rwC1 <- rwC(sigma[1], x0, N)
rwC2 <- rwC(sigma[2], x0, N)
rwC3 <- rwC(sigma[3], x0, N)
rwC4 <- rwC(sigma[4], x0, N)
rwC <- cbind(rwC1$x, rwC2$x, rwC3$x, rwC4$x)

for (j in 1:4){
  qqplot(rwR[,j],rwC[,j],xlab='rwR',ylab='rwC')
  abline(0,b=1)
  }

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- 'C:/Users/Administrator/Desktop/Rcpp/'
sourceCpp(paste0(dir_cpp,"rwC.cpp"))
library(microbenchmark)
N <- 1000
sigma <- .05
x0 <- 25
ts <- microbenchmark(R=rw.Metropolis(sigma,x0,N),cpp=rwC(sigma,x0,N))
summary(ts)[,c(1,3,5,6)]

