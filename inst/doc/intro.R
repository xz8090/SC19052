## -----------------------------------------------------------------------------
ParaApp <- function(s,mu,sigma,k,method){
  if(method=="normal"){
    Fs <- pnorm((s-mu)/sigma)
  }
  if(method=="Tgamma"){
    alpha <- 4/k^2
    beta <- 2/(k*sigma) 
    x0 <- mu-2*sigma/k
    Fs <- pgamma(s-x0,shape=alpha,rate=beta)
  }
  if(method=="normal power"){
    Fs <- pnorm(-3/k+sqrt(9/k^2+1+6/k*(s-mu)/sigma))
  }
  if(method=="Wilson Hilfery"){
    Fs <- pnorm(3*(2/k)^(2/3)*((s-mu)/sigma+2/k)^(1/3)-6/k+k/6)
  }
  return(Fs)
}

## -----------------------------------------------------------------------------
s <- c(0:1000)
mu <- 300
sigma <- sqrt(60000)
k <- 1.2247
Fs1 <- ParaApp(s,mu,sigma,k,method="normal")
Fs2 <- ParaApp(s,mu,sigma,k,method="Tgamma")
Fs3 <- ParaApp(s,mu,sigma,k,method="normal power")
Fs4 <- ParaApp(s,mu,sigma,k,method="Wilson Hilfery")
plot(s,Fs1,pch="",ylab="Fs")
lines(s,Fs1,col='red')
lines(s,Fs2,col='blue')
lines(s,Fs3,col='green')
lines(s,Fs4,col='yellow')

## -----------------------------------------------------------------------------
Panjer.Poisson <- function(p,lambda){
  if(sum(p)>1||any(p<0)) stop("p parameter not a density")
  if(lambda*sum(p)>727) stop("Underflow")
  cumul <- f <- exp(-lambda*sum(p))
  r <- length(p)
  s <- 0
  repeat{
    s <- s+1
    m <- min(s,r)
    last <- lambda/s*sum(1:m*head(p,m)*rev(tail(f,m)))
    f <- c(f,last)
    cumul <- cumul+last
    if(cumul>0.99999999)
    break
  }
  return(list(f=f,s=s))
}

## -----------------------------------------------------------------------------
fs <- Panjer.Poisson(c(0.25,0.5,0.25),4)
s <- c(0:fs$s)
print(cbind(s,fs$f))

## -----------------------------------------------------------------------------
SparseVec <- function(p,lambda){
  freq <- p*lambda
  if (any(freq<0)) stop("negative frequency")
  M <- length(freq)
  mu <- sum((1:M)*freq)
  sigma2 <- sum((1:M)^2*freq)
  MM <- ceiling(mu+10*sqrt(sigma2))+6
  fs <- dpois(0:(MM-1),freq[1])
  for(j in 2:M){
    MMM <- trunc((MM-1)/j)
    fj <- rep(0,MM)
    fj[(0:MMM)*j+1] <- dpois(0:MMM,freq[j])
    fs <- convolve(fs,rev(fj),type="o")
  }
  return(list(f=fs,s=length(fs)-1))
}

## -----------------------------------------------------------------------------
f <- SparseVec(c(0.25,0.5,0.25),4)
print(f)

## -----------------------------------------------------------------------------
VaRplot <- function(x){
p <- seq(0.01,0.99,0.001)
n <- length(p)

VaR=TVaR=NULL
for(j in 1:n){
  VaR[j] <- quantile(x,p[j])
  TVaR[j] <- VaR[j]+mean((x-VaR[j])*(x>VaR[j]))/(1-p[j])
}

plot(p,VaR,type="s",ylab="",lty=1)
lines(p,TVaR,type="s",lty=2)
legend("topleft",c("VaR","TVaR"),lty=c(1,2),bty="n")
}

## -----------------------------------------------------------------------------
x <- rnorm(100,2,1000)
VaRplot(x)

