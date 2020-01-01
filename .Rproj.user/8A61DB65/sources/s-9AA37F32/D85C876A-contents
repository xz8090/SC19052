
#' @title Calculation of cumulative loss by parametric approximation.
#' @description The prediction model is described in <Risk model: insurance loss prediction based on R>.
#' @param s the possible value of cumulative loss (numeric)
#' @param mu The mean of cumulative loss (numeric)
#' @param sigma The standard deviation of the cumulative loss (numeric)
#' @param k The skewness coefficient of cumulative loss (numeric)
#' @param method The method of parametric approximation (String)
#' @return The probability that the cumulative loss is less than or equal to s (numeric)
#' @examples
#' \dontrun{
#' s <- 1000
#' mu <- 300
#' sigma <- sqrt(60000)
#' k <- 1.2247
#' Fs1 <- ParaApp(s,mu,sigma,k,method="normal")
#' Fs2 <- ParaApp(s,mu,sigma,k,method="Tgamma")
#' Fs3 <- ParaApp(s,mu,sigma,k,method="normal power")
#' Fs4 <- ParaApp(s,mu,sigma,k,method="Wilson Hilfery")
#' }
#' @export
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
