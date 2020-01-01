
#' @title VaR and TVaR plot.
#' @description The prediction model is described in <Risk model: insurance loss prediction based on R>.
#' @param x The possible value of the random variable of loss (NumericVector)
#' @return  A plot
#' @examples
#' \dontrun{
#' x <- rnorm(100,2,1000)
#' VaRplot(x)
#' }
#' @export
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