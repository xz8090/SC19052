
#' @title Calculation of cumulative loss by Panjer recursion.
#' @description The prediction model is described in <Modern actuarial risk theory - based on R>.
#' @param p the probability distribution (NumericVector)
#' @param lambda The parameters of the poisson distribution (numeric)
#' @return The maximum possible cumulative loss and the probability of each value (list)
#' @examples
#' \dontrun{
#' p <- c(0.25,0.5,0.25)
#' lambda <- 4
#' fs <- Panjer.Poisson(p,lambda)
#' f <- Panjer.Poisson(p,lambda)$f
#' s <- Panjer.Poisson(p,lambda)$s
#' }
#' @export
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
