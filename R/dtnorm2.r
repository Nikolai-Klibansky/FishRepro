#' Probability density of a variable following a truncated-normal distribution
#'
#' Probability density of a variable following a truncated-normal distribution. This function calls msm::dtnorm internally but has different behavior.
#' @param n number of quantiles to r
#' @param a coefficient
#' @param b exponent
#' @author Nikolai Klibansky
#' @export
#' @examples
#'
#' plot(dtnorm2(),type="o")
#' plot(dtnorm2(integer=FALSE),type="o")
#'

dtnorm2 <- function(n=101,
                     mean=0,
                     sd=1,
                     alpha=0.05,
                     integers=TRUE){
  CL <- qnorm(p=c(alpha/2,1-alpha/2),mean=mean,sd=sd)
  if(integers){
    x <- round(CL[1]):round(CL[2])
  }else{
    x <- seq(CL[1],CL[2],length=n)}
  y <- pmin(msm::dtnorm(x,mean,sd,lower=CL[1],upper=CL[2]),1) # Constrain to probability of <=1
  data <- data.frame("x"=x,"y"=y)
  return(data)
}
