#' Logistic model with four parameters
#'
#' Logistic model with four parameters (\eqn{y = f+(a-f)/(1+exp(-(x-m)/s))})
#' @param x numeric vector
#' @param m point of inflection parameter
#' @param s scale parameter
#' @param a upper asymptote parameter
#' @param f lower asymptote parameter
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#' x <- 10:100
#' y <- lgs_msaf(x=x,m=50,s=2,a=0.8,f=0.1)
#' y2 <- lgs_msaf(x=x,m=50,s=2,a=0.5,f=0.2)
#' plot(x,y,type="o",ylim=c(0,1))
#' points(x,y2,type="o",col="blue")
#'

lgs_msaf <- function(m,s,a,f,x) {f+(a-f)/(1+exp(-(x-m)/s))}
