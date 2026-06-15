#' Logistic model with three parameters
#'
#' Logistic model with three parameters (\eqn{y = a/(1+exp(-(x-m)/s))})
#' @param x numeric vector
#' @param m point of inflection parameter
#' @param s scale parameter
#' @param a upper asymptote parameter
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#' x <- 10:100
#' y <- lgs_msa(x=x,m=50,s=2,a=0.8)
#' y2 <- lgs_msa(x=x,m=50,s=2,a=0.5)
#' plot(x,y,type="o",ylim=c(0,1))
#' points(x,y2,type="o",col="blue")
#'

lgs_msa <- function(x,m,s,a) {a/(1+exp(-(x-m)/s))}
