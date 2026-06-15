#' Logistic model with two parameters
#'
#' Logistic model with two parameters (\eqn{y = 1/(1+exp(-(x-m)/s))})
#' @param x numeric vector
#' @param m point of inflection parameter
#' @param s scale parameter
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#' x <- 10:100
#' y <- lgs_ms(x=x,m=50,s=2)
#' y2 <- lgs_ms(x=x,m=50,s=5)
#' plot(x,y,type="o")
#' points(x,y2,type="o",col="blue")
#'

lgs_ms <- function(x,m,s) {1/(1+exp(-(x-m)/s))}
