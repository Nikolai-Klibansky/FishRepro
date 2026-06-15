#' Power model with three parameters
#'
#' Power model with three parameters (a+b*x^c)
#' @param x independent variable
#' @param a intercept
#' @param b coefficient
#' @param c exponent
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#'

pow3 <- function(x,a,b,c){a+b*x^c}
