#' Exponentially transformed linear model
#'
#' Exponentially transformed linear model exp(a+b*x)
#' @param x numeric vector
#' @param a intercept
#' @param b coefficient
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#'

exp_lin <- function(x,a,b){exp(a+b*x)}
