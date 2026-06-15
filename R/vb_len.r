#' Von Bertalanffy (VB) age to length function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param a ages to convert to lengths
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param aP proportion of age (value between 0 and 1) at which to compute length (e.g. aP=0.5 to compute length at midyear)
#' @author Nikolai Klibansky
#' @export
#' @examples
#' a <- 1:20
#' L <- vb_len(a=a,Linf=1200,K=0.3,t0=-0.5)
#' plot(a,L,type="o")
#'
vb_len <- function(a,Linf,K,t0,aP=0)  {Linf*(1-exp(-K*(a+aP-t0)))}
