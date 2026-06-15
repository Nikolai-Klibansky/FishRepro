#' Von Bertalanffy (VB) length to age function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param L lengths to convert to ages
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @author Nikolai Klibansky
#' @export
#' @examples
#' L <- seq(100,1000, by=100)
#' a <- vb_age(L=L,Linf=1200,K=0.3,t0=-0.5)
#' plot(a,L,type="o")
#'
vb_age <- function(L,Linf,K,t0) {(log(1-L/Linf)/-K)-t0}
