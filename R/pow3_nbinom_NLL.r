#' Objective function for fitting three parameter power function with negative binomial error
#'
#' Objective function for \code{\link[FishRepro]{pow3}} model with three parameters, when residual error of y follows a negative binomial distribution. Computes negative log likelihood value. Parameter values are input in a manner compatible with \code{\link[stats]{optim}}.
#' @param params Named vector of parameters `a`, `b`, and `c` in \code{\link[FishRepro]{pow3}} function and dispersion parameter `k` passed to the `size` argument of \code{\link[stats]{dnbinom}}. `params` must be in order: `a`, `b`, `c`, `k`
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @author Nikolai Klibansky
#' @export
#' @examples
#' x_ob <- data_fb_L$TL
#' y_ob <- data_fb_L$fb
#'
#' fit <- lm(log(y_ob)~log(x_ob))
#'
#' # Generate guesses for a, b, and c parameters
#' gs.a <- 1
#' gs.b <- exp(coef(fit)[[1]])
#' gs.c <- coef(fit)[[2]]
#'
#' # Compute negative log-likelihood for a particular set of parameter values
#' pow3_nbinom_NLL(params=c("a"=gs.a,"b"=gs.b,"c"=gs.c,"k"=2),x=x_ob,y=y_ob)
#'

pow3_nbinom_NLL <- function(params,x,y){
  a <- params[1]
  b <- params[2]
  c <- params[3]
  k <- params[4]
  y.pred <- FishRepro::pow3(x=x,a=a,b=b,c=c)
  -sum(dnbinom(y,mu=y.pred,size=k,log=TRUE))
}
