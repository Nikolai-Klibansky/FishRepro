#' Objective function for fitting power function with two parameters and negative binomial error
#'
#' Objective function for \code{\link[FishRepro]{pow}} model with two parameters, when residual error follows a negative binomial distribution. Computes negative log likelihood value. Parameter values are input in a manner compatible with \code{\link[stats]{optim}}.
#' @param params Named vector of parameters `a` and `b` in \code{\link[FishRepro]{pow}} function and dispersion parameter `k` passed to the `size` argument of \code{\link[stats]{dnbinom}}. `params` must be in order: `a`, `b`, `k`
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
#' # Generate guesses for a and b parameters based on linear fit of transformed data
#' gs.a <- exp(coef(fit)[[1]])
#' gs.b <- coef(fit)[[2]]
#'
#' # Compute negative log-likelihood for a particular set of parameter values
#' pow_nbinom_NLL(params=c("a"=gs.a,"b"=gs.b,"k"=2),x=x_ob,y=y_ob)
#'

pow_nbinom_NLL <- function(params,x,y){
  a <- params[1]
  b <- params[2]
  k <- params[3]
  y.pred <- FishRepro::pow(x=x,a=a,b=b)
  -sum(dnbinom(y,mu=y.pred,size=k,log=TRUE))
}
