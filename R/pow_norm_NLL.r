#' Objective function for fitting power function with two parameters and normal error
#'
#' Objective function for \code{\link[FishRepro]{pow}} model with two parameters, when residual error of y follows a normal distribution. Computes negative log likelihood value. Parameter values are input in a manner compatible with \code{\link[stats]{optim}}.
#' @param params Named vector of parameters `a` and `b` in \code{\link[FishRepro]{pow}} function. `params` must be in order: `a`, `b`
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
#' pow_norm_NLL(params=c("a"=gs.a,"b"=gs.b),x=x_ob,y=y_ob)
#'
#' # For fun, you could try looping through some values of the `a` parameter to do your own simple search for the best fit
#' # and plot a likelihood profile
#' a <- seq(gs.a*0.5, gs.a*1.5,length=50)
#' NLL <- a*NA # Initialize results vector
#' for(i in seq_along(a)){
#'  ai <- a[i]
#'  NLL[i] <- pow_norm_NLL(params=c("a"=ai,"b"=gs.b),x=x_ob,y=y_ob)
#'  }
#'  plot(a,NLL,type="o")

pow_norm_NLL <- function(params,x,y){
  a <- params[1]
  b <- params[2]
  y.pred <- pow(x=x,a=a,b=b)
  sd.resid <- sd(y-y.pred)
  -sum(dnorm(y,mean=y.pred,sd <- sd.resid,log=TRUE))
}
