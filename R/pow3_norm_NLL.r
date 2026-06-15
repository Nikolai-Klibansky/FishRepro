#' Objective function for fitting three parameter power function with normal error
#'
#' Objective function for \code{\link[FishRepro]{pow3}} model with three parameters, when residual error of y follows a normal distribution. Computes negative log likelihood value. Parameter values are input in a manner compatible with \code{\link[stats]{optim}}.
#' @param params Named vector of parameters `a`, `b`, and `c` in \code{\link[FishRepro]{pow3}} function. `params` must be in order: `a`, `b`, `c`
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
#' pow3_norm_NLL(params=c("a"=gs.a,"b"=gs.b,"c"=gs.c),x=x_ob,y=y_ob)
#'

pow3_norm_NLL <- function(params,x,y){
  a <- params[1]
  b <- params[2]
  c <- params[3]
  y.pred <- pow3(x=x,a=a,b=b,c=c)
  sd.resid <- sd(y-y.pred)
  -sum(dnorm(y,mean=y.pred,sd <- sd.resid,log=TRUE))
}
