#' Fit power model with two parameters and negative binomial error to y~x data.
#'
#' Fit power model \code{\link[FishRepro]{pow}} with two parameters, where residual error follows a negative binomial distribution.
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @param klim Interval to search identify starting value of dispersion parameter k. Passed to the interval argument in \code{\link[stats]{optimize}}.
#' @param args_optim Additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Fit model to batch fecundity~length data
#' x_ob <- data_fb_L$TL
#' y_ob <- data_fb_L$fb
#' fit <- fit_pow_nbinom(x=x_ob,y=y_ob)
#' x_pr <- seq(min(x_ob),max(x_ob),length=100)
#' y_pr <- do.call(pow,c(list(x=x_pr),as.list(fit$par[c("a","b")])))
#' plot(x_ob,y_ob)
#' points(x_pr,y_pr,type="l",lwd=2)
#'

fit_pow_nbinom <- function(x,y,klim=c(0,10),args_optim=list()){
  # Define objective function to use with optimize() to generate starting value of k
  pow_nbinom_NLL_optimize <- function(k,x, y, a, b){
    y.pred <- pow(x=x, a=a, b=b)
    out <- c()
    for(k.i in k){
      out <- c(out,-sum(dnbinom(x=y,mu=y.pred,size=k.i,log=TRUE)))
    }
    return(out)
  }

  # Fit model to data with transformation and linear regression
  fit <- lm(log(y)~log(x))

  # Generate guesses for a and b parameters based on linear fit of transformed data
  gs.a <- exp(coef(fit)[[1]])
  gs.b <- coef(fit)[[2]]
  # Generate guess for k parameter fixing a and b parameters with guesses
  gs.k <- optimize(f = pow_nbinom_NLL_optimize, interval=klim,
                   x=x, y=y, a=gs.a, b=gs.b)$minimum

  # Setup arguments
  args_optim_user <- args_optim
  args_optim_default <- list(fn=pow_nbinom_NLL,
                             par=c(a=gs.a, b=gs.b, k=gs.k),
                             x=x,
                             y=y,
                             control=list(parscale=c(gs.a, gs.b, gs.k))
  )

  args_optim <- modifyList(args_optim_default,args_optim_user)

  # Fit model
  do.call(optim,args_optim)

}
