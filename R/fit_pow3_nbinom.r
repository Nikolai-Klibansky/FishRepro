#' Fit power model with three parameters and negative binomial error to y~x data.
#'
#' Fit power model \code{\link[FishRepro]{pow3}} with three parameters, where residual error follows a negative binomial distribution.
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @param klim Interval to search identify starting value of overdispersion parameter k. Passed to the interval argument in \code{\link[stats]{optimize}}.
#' @param args_optim Additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Fit model to batch fecundity~length data
#' x_ob <- data_fb_L$TL
#' y_ob <- data_fb_L$fb
#' fit <- fit_pow3_nbinom(x=x_ob,y=y_ob)
#' x_pr <- seq(min(x_ob),max(x_ob),length=100)
#' y_pr <- do.call(pow3,c(list(x=x_pr),as.list(fit$par[c("a","b","c")])))
#' plot(x_ob,y_ob)
#' points(x_pr,y_pr,type="l",lwd=2)
#'

fit_pow3_nbinom <- function(x,y,klim=c(0,10),args_optim=list()){
  # Define objective function to use with optimize() to generate starting value of k
  pow3_nbinom_NLL_optimize <- function(k, x, y, a, b, c){
    y.pred <- pow3(x=x, a=a, b=b, c=c)
    out <- c()
    for(k.i in k){
      out <- c(out,-sum(dnbinom(x=y,mu=y.pred,size=k.i,log=TRUE)))
    }
    return(out)
  }

  # Fit model to data with transformation and linear regression
  fit <- lm(log(y)~log(x))

  # Generate guesses for a and b parameters based on linear fit of transformed data
  gs.a <- 1 # Intercept
  gs.b <- exp(coef(fit)[[1]])
  gs.c <- coef(fit)[[2]]
  # Generate guess for k parameter fixing a and b parameters with guesses
  gs.k <- optimize(f = pow3_nbinom_NLL_optimize, interval=klim,
                   x=x, y=y, a=gs.a, b=gs.b, c=gs.c)$minimum

  # Setup arguments
  args_optim_user <- args_optim
  args_optim_default <- list(fn=pow3_nbinom_NLL,
                             par=c(a=gs.a, b=gs.b, c=gs.c, k=gs.k),
                             x=x,
                             y=y,
                             control=list(parscale=c(a=gs.a, b=gs.b, c=gs.c, k=gs.k))
  )

  args_optim <- modifyList(args_optim_default,args_optim_user)

  # Fit model
  do.call(optim,args_optim)
}
