#' Fit power model with with two parameters and normal error to y~x data
#'
#' Fit power model \code{\link[FishRepro]{pow}} with two parameters, where residual error follows a normal distribution.
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @param args_optim Additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Fit model to batch fecundity~length data
#' x_ob <- data_fb_L$TL
#' y_ob <- data_fb_L$fb
#' fit <- fit_pow_norm(x=x_ob,y=y_ob,args_optim = list(control=list(maxit=1000)))
#' x_pr <- seq(min(x_ob),max(x_ob),length=100)
#' y_pr <- do.call(pow,c(list(x=x_pr),as.list(fit$par)))
#' plot(x_ob,y_ob)
#' points(x_pr,y_pr,type="l",lwd=2)
#'
fit_pow_norm <- function(x,y,args_optim=list()){
  # Fit model to data with transformation and linear regression
  fit <- lm(log(y)~log(x))
  # Generate guesses for a and b parameters based on linear fit of transformed data
  gs.a <- exp(coef(fit)[[1]])
  gs.b <- coef(fit)[[2]]

  # Setup arguments
  args_optim_user <- args_optim
  args_optim_default <- list(fn=pow_norm_NLL,
                             par=c(a=gs.a, b=gs.b),
                             x=x,
                             y=y,
                             control=list(parscale=c(gs.a, gs.b))
  )
  args_optim <- modifyList(args_optim_default,args_optim_user)

  # Fit model
  do.call(optim,args_optim)
}

