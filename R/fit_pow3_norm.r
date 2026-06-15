#' Fit power model with with three parameters and normal error to y~x data
#'
#' Fit power model \code{\link[FishRepro]{pow3}} with three parameters, where residual error follows a normal distribution.
#' @param x independent variable
#' @param y dependent variable
#' @param args_optim Additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Fit model to batch fecundity~length data
#' x_ob <- data_fb_L$TL
#' y_ob <- data_fb_L$fb
#' fit <- fit_pow3_norm(x=x_ob,y=y_ob,args_optim = list(control=list(maxit=1000)))
#' x_pr <- seq(min(x_ob),max(x_ob),length=100)
#' y_pr <- do.call(pow3,c(list(x=x_pr),as.list(fit$par)))
#' plot(x_ob,y_ob)
#' points(x_pr,y_pr,type="l",lwd=2)
#'
fit_pow3_norm <- function(x,
                          y,
                          args_optim=list()
                          )
{
  # Fit model to data with transformation and linear regression
  fit <- lm(log(y)~log(x))

  gs.a <- 1 # Intercept
  # Generate guesses for b and c parameters based on linear fit of transformed data
  gs.b <- exp(coef(fit)[[1]])
  gs.c <- coef(fit)[[2]]

  # Setup arguments
  args_optim_user <- args_optim
  args_optim_default <- list(fn=pow3_norm_NLL,
                             par=c(a=gs.a, b=gs.b, c=gs.c),
                             x=x,
                             y=y,
                             control=list(parscale=c(a=gs.a, b=gs.b, c=gs.c))
  )
  args_optim <- modifyList(args_optim_default,args_optim_user)

  # Fit model
  do.call(optim,args_optim)
}

