#' Fit plateau model with four parameters
#'
#' Fit plateau model with four parameters \code{\link[FishRepro]{pla_msda}}, where residual error follows a negative binomial distribution.
#' @param x numeric vector. A continuous variable
#' @param y numeric vector. A binary variable (0 or 1)
#' @param args_optim list of additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @note A plateau model is type of double logistic model where the slope of the first (left) logistic function is positive and the slope of the second (right) logistic function is negative, resulting in a domed shape. If the slopes are steep enough and the points of inflection are far enough apart, then the model will have a flat-topped section between the points of inflection.
#' @export
#' @examples
#' data <- FishRepro::data_SI_d_bin
#' data$Start.date <- as.Date(data$Start.date)
#' data$doy <- as.numeric(format(data$Start.date,"%j"))
#' x <- data$doy
#' y <- data$SI
#'
#' # Compute some reasonable initial parameter values
#' x2 <- x[as.logical(y)]
#' tab <- table(x,y)
#' ptab <- prop.table(tab,1)
#' Q <- quantile(x2,probs=c(0.25,0.75))
#'
#' mn_x <- mean(x2)
#' sd_x <- sd(x2)
#' m.init <- Q[[1]]
#' s.init <- sqrt((3*(sd_x^2))/(pi^2)) # use relationship between scale parameter and sd of a logistic distribution
#' d.init <- as.numeric(diff(Q))
#' a.init <- quantile(ptab[,"1"],probs=0.95)[[1]]
#'
#' plot(rownames(ptab),ptab[,"1"],xlab="doy",ylab="proportion with SI")
#' abline(v=m.init+c(0,d.init),lty=2)
#' abline(h=a.init,lty=2)
#'
#' x_pr <- seq(min(x),max(x),by=1)
#' y_pr <- pla_msda(x=x_pr,m=m.init,s=s.init,d=d.init,a=a.init)
#'
#' # plot predicition based on guesses
#' points(x_pr,y_pr,type="l",col="green")
#' # fit model
#' fit <- fit_pla_msda(x=x,y=y,args_optim=list(par=c("m"=m.init, "s"=s.init, "d"=d.init, "a"=a.init)))
#' # plot predictions
#' curve(expr=do.call(pla_msda,c(list(x=x),as.list(fit$par))),
#'       add=TRUE,lwd=2,lty=1,col="blue")
#' legend("topright",legend=c("observed","predicted"),
#'        pch=c(16,NA),lwd=2,col="blue",bty="n")
#' #'

fit_pla_msda <- function(x,
                         y,
                        args_optim=list()
){
  # Compute some reasonable initial parameter values
  x2 <- x[as.logical(y)]
  tab <- table(x,y)
  ptab <- prop.table(tab,1)
  Q <- quantile(x2,probs=c(0.25,0.75))

  mn_x <- mean(x2)
  sd_x <- sd(x2)
  m.init <- Q[[1]]
  s.init <- sqrt((3*(sd_x^2))/(pi^2)) # use relationship between scale parameter and sd of a logistic distribution
  d.init <- as.numeric(diff(Q))
  a.init <- quantile(ptab[,"1"],probs=0.95)[[1]]

  par.init <- c(m=m.init, s=s.init, d=d.init, a=a.init)

#  Setup arguments
args_optim_user <- args_optim

# If initial par values were given by the user but they didn't give new values
# for control$parscale, use the par values
if("par"%in%names(args_optim_user)&is.null(args_optim$control$parscale)){
  args_optim_user$control$parscale <- args_optim_user$par
}
args_optim_default <- list(fn=pla_msda_NLL,
                           par=par.init,
                           x=x,
                           y=y,
                           control=list(parscale=par.init),
                           method="BFGS",
                           hessian=TRUE
)

args_optim <- modifyList(args_optim_default,args_optim_user)

# Fit model
do.call(optim,args_optim)
}
