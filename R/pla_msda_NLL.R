#' Objective function for fitting plateau model with four parameters
#'
#' Objective function for fitting plateau model with four parameters \code{\link[FishRepro]{pla_msda}}.
#' @param params A vector of initial values of the parameters of \code{\link[FishRepro]{pla_msda}} to be optimized over: `m`, `s`, `d`, and `a`. Passed to the \code{par} argument in \code{\link[stats]{optim}}.
#' @param x independent variable. Should be continuous and numeric (e.g. day of year, hour of the day)
#' @param y dependent variable. Should be binary with values of 0 and 1.
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
#' plot(rownames(ptab),ptab[,"1"],xlab="doy",ylab="proportion with SI",pch=16)
#' abline(v=m.init+c(0,d.init),lty=2)
#' abline(h=a.init,lty=2)
#'
#' x_pr <- seq(min(x),max(x),by=1)
#' y_pr <- pla_msda(x=x_pr,m=m.init,s=s.init,d=d.init,a=a.init)
#'
#' # plot prediction based on initial guesses
#' points(x_pr,y_pr,type="l",col="green")
#'
#' # Calculate negative log-likelihood associated with the initial values
#' NLL_init <- pla_msda_NLL(par=c("m"=m.init,"s"=s.init,"d"=d.init,"a"=a.init),
#'                          x=data$doy,
#'                          y=data$SI)
#' NLL_init
#'
#' # fit model
#' fit <- fit_pla_msda(x=x,y=y,args_optim=list(par=c("m"=m.init, "s"=s.init, "d"=d.init, "a"=a.init)))
#'
#' # plot predictions
#' curve(expr=do.call(pla_msda,c(list(x=x),as.list(fit$par))),
#'       add=TRUE,lwd=2,lty=1,col="blue")
#' legend("topright",legend=c("observed","initial","predicted"),
#'        pch=c(16,NA,NA),lty=c(NA,1,1),lwd=c(0,1,2),col=c("black","green","blue"),bty="n")
#'
#' # Get NLL associated with fit
#' NLL_fit <- fit$value
#'
#' # How much did the NLL improve during the optimization?
#' NLL_init-NLL_fit
#'
#' # This should equal NLL_fit
#' pla_msda_NLL(par=c("m"=fit$par[["m"]],"s"=fit$par[["s"]],"d"=fit$par[["d"]],"a"=fit$par[["a"]]),
#'              x=data$doy,
#'              y=data$SI)
#'
#'
#' #### Likelihood profiles
#'
#' # Generate likelihood profile for m
#' Lprof <- data.frame(m=seq(fit$par[["m"]]*.9,fit$par[["m"]]*1.1,length=11),
#'                     NLL=NA)
#' for(i in seq_along(Lprof$m)){
#'   mi <- Lprof$m[i]
#'   Lprof[i,"NLL"] <- pla_msda_NLL(par=c("m"=mi,"s"=fit$par[["s"]],"d"=fit$par[["d"]],"a"=fit$par[["a"]]),
#'                                  x=data$doy,
#'                                  y=data$SI)
#' }
#' par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' with(Lprof,plot(m,NLL,type="o"))
#' points(fit$par[["m"]],NLL_fit,pch=16,cex=2)
#'
#' # Generate likelihood profile for s
#' Lprof <- data.frame(s=seq(fit$par[["s"]]*.5,fit$par[["s"]]*2,length=11),
#'                     NLL=NA)
#' for(i in seq_along(Lprof$s)){
#'   si <- Lprof$s[i]
#'   Lprof[i,"NLL"] <- pla_msda_NLL(par=c("m"=fit$par[["m"]],"s"=si,"d"=fit$par[["d"]],"a"=fit$par[["a"]]),
#'                                  x=data$doy,
#'                                  y=data$SI)
#' }
#' with(Lprof,plot(s,NLL,type="o"))
#' points(fit$par[["s"]],NLL_fit,pch=16,cex=2)
#'
#' # Generate likelihood profile for d
#' Lprof <- data.frame(d=seq(fit$par[["d"]]*.9,fit$par[["d"]]*1.1,length=11),
#'                     NLL=NA)
#' for(i in seq_along(Lprof$d)){
#'   di <- Lprof$d[i]
#'   Lprof[i,"NLL"] <- pla_msda_NLL(par=c("m"=fit$par[["m"]],"s"=fit$par[["s"]],"d"=di,"a"=fit$par[["a"]]),
#'                                  x=data$doy,
#'                                  y=data$SI)
#' }
#' with(Lprof,plot(d,NLL,type="o"))
#' points(fit$par[["d"]],NLL_fit,pch=16,cex=2)
#'
#' # Generate likelihood profile for a
#' Lprof <- data.frame(a=seq(fit$par[["a"]]*.9,fit$par[["a"]]*1.1,length=11),
#'                     NLL=NA)
#' for(i in seq_along(Lprof$a)){
#'   ai <- Lprof$a[i]
#'   Lprof[i,"NLL"] <- pla_msda_NLL(par=c("m"=fit$par[["m"]],"s"=fit$par[["s"]],"d"=fit$par[["d"]],"a"=ai),
#'                                  x=data$doy,
#'                                  y=data$SI)
#' }
#' with(Lprof,plot(a,NLL,type="o"))
#' points(fit$par[["a"]],NLL_fit,pch=16,cex=2)
#'


pla_msda_NLL <- function(params=c("m"=NA,"s"=NA,"d"=NA,"a"=NA),x,y){
  m <- params[1]
  s <- params[2]
  d <- params[3]
  a <- params[4]
  P <-  pla_msda(x=x,m=m,s=s,d=d,a=a)
  -sum(log(pmax(dbinom(y, size = 1, prob = P, log = FALSE),1e-15))) # The pmax helps avoid -Inf values and thus errors
}

