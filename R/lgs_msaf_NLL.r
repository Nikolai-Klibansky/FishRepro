#' Objective function for logistic model with four parameters
#'
#' Objective function for \code{\link[FishRepro]{lgs_msaf}} logistic model with four parameters. Computes negative log likelihood value. Parameter values are input in a manner compatible with \code{\link[stats]{optim}}.
#' @param params Named vector of parameters m, s, a, and f in \code{\link[FishRepro]{lgs_msaf}} function. `params` must be in order: `m`, `s`, `a`, `f`
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Spawning activity-size----------------
#' # Generate raw version of data set
#' data_SI_L_raw <- data.frame(
#'   "L"=rep(data_SI_L$L,data_SI_L$N.f),
#'   "SI"=local({ SI <- c()           # Spawning indicators (presence/absence)
#'   for(row.i in 1:nrow(data_SI_L)){
#'     D <- data_SI_L[row.i,]
#'     x <- rep(c(0,1),c(D$N.f-D$N.f.SI,D$N.f.SI))
#'     SI <- c(SI,x)
#'   }
#'   return(SI)
#'   }))
#'
#'
#' data_SI_L_raw_tab <- prop.table(table(data_SI_L_raw$L,data_SI_L_raw$SI),1)
#'
#' # Guess initial parameter values
#' gs_SI_L <- list(m=median(data_SI_L_raw$L),
#'                 s=sd(data_SI_L_raw$L)*sqrt(3)/pi,
#'                 a=.99,
#'                 f=0.01,
#'                 x=seq(min(data_SI_L_raw$L),max(data_SI_L_raw$L),length=100))
#'
#'
#' lower <- c(min(data_SI_L_raw$L),-Inf,0,0)
#' upper <- c(max(data_SI_L_raw$L),Inf,1,1)
#' fit_SI_L <- optim(fn=lgs_msaf_NLL, par=c(m=gs_SI_L$m, s=gs_SI_L$s, a=gs_SI_L$a, f=gs_SI_L$f),
#'                 x=data_SI_L_raw$L, y=data_SI_L_raw$SI,
#'                 method="L-BFGS-B", lower=lower, upper=upper,
#'                 control=list(parscale=c(gs_SI_L$m, gs_SI_L$s, gs_SI_L$a, gs_SI_L$f)),
#'                 hessian=TRUE)
#'
#' covmat_SI_L <- solve(fit_SI_L$hessian) # solve Hessian matrix, which gives us the covariance matrix
#'
#' par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-.02,lend="butt")
#'
#' # Plot observed P.SI by length
#' plot(as.numeric(rownames(data_SI_L_raw_tab)),data_SI_L_raw_tab[,"1"],
#'      col="blue3", lwd=2,pch=16, type="p",
#'      xlab="Length",
#'      ylab="Proportion",
#'      ylim=c(0,1))
#'   # Plot guesses so you can see how much the function changed during the fitting process
#'   points(x=gs_SI_L$x,
#'          y=lgs_msaf(m=gs_SI_L$m, s=gs_SI_L$s, a=gs_SI_L$a, f=gs_SI_L$f, x=gs_SI_L$x),
#'          col="blue3",lty=2,type="l")
#' # Plot fitted value
#' points(x=gs_SI_L$x,
#'        y=lgs_msaf(m=fit_SI_L$par[["m"]], s=fit_SI_L$par[["s"]], a=fit_SI_L$par[["a"]], f=fit_SI_L$par[["f"]], x=gs_SI_L$x),
#'        lwd=3, type="l", col="blue3")


lgs_msaf_NLL <- function(params,x,y){
  m=params[1]; s=params[2]; a=params[3]; f=params[4]
  P=FishRepro::lgs_msaf(m=m,s=s,a=a,f=f,x=x)
  -sum(log(pmax(dbinom(y,size=1,prob=P),0.000001)))
}
