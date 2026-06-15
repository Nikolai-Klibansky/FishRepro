#' Fit plateau model with three parameters
#'
#' Fit plateau model with three parameters \code{\link[FishRepro]{pla_msd}}, where residual error follows a negative binomial distribution.
#' @param x numeric vector. A continuous variable
#' @param y numeric vector. A binary variable (0 or 1)
#' @param m point of inflection for left logistic function
#' @param s scale parameter shared by both logistic functions
#' @param d difference between left and right points of inflection. When fit to binary event data, average duration of events.
#' @param parLo Lower bounds for parameters
#' @param palUp Upper bounds for parameters
#' @author Nikolai Klibansky
#' @note A plateau model is type of double logistic model where the slope of the first (left) logistic function is positive and the slope of the second (right) logistic function is negative, resulting in a domed shape. If the slopes are steep enough and the points of inflection are far enough apart, then the model will have a flat-topped section between the points of inflection.
#' @export
#' @examples
#'  # Setup for some pretty colors
#'  Tr <- 0.75 # Transparency
#'  cols <- c(HO = rgb(0.33,0.10,0.55,Tr),      # purple4
#'            POF = rgb(1.00,0.08,0.58,Tr)      # deeppink
#'  )
#'
#'  # Example 1: fit to HO presence
#'  data <- data_SI_h_bin
#'
#'  # HO
#'  hour <- data$Hour
#'  HO <- data$HO
#'
#'  tab <- table(hour,HO)
#'  ptab <- prop.table(tab,1)
#'
#'  par(mar=c(3,3,1,1),mgp=c(1,.2,0),tck=-0.01)
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],pch=16,
#'       xlab="hour",
#'       ylab="proportion with HO",
#'       col=cols["HO"],
#'  )
#'
#'  fit <- fit_pla_msd(x,HO)
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["HO"])
#'  legend("topright",legend=c("observed","predicted"),
#'         pch=c(16,NA),lwd=2,col=cols["HO"],bty="n")
#'
#'  # Example 2: fit to POF presence
#'  # POF
#'  # NOTE: when fitting to the POF presence data, we need to fit to shifted hour (Hour.shf)
#'  # so that the domed shaped of the distribution is apparent
#'  hour.shf <- data$Hour.shf
#'  POF <- data$POF
#'
#'  # Plot POF presence by (unshifted) hour and see what it looks like
#'  # If you try to fit the pla_msd model to this data it does not work well.
#'  tab <- table(hour,POF)
#'  ptab <- prop.table(tab,1)
#'
#'  par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1,.2,0),tck=-0.01)
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],
#'       xlab="hour",
#'       ylab="proportion with POF",
#'       main="poor fit when x=hour",
#'       col=cols["POF"],pch=16)
#'
#'  fit <- fit_pla_msd(hour,POF)
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["POF"])
#'  legend("topright",legend=c("observed","predicted"),
#'         pch=c(16,NA),lwd=2,col=cols["POF"],bty="n")
#'
#'  # Now plot and fit by shifted hour
#'  tab <- table(hour.shf,POF)
#'  ptab <- prop.table(tab,1)
#'
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],
#'       xlab="hour (shifted)",
#'       ylab="proportion with POF",
#'       main="good fit when x=Hour.shf",
#'       col=cols["POF"],
#'       pch=16)
#'
#'  fit <- fit_pla_msd(hour.shf,POF)
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["POF"])
#'  legend("topright",legend=c("observed","predicted"),
#'         pch=c(16,NA),lwd=2,col=cols["POF"],bty="n")
#'
#'  # Example 3: fit to HO presence supplying your own initial parameter values
#'
#'  # Try good initial values
#'  tab <- table(hour,HO)
#'  ptab <- prop.table(tab,1)
#'
#'  par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1,.2,0),tck=-0.01)
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],
#'       xlab="hour",
#'       ylab="proportion with HO",
#'       main="good initial values",
#'       col=cols["HO"],
#'       pch=16
#'  )
#'
#'  # Use your own initial parameter values when fitting, which will override
#'  # internally generated initial values
#'  m.init <- 6
#'  s.init <- 1
#'  d.init <- 12 # try HO duration of 12 hours
#'
#'  # Plot the shape of the model with your initial guesses
#'  curve(expr=pla_msd(x,m=m.init,s=s.init,d=d.init),
#'        add=TRUE,
#'        lwd=2,lty=2,col=cols["HO"])
#'
#'  # Now fit the model with your initial parameter values and plot the fit
#'  fit <- fit_pla_msd(hour,HO,args_optim = list(par=c(m=m.init, s=s.init, d=d.init)))
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["HO"])
#'
#'  legend("topright",legend=c("observed","initial predictions","fitted predictions"),
#'         pch=c(16,NA),lwd=2,lty=c(NA,2,1),col=cols["HO"],bty="n")
#'
#'  # Try some poor starting values
#'  # In this case, `fit_pla_msd` still finds the best fit, but if the starting values
#'  # are too extreme it will fail.
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],
#'       xlab="hour",
#'       ylab="proportion with HO",
#'       main="poor starting values",
#'       col=cols["HO"],
#'       pch=16
#'  )
#'
#'  # Try some poor initial parameter values
#'  m.init <- 2
#'  s.init <- 10
#'  d.init <- 24
#'
#'  # Plot the shape of the model with your initial guesses
#'  curve(expr=pla_msd(x,m=m.init,s=s.init,d=d.init),
#'        add=TRUE,
#'        lwd=2,lty=2,col=cols["HO"])
#'
#'  # Now fit the model with your initial parameter values and plot the fit
#'  fit <- fit_pla_msd(hour,HO,args_optim = list(par=c(m=m.init, s=s.init, d=d.init)))
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["HO"])
#'  legend("topright",legend=c("observed","initial predictions","fitted predictions"),
#'         pch=c(16,NA),lwd=2,lty=c(NA,2,1),col=cols["HO"],bty="n")
#'
#'
#'  # Example 4: Add confidence intervals to observed values
#'  # and plot HO and POF on the same plot
#'
#'  tab <- table(hour,HO)
#'  ptab <- prop.table(tab,1)
#'
#'  par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1,.2,0),tck=-0.01)
#'  plot(x=as.numeric(rownames(ptab)),
#'       y=ptab[,"1"],
#'       xlab="hour",
#'       ylab="proportion with spawning indicator",
#'       col=cols["HO"],
#'       pch=16,
#'       xlim=c(0,36)
#'  )
#'
#'  CI <- apply(tab,1,FUN=function(xi){prop.test(x=xi["1"],n=sum(xi))}$conf.int)
#'  arrows(x0=as.numeric(rownames(ptab)),
#'         y0=CI[1,],y1=CI[2,],
#'         code=3,length=0.05,angle=90,col=cols["HO"])
#'
#'  fit <- fit_pla_msd(x,HO)
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["HO"])
#'  legend("topleft",legend=c("observed HO","predicted HO"),
#'         pch=c(16,NA),lwd=2,col=cols["HO"],bty="n")
#'
#'  # POF
#'  tab <- table(hour.shf,POF)
#'  ptab <- prop.table(tab,1)
#'
#'  points(x=as.numeric(rownames(ptab)),
#'         y=ptab[,"1"],
#'         xlab="hour (shifted)",
#'         ylab="proportion with POF",
#'         col=cols["POF"],
#'         pch=16)
#'  CI <- apply(tab,1,FUN=function(xi){prop.test(x=xi["1"],n=sum(xi))}$conf.int)
#'  arrows(x0=as.numeric(rownames(ptab)),
#'         y0=CI[1,],y1=CI[2,],
#'         code=3,length=0.05,angle=90,col=cols["POF"])
#'
#'  fit <- fit_pla_msd(hour.shf,POF)
#'  curve(expr=do.call(pla_msd,c(list(x=x),as.list(fit$par))),
#'        add=TRUE,lwd=2,lty=1,col=cols["POF"])
#'  legend("topright",legend=c("observed POF","predicted POF"),
#'         pch=c(16,NA),lwd=2,col=cols["POF"],bty="n")
#'

fit_pla_msd <- function(x,
                        y,
                        args_optim=list()
){
  # set initial parameter values
    x2 <- x[as.logical(y)]
    mn_x <- mean(x2)
    sd_x <- sd(x2)
    m.init <- mn_x - sd_x
    s.init <- sqrt((3*(sd_x^2))/(pi^2)) # use relationship between scale parameter and sd of a logistic distribution
    d.init <- 2*sd_x

#  Setup arguments
args_optim_user <- args_optim
args_optim_default <- list(fn=pla_msd_NLL,
                           par=c(m=m.init, s=s.init, d=d.init),
                           x=x,
                           y=y,
                           control=list(parscale=c(m=m.init, s=s.init, d=d.init)),
                           method="BFGS",
                           hessian=TRUE
)
args_optim <- modifyList(args_optim_default,args_optim_user)

# Fit model
do.call(optim,args_optim)
}
