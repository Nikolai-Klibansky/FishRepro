#' Fit logistic model with four parameters to y~x data
#'
#' Fit logistic model \code{\link[FishRepro]{lgs_msaf}} with four parameters.
#' @param x numeric vector. independent variable
#' @param y numeric vector. dependent variable
#' @param args_optim Additional arguments to pass to \code{\link[stats]{optim}}
#' @author Nikolai Klibansky
#' @export
#' @examples
#'  ### Spawning activity-size----------------
#'  # Convert spawning indicator~length data set to raw form, where SI is binary (absence=0, presence=1)
#'  data_bin <- data.frame(
#'    "L"=rep(data_SI_L$L,data_SI_L$N.f),
#'    "SI"=local({ SI <- c()
#'    for(row.i in 1:nrow(data_SI_L)){
#'      D <- data_SI_L[row.i,]
#'      x <- rep(c(0,1),c(D$N.f-D$N.f.SI,D$N.f.SI))
#'      SI <- c(SI,x)
#'    }
#'    return(SI)
#'    }))
#'  data_bin_tab <- prop.table(table(data_bin$L,data_bin$SI),1)
#'
#'  # Guess initial parameter values
#'  gs.m <- median(data_bin$L)
#'  gs.s <- sd(data_bin$L)*sqrt(3)/pi
#'  gs.a <- .99
#'  gs.f <- 0.01
#'  x_pr <- seq(min(data_bin$L),max(data_bin$L),length=100)
#'  # Fit model
#'  fit <- fit_lgs_msaf(x=data_bin$L, y=data_bin$SI, args_optim=list(control=list(parscale=c(gs.m, gs.s, gs.a, gs.f))))
#'  y_pr <- do.call(lgs_msaf,c(list(x=x_pr),as.list(fit$par)))
#'  # Plot results
#'  par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-.02,lend="butt")
#'
#'  # Plot observed proportion of females with spawning indicators by length
#'  data_bin_tab <- prop.table(table(data_bin$L,data_bin$SI),1)
#'  plot(as.numeric(rownames(data_bin_tab)),data_bin_tab[,"1"],
#'       col="blue3", lwd=2,pch=16, type="p",
#'       xlab="Length",
#'       ylab="Proportion",
#'       ylim=c(0,1))
#'  # Plot guesses so you can see how much the function changed during the fitting process
#'  points(x=x_pr,
#'         y=lgs_msaf(m=gs.m, s=gs.s, a=gs.a, f=gs.f, x=x_pr),
#'         col="blue3",lty=2,type="l")
#'  # Plot fitted value
#'  points(x=x_pr,
#'         y=y_pr,
#'         lwd=3, type="l", col="blue3")
fit_lgs_msaf <- function(x,y,args_optim=list()){
  # Guess initial parameter values
  gs.m <- median(x)
  gs.s <- sd(x)*sqrt(3)/pi
  gs.a <- .99
  gs.f <- 0.01

  lower <- c(min(x),-Inf,0,0)
  upper <- c(max(x),Inf,1,1)

  # Setup arguments
  args_optim_user <- args_optim
  args_optim_default <- list(fn=lgs_msaf_NLL,
                             par=c(m=gs.m, s=gs.s, a=gs.a, f=gs.f),
                             x=x,
                             y=y,
                             method="L-BFGS-B", lower=lower, upper=upper,
                             control=list(parscale=c(gs.m, gs.s, gs.a, gs.f)),
                             hessian=TRUE
  )
  args_optim <- modifyList(args_optim_default,args_optim_user)

  # Fit model
  do.call(optim,args_optim)
}

