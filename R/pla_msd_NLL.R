#' Objective function for fitting plateau model with three parameters
#'
#' Objective function for fitting plateau model with three parameters \code{\link[FishRepro]{pla_msd}}.
#' @param x numeric vector
#' @param m point of inflection for left logistic function
#' @param s scale parameter shared by both logistic functions
#' @param d difference between left and right points of inflection. When fit to binary event data, average duration of events.
#' @param parLo Lower bounds for parameters
#' @param palUp Upper bounds for parameters
#' @author Nikolai Klibansky
#' @note A plateau model is type of double logistic model where the slope of the first (left) logistic function is positive and the slope of the second (right) logistic function is negative, resulting in a domed shape. If the slopes are steep enough and the points of inflection are far enough apart, then the model will have a flat-topped section between the points of inflection.
#' @export
#' @examples
#' # Example
#' dt <- FishRepro::data_SI_d_raw
#' dt$Start.date <- as.Date(dt$Start.date)
#' dt$doy <- as.numeric(format(dt$Start.date,"%j"))
#' pla_msd_NLL(params=c("m"=round(min(dt$doy)+diff(range(dt$doy))*0.25),
#'                      "s"=1,
#'                      "d"=diff(range(dt$doy))/2),
#'             x=dt$doy,
#'             y=dt$SI
#' )

pla_msd_NLL <- function(params,x,y,
                        parLo=list("m"=1,"s"=0,"d"=1),
                        parUp=list("m"=365,"s"=100,"d"=365)
){
  bs1 <- 0.001
  bs2 <- 0.001
  pen_big <- 1000 # large penalty to add when parameter values go out of bounds

  m <- params[1]; s <- params[2]; d <- params[3]
  P <-  pla_msd(m=m,s=s,d=d,x=x)
  P <- pmin(pmax(P,0),1) # keep P between zero and 1


  # Add penalties to NLL for each parameter to keep each away from bounds
  pen_m <- dbeta((m-parLo$m)/(parUp$m-parLo$m),shape1=bs1,shape2=bs2)
  pen_s <- dbeta((s-parLo$s)/(parUp$s-parLo$s),shape1=bs1,shape2=bs2)
  pen_d <- dbeta((d-parLo$d)/(parUp$d-parLo$d),shape1=bs1,shape2=bs2)

  # If parameters get outside of bounds, add a large penalty
  if(m<parLo$m|m>parUp$m){pen_m <- pen_big}
  if(s<parLo$s|s>parUp$s){pen_s <- pen_big}
  if(d<parLo$d|d>parUp$d){pen_d <- pen_big}

  #The pmax helps avoid -Inf values and thus errors
  -sum(log(pmax(dbinom(y, size = 1, prob = P, log = FALSE),1e-15))) +
    pen_m + pen_s + pen_d
}

