#' Diel sampling
#'
#' A function which randomly chooses hours from the first or second half of the day
#' @param freq frequency of sampling (e.g. every 1 hour, every 2 hours)
#' @param n number of days to draw hours for
#' @param h_sunrise hour of sunrise
#' @param h_sunset hour of sunset
#' @param  pattern sample pattern: "rand"= randomly chooses whether to sample during the AM or PM, "alt"= alternate periods (e.g. 1,2,1,2,1,2,..), "all"= sample every hour of the day, "day"=sample every hour during daylight hours (e.g. 0600-1700h), "day_rand"=sample only during daylight hours and randomly choose 12 sampling times (i.e. repeat some hours and skip some hours)
#' @returns Character vector of times of day in 24h format: 00:00 to 23:00
#' @author Nikolai Klibansky
#' @export
#' @examples
#' diel_sample(n=1,pattern="day")
#' diel_sample(n=1,pattern="rand")
#' diel_sample(n=1,pattern="alt")
#' diel_sample(n=1,pattern="day_rand")
#' diel_sample(n=1,freq=1,pattern="all")
#' diel_sample(n=1,freq=2,pattern="all")
#' diel_sample(n=1,freq=4,pattern="all")

diel_sample <- function(n,
                    freq=1,
                    h_sunrise=6,
                    h_sunset=17,
                    pattern="day") {
  fs <- c(TRUE,rep(FALSE,freq-1))  # Select times at given frequency
  dh <- list("1st"=c(0:11)[fs],
          "2nd"=c(12:23)[fs])
  if(pattern=="rand"){
    return(paste(unlist(dh[sample(1:2,n,replace=TRUE)]),"00",sep=":"))}
  if(pattern=="alt"){
    return(paste(unlist(dh[rep(1:2,length=n)]),"00",sep=":")[fs])}
  if(pattern=="all"){
    return(paste(rep(0:23,n),"00",sep=":")[fs])}
  if(pattern=="day"){
    return(paste(rep(h_sunrise:h_sunset,n),"00",sep=":")[fs])}
  if(pattern=="day_rand"){
    out <- c() #initialize
    # For each of n sampling days, randomly draw 12 values during daylight hours, sort by sampling day
    # and combine into one long vector.
    for(i in 1:n){out <- c(out,sort(x=sample(h_sunrise:h_sunset,size=12,replace=T)))}
    return(paste(out,"00",sep=":")[fs])
  }
}
