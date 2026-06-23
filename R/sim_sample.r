#' Simulate sampling of a population of fish
#'
#' Simulate sampling of a population of fish
#' !!! FUNCTION UNDER CONSTRUCTION!!!
#' @param dates.smp.n Number of dates on which to sample the population
#' @param Sn0.m Start of known spawning season: month
#' @param Sn0.d Start of known spawning season: day
#' @param Sn1.m End of known spawning season: month
#' @param Sn1.d End of known spawning season: day
#' @param CrSn0.m Start of cruise (sampling) season: month
#' @param CrSn0.d Start of cruise (sampling) season: day
#' @param CrSn1.m End of cruise (sampling) season: month
#' @param CrSn1.d End of cruise (sampling) season: day
#' @param Ll.n Number of longlines to pull at each sampling event. Multiplied by empirical mean cpue to allow simulated manipulation of catch rates.
#' @param Ll.cpue.mn Mean longline catch rate (both sexes) CPUE (number of fish per longline set; 100 hooks per 1 hour set)
#' @param Ll.cpue.sd Standard deviation of longline catch rate (both sexes) CPUE (number of fish per longline set [100 hooks, 1 hour set])
#' @param D.Ll.TL.mean Mean total length of fish in longline catches
#' @param D.Ll.TL.sd Standard deviation of fish in longline catches
#' @param n_Lbinpsc Number of females to select per size class per longline. Only used when sampling_sizeselect is 4 or 5.
#' @param sampling_annual Annual sampling type: 1=sample all year, 2=sample during known spawning months only, 3=during cruise season
#' @param sampling_pattern Temporal sampling pattern: 2=random, 3=standard survey, 4=directed research
#' @param sampling_sizeselect Size-based sampling type. Types 1 and 2 both determine the number of females to sample as the number of length bins in the model (n_Lbin). Type 1 selects a perfectly uniform size distribution (one female per length bin), while type 2 selects n_Lbin  lengths in a random uniform manner. Types 4 and 5 both determine the number of females to sample based on empirical catch data by making a random draw from a normal distribution with mean = Ll.cpue.mn*Ll.n and sd = Ll.cpue.sd. When type=4 the function samples randomly from a normal size-frequency distribution with mean = D.Ll.TL.mean and sd = D.Ll.TL.sd. If sampling_sizeselect==5, sample from a
# uniform size-frequency distribution.
#' @param weather_ratio unfishable days/fishable days  (i.e. 0 indicates all days are fishable)
#' @param args_diel_sample Arguments passed to \code{\link[FishRepro]{diel_sample}} function.
#' @author Nikolai Klibansky
#' @examples
#' \dontrun{
#' # Write example here

#'}
#'
sim_sample <- function(d_sim=365,
                       dates.smp.n=30,
                       Sn0.m=4,
                       Sn0.d=1,
                       Sn1.m=10,
                       Sn1.d=31,
                       CrSn0.m=5,
                       CrSn0.d=1,
                       CrSn1.m=9,
                       CrSn1.d=30,
                       Ll.n=4,
                       Ll.cpue.mn=0.698,
                       Ll.cpue.sd=0.537,
                       D.Ll.TL.mean=582.2029,
                       D.Ll.TL.sd=100.6291,
                       n_Lbinpsc=1,
                       sampling_annual=3,
                       sampling_pattern=4,
                       sampling_sizeselect=4,
                       weather_ratio=0.135,
                       args_diel_sample=list(freq=1,pattern="day")
                      ) {
  # Model structure
  h_sim <- 24*d_sim        # number of hours in simulation
  t_sim <- 1:h_sim   # Time points in simulation (in hours)
  POSIX.t <- as.POSIXct(x=(0:(h_sim-1))*3600,origin="1970-01-01",tz="UTC") # Time points in simulation (in POSIX date and hour)

  # SAMPLING MODEL----------------------------------------------------------------

  # THIS MAY NOT BE NEEDED (2026-06-15)
  # # If dates smp.n is given as a percentage, than calculate it from the
  # # number of days in known spawning season
  # if(is.character(dates.smp.n)){
  #   dates.smp.n=round((as.numeric(dates.smp.n)/100)*as.numeric(as.Date(paste("1970",paste(Sn1.m),paste(Sn1.d),sep="-"))-
  #                                                                as.Date(paste("1970",paste(Sn0.m),paste(Sn0.d),sep="-"))+1))
  # }

  # THIS MAY NOT BE NEEDED (2026-06-15)
  # # Time(s) of day to sample at
  # if(sampling.diel==1){tod.smp <- paste(rep(h.sampling.diel,dates.smp.n),"00",sep=":")}
  # if(sampling.diel==2){tod.smp <- paste(sample(x=0:23, size=dates.smp.n, replace=TRUE),"00",sep=":")}
  # if(sampling.diel==3){tod.smp <- paste(0:23,"00",sep=":")}
  # if(sampling.diel==4){
  #   tod.smp <- diel_sample
  # }

  # Range of simulation time to sample during
  if(sampling_annual==1){
    POSIX.Sn0.h00=POSIX.t[1]; POSIX.Sn1.h23=POSIX.t[hs]
    h_smp_min=1;         h_smp_max=hs}

  if(sampling_annual==2){
    # Calculate hs values of first and last days of known spawning months
    # Determine the POSIXct value for the first hour of the first month
    POSIX.Sn0.h00=as.POSIXct(paste("1970",paste(Sn0.m),paste(Sn0.d),"00",sep="-"),tz="UTC")
    h_smp_min=which(POSIX.t==POSIX.Sn0.h00)

    # Determine the POSIXct value for the last hour of the last month
    POSIX.Sn1.h00=as.POSIXct(paste("1970",paste(Sn1.m),paste(Sn1.d),"00",sep="-"),tz="UTC") #  first hour of that month
    POSIX.Sn1.h23=POSIX.Sn1.h00+(3600*23)                                                   #  last hour of that month
    h_smp_max=which(POSIX.t==POSIX.Sn1.h23)}

  if(sampling_annual==3){
    # Calculate hs values of first and last days of cruise season
    # Determine the POSIXct value for the first hour of the first month
    POSIX.Sn0.h00=as.POSIXct(paste("1970",paste(CrSn0.m),paste(CrSn0.d),"00",sep="-"),tz="UTC")
    h_smp_min=which(POSIX.t==POSIX.Sn0.h00)

    # Determine the POSIXct value for the last hour of the last month
    POSIX.Sn1.h00=as.POSIXct(paste("1970",paste(CrSn1.m),paste(CrSn1.d),"00",sep="-"),tz="UTC") #  first hour of that month
    POSIX.Sn1.h23=POSIX.Sn1.h00+(3600*23)                                                       #  last hour of that month
    h_smp_max=which(POSIX.t==POSIX.Sn1.h23)
  }

  # Calculate Date values of first and last dates in the first and last known spawning months, respectively
  date.Sn0=as.Date(POSIX.Sn0.h00) # First date (first hour)
  date.Sn1=as.Date(POSIX.Sn1.h23) # Last date (first hour)

  # Choose sample dates:
  # Randomly (random dates and times of day)
  if (sampling_pattern==2) {
    h_smp=floor(sort(sample(x=h_smp_min:h_smp_max, size=dates.smp.n, replace=FALSE))) # Hours of simulation to sample at
    h_smp_POSIX.t=POSIX.t[h_smp]    # POSIX date-times to sample at
  }

  # standard survey
  if (sampling_pattern==3) {
    source("survey.schedule.r") # Run standard survey Scheduling program
    h_smp=D.traps$Time.Haul.hs       # Hours of simulation to sample at
    h_smp_POSIX.t=POSIX.t[h_smp]    # POSIX date-times to sample at
  }

  # DIRECTED RESEARCH
  # Attempt to sample on evenly spaced dates within the sampling period
  if (sampling_pattern==4) {
    dates.smp=local({
      # Determine the dates in the sample period
      dates.smp.period=seq(date.Sn0,date.Sn1,1)
      # Determine the number of dates in the sample period
      dates.smp.period.n=length(dates.smp.period)
      # Generate weather vector for each date (TRUE=weather good enough to sample)
      weather=sample(x=c(FALSE,TRUE), size = dates.smp.period.n, replace = TRUE, prob = c(weather_ratio,1))
      # Calculate the approximate length of each interval (i.e. sub-period)
      dates.smp.int=floor(dates.smp.period.n/dates.smp.n)
      # Divide the sample period into sub-periods based on dates.smp.n
      dates.smp.period.sub=ceiling(1:dates.smp.period.n/dates.smp.int)
      dates.smp=c()
      for(i in 1:length(unique(dates.smp.period.sub))){
        index=which(dates.smp.period.sub==unique(dates.smp.period.sub)[i])
        good.days=dates.smp.period[index][weather[index]]
        if(length(good.days)>0){
          dates.smp=c(dates.smp,min(good.days))
        }
      }
      dates.smp=as.Date(dates.smp,origin=as.Date("1970-01-01"))
      dates.smp
    })

    h_smp_POSIX.t=local({n.hrs=length(do.call(diel_sample,c(list(n=1),args_diel_sample))) # n hours per day
    as.POSIXct(paste(rep(dates.smp,each=n.hrs),
                     do.call(diel_sample,c(list(n=dates.smp.n),args_diel_sample))
                     ,sep=" "),tz="UTC")
    })
    h_smp=which(POSIX.t%in%h_smp_POSIX.t) # Hours of simulation to sample at
  }


  # SIZE-SELECTIVE FISH SAMPLING
  # Define what function to use to conduct size-selective fish sampling

  # Go fishing and conduct fish sampling on each predetermined sampling date and time
  if(sampling_sizeselect==1){ # Perfectly even uniform sampling: no size selectivity
    L.f.sampling_sizeselect=function(L.f){L.f}
  }
  if(sampling_sizeselect==2){ # Random-uniform sampling (with replacement) by size
    L.f.sampling_sizeselect=function(L.f){sample(x=L.f,size=n_Lbin,replace=T)}
  }

  # Directed research
  if(sampling_sizeselect%in%c(4,5)){
    # Calculate the number of fish to sample based on empirical catch data
    # and, if sampling_sizeselect==4 sample uniformly from empirically based
    # size-frequency distribution. If sampling_sizeselect==5, sample from a
    # uniform size-frequency distribution.

    L.f.sampling_sizeselect <- function(L.f){
      ### Calculate number of females caught in each longline (n=1 longline)
      #   NOTE: the inclusion of sex ratio (sex.rat) in the calculation
      Ll.catch.f <- pmax(round(rnorm(n = 1, mean = Ll.cpue.mn * Ll.n, sd = Ll.cpue.sd) * (1-sex.rat)),0)
      ### Size frequency distribution
      #   Since size frequencies appear very normally distributed, generate
      # capture probabilities based on a normal distribution defined by
      # mean and sd from empirical data
      if(sampling_sizeselect==4){
        Ll.prob.f <- dnorm(L.f,mean=D.Ll.TL.mean,sd=D.Ll.TL.sd)
      }
      #   Uniform size distribution
      if(sampling_sizeselect==5){
        Ll.prob.f <- rep(1,length(L.f))
      }
      #   to generate size frequencies in each longline catch
      smp=sample(x=L.f, size=Ll.catch.f, prob=Ll.prob.f,replace=TRUE)
      tab.smp=table(smp)
      as.numeric(rep(names(tab.smp),pmin(as.numeric(tab.smp),n_Lbinpsc)))
    }
  }



  # Create data frame of spawning indicator data of all size classes at each sampling event
  D.smp.raw <-    data.frame(cbind(L.f,f.M, ov_SI_type[,h_smp]))
  D.smp.raw.SC <- data.frame(cbind(L.f,f.M, ov_SC[,h_smp]))
  D.smp.raw.HO <- data.frame(cbind(L.f,f.M, ov_HO[,h_smp]))
  D.smp.raw.PF <- data.frame(cbind(L.f,f.M, ov_PF[,h_smp]))
  D.smp.raw.fb <- data.frame(cbind(L.f,f.M, ov_fb[,h_smp]))

  # Select fish from each sampling event (time) then compile the data for the season
  # Initialize empty data frame (i.e. data sheet)
  D.smp=data.frame("SamplingEvent"=integer(0),"Date"=character(0),"Time"=character(0),
                   "hs"=integer(0),"L.f"=integer(0),"f"=integer(0), "f.M"=integer(0),"HO"=integer(0),"PF"=integer(0),
                   "fb"=integer(0))

#   # Sample fish within the catch for each sampling event i
#   for(i in 1:length(h_smp)){
#     hs.i=h_smp[i]
#
#     # Choose lengths of females to sample
#     L.f.i=L.f.sampling_sizeselect(L.f)
#
#     # Number of observations per sampling event
#     n.obs.i=length(L.f.i)
#
#     if(n.obs.i>0){ # If any females were caught in a sampling event, "work them up" otherwise proceed
#       # to the next sampling event
#       # Sample females of appropriate lengths
#       # Female maturity
#       f.M.i=D.smp.raw[match(L.f.i,D.smp.raw$L.f),"f.M"]
#       # Female spawning capability
#       SC.i=local({D.smp.raw.SC[match(L.f.i,D.smp.raw.SC$L.f),paste("t",h_smp[i],sep=".")]})
#       # Batch1 structure presence
#       HO.i=local({HO.i=D.smp.raw.HO[match(L.f.i,D.smp.raw.HO$L.f),paste("t",h_smp[i],sep=".")]
#       return(as.numeric(HO.i>0))})
#       PF.i=local({PF.i=D.smp.raw.PF[match(L.f.i,D.smp.raw.PF$L.f),paste("t",h_smp[i],sep=".")]
#       return(as.numeric(PF.i>0))})
#       # Batch1 batch fecundity
#       fb.i=local({fb.i=D.smp.raw.fb[match(L.f.i,D.smp.raw.fb$L.f),paste("t",h_smp[i],sep=".")]
#       return(fb.i)})
#
#       # Compile data from sampling event i
#       data.i=data.frame(SamplingEvent=rep(i,n.obs.i),
#                         Date=rep(as.Date(format(POSIX.t[hs.i],"%Y-%m-%d")),n.obs.i),
#                         Time=rep(format(POSIX.t[hs.i],"%H:%M"),n.obs.i),
#                         hs=rep(hs.i,n.obs.i),
#                         L.f=L.f.i,
#                         f=1,
#                         f.M=f.M.i,
#                         SC=SC.i,
#                         HO=HO.i,
#                         PF=PF.i,
#                         fb=fb.i)
#
#       # Add data from sampling event i to sampling data
#       # (NOTE: each row represents the sampling of a single female)
#       D.smp=rbind(D.smp,data.i)
#     }
#   }
#   rm(L.f.i)
#
#   my.set.seed(iter=iter.i) # Reset random seed
#
#   # Add column indicating if any spawning indicators were present (i.e. HO and/or PF)
#   D.smp$SI <- pmin(rowSums(D.smp[,c("HO","PF")]),1)
#
#   # Sampling warnings
#   any.HO=sum(D.smp$HO)>0 # Were any HO observed in samples?
#   any.PF=sum(D.smp$PF)>0 # Were any PF observed in samples?
return(h_smp_POSIX.t)
}
