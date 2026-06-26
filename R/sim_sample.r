#' Simulate sampling of a population of fish
#'
#' Simulate sampling of a population of fish
#' !!! FUNCTION UNDER CONSTRUCTION!!!
#' @param sim_spawn_out Object returned by \code{\link[FishRepro]{sim_spawn}}
#' @param dates_smp_n Number of dates on which to sample the population
#' @param Sn0_m Start of known spawning season: month
#' @param Sn0_d Start of known spawning season: day
#' @param Sn1_m End of known spawning season: month
#' @param Sn1_d End of known spawning season: day
#' @param CrSn0_m Start of cruise (sampling) season: month
#' @param CrSn0_d Start of cruise (sampling) season: day
#' @param CrSn1_m End of cruise (sampling) season: month
#' @param CrSn1_d End of cruise (sampling) season: day
#' @param Ll_n Number of longlines to pull at each sampling event. Multiplied by empirical mean cpue to allow simulated manipulation of catch rates.
#' @param Ll_cpue_mn Mean longline catch rate (both sexes) CPUE (number of fish per longline set; 100 hooks per 1 hour set)
#' @param Ll_cpue_sd Standard deviation of longline catch rate (both sexes) CPUE (number of fish per longline set. 100 hooks, 1 hour set)
#' @param Ll_TL_mean Mean total length of fish in longline catches
#' @param Ll_TL_sd Standard deviation of fish in longline catches
#' @param sex_rat Sex ratio in population
#' @param n_Lbinpsc Number of females to select per size class per longline. Only used when sampling_sizeselect is 4 or 5.
#' @param sampling_annual Annual sampling type: 1=sample all year, 2=sample during known spawning months only, 3=during cruise season
#' @param sampling_pattern Temporal sampling pattern: 2=random, 3=standard survey, 4=directed research
#' @param sampling_sizeselect Size-based sampling type. Types 1 and 2 both determine the number of females to sample as the number of length bins in the model (n_Lbin). Type 1 selects a perfectly uniform size distribution (one female per length bin), while type 2 selects n_Lbin  lengths in a random uniform manner. Types 4 and 5 both determine the number of females to sample based on empirical catch data by making a random draw from a normal distribution with mean = Ll_cpue_mn*Ll_n and sd = Ll_cpue_sd. When type=4 the function samples randomly from a normal size-frequency distribution with mean = Ll_TL_mean and sd = Ll_TL_sd. If sampling_sizeselect==5, sample from a uniform size-frequency distribution.
#' @param weather_ratio unfishable days/fishable days  (i.e. 0 indicates all days are fishable)
#' @param args_diel_sample Arguments passed to \code{\link[FishRepro]{diel_sample}} function.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Simulate spawning with defaults
#' set.seed(23456)
#' out_spawn <- sim_spawn()
#'
#' # Simulate sampling with defaults
#' out_sample <- sim_sample(sim_spawn_out = out_spawn)
#'
#' # Set dh_HO and dh_PF to match values in \code{sim_spawn}
#' par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' plot_spawn(data_pop=  out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "month",
#'            dh_HO = 6,
#'            dh_PF = 8
#' )
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "week",
#'            dh_HO = 6,
#'            dh_PF = 8)
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "day",
#'            dh_HO = 6,
#'            dh_PF = 8)
#' par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "hour",
#'            dh_HO = 6,
#'            dh_PF = 8,
#'            nm_trt_plot  = c("Bs","HO","PF"))
#'
#' # Simulate sampling with some minor improvements to sampling
#' out_sample <- sim_sample(sim_spawn_out = out_spawn,
#'                          args_diel_sample=list(freq=1,pattern="all")
#' )
#'
#' # Set dh_HO and dh_PF to match values in \code{sim_spawn}
#' par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' plot_spawn(data_pop=  out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "month",
#'            dh_HO = 6,
#'            dh_PF = 8
#' )
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "week",
#'            dh_HO = 6,
#'            dh_PF = 8)
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "day",
#'            dh_HO = 6,
#'            dh_PF = 8)
#' par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp,
#'            time_unit = "hour",
#'            dh_HO = 6,
#'            dh_PF = 8,
#'            nm_trt_plot  = c("Bs","HO","PF"))
#'
#' }
#'

sim_sample <- function(sim_spawn_out=NULL,
                       d_sim=365,
                       dates_smp_n=30,
                       Sn0_m=4,
                       Sn0_d=1,
                       Sn1_m=10,
                       Sn1_d=31,
                       CrSn0_m=5,
                       CrSn0_d=1,
                       CrSn1_m=9,
                       CrSn1_d=30,
                       Ll_n=4,
                       Ll_cpue_mn=0.698,
                       Ll_cpue_sd=0.537,
                       Ll_TL_mean=582.2029,
                       Ll_TL_sd=100.6291,
                       sex_rat=0.5,
                       n_Lbinpsc=1,
                       sampling_annual=2,
                       sampling_pattern=4,
                       sampling_sizeselect=4,
                       weather_ratio=0.135,
                       args_diel_sample=list(freq=1,pattern="day")
) {
  sso <- sim_spawn_out

  # Model structure
  h_sim <- 24*d_sim        # number of hours in simulation
  t_sim <- 1:h_sim   # time points in simulation (in hours)
  POSIX.t <- as.POSIXct(x=(0:(h_sim-1))*3600,origin="1970-01-01",tz="UTC") # time points in simulation (in POSIX date and hour)

  Lbin <- sso$data_pop_par$Lbin
  fm <- sso$data_pop_par$fm


  # SAMPLING MODEL----------------------------------------------------------------
  # Range of simulation time to sample during
  if(sampling_annual==1){
    POSIX.Sn0_h00 <- POSIX.t[1]
    POSIX.Sn1_h23 <- POSIX.t[hs]
    h_smp_min <- 1
    h_smp_max <- hs
  }

  if(sampling_annual==2){
    # Calculate hs values of first and last days of known spawning months
    # Determine the POSIXct value for the first hour of the first month
    POSIX.Sn0_h00 <- as.POSIXct(paste("1970",paste(Sn0_m),paste(Sn0_d),"00",sep="-"),tz="UTC")
    h_smp_min <- which(POSIX.t==POSIX.Sn0_h00)

    # Determine the POSIXct value for the last hour of the last month
    POSIX.Sn1_h00 <- as.POSIXct(paste("1970",paste(Sn1_m),paste(Sn1_d),"00",sep="-"),tz="UTC") #  first hour of that month
    POSIX.Sn1_h23 <- POSIX.Sn1_h00+(3600*23)                                                   #  last hour of that month
    h_smp_max <- which(POSIX.t==POSIX.Sn1_h23)}

  if(sampling_annual==3){
    # Calculate hs values of first and last days of cruise season
    # Determine the POSIXct value for the first hour of the first month
    POSIX.Sn0_h00 <- as.POSIXct(paste("1970",paste(CrSn0_m),paste(CrSn0_d),"00",sep="-"),tz="UTC")
    h_smp_min <- which(POSIX.t==POSIX.Sn0_h00)

    # Determine the POSIXct value for the last hour of the last month
    POSIX.Sn1_h00 <- as.POSIXct(paste("1970",paste(CrSn1_m),paste(CrSn1_d),"00",sep="-"),tz="UTC") #  first hour of that month
    POSIX.Sn1_h23 <- POSIX.Sn1_h00+(3600*23)                                                       #  last hour of that month
    h_smp_max <- which(POSIX.t==POSIX.Sn1_h23)
  }

  # Calculate date values of first and last dates in the first and last known spawning months, respectively
  date.Sn0 <- as.Date(POSIX.Sn0_h00) # First date (first hour)
  date.Sn1 <- as.Date(POSIX.Sn1_h23) # Last date (first hour)

  # Choose sample dates:
  # Randomly (random dates and times of day)
  if (sampling_pattern==2) {
    h_smp <- floor(sort(sample(x=h_smp_min:h_smp_max, size=dates_smp_n, replace=FALSE))) # Hours of simulation to sample at
    h_smp_POSIX.t <- POSIX.t[h_smp]    # POSIX date-times to sample at
  }

  # standard survey
  if (sampling_pattern==3) {
    source("survey.schedule.r") # Run standard survey Scheduling program
    h_smp <- D.traps$Time.Haul.hs       # Hours of simulation to sample at
    h_smp_POSIX.t <- POSIX.t[h_smp]    # POSIX date-times to sample at
  }

  # DIRECTED RESEARCH
  # Attempt to sample on evenly spaced dates within the sampling period
  if (sampling_pattern==4) {
    dates_smp <- local({
      # Determine the dates in the sample period
      dates_smp_period <- seq(date.Sn0,date.Sn1,1)
      # Determine the number of dates in the sample period
      dates_smp_period.n <- length(dates_smp_period)
      # Generate weather vector for each date (TRUE=weather good enough to sample)
      weather <- sample(x=c(FALSE,TRUE), size = dates_smp_period.n, replace = TRUE, prob = c(weather_ratio,1))
      # Calculate the approximate length of each interval (i.e. sub-period)
      dates_smp_int <- floor(dates_smp_period.n/dates_smp_n)
      # Divide the sample period into sub-periods based on dates_smp_n
      dates_smp_period.sub <- ceiling(1:dates_smp_period.n/dates_smp_int)
      dates_smp <- c()
      for(i in seq_along(unique(dates_smp_period.sub))){
        index <- which(dates_smp_period.sub==unique(dates_smp_period.sub)[i])
        good.days <- dates_smp_period[index][weather[index]]
        if(length(good.days)>0){
          dates_smp <- c(dates_smp,min(good.days))
        }
      }
      dates_smp <- as.Date(dates_smp,origin=as.Date("1970-01-01"))
      dates_smp
    })

    h_smp_POSIX.t <- local({n.hrs <- length(do.call(diel_sample,c(list(n=1),args_diel_sample))) # n hours per day
    as.POSIXct(paste(rep(dates_smp,each=n.hrs),
                     do.call(diel_sample,c(list(n=dates_smp_n),args_diel_sample))
                     ,sep=" "),tz="UTC")
    })
    h_smp <- which(POSIX.t%in%h_smp_POSIX.t) # Hours of simulation to sample at
  }


  # SIZE-SELECTIVE FISH SAMPLING
  # Define what function to use to conduct size-selective fish sampling

  # Go fishing and conduct fish sampling on each predetermined sampling date and time
  if(sampling_sizeselect==1){ # Perfectly even uniform sampling: no size selectivity
    Lbin_sampling_sizeselect <- function(Lbin){Lbin}
  }
  if(sampling_sizeselect==2){ # Random-uniform sampling (with replacement) by size
    Lbin_sampling_sizeselect <- function(Lbin){
      n_Lbin <- length(Lbin)
      sample(x=Lbin,size=n_Lbin,replace=T)
    }
  }

  # Directed research
  if(sampling_sizeselect%in%c(4,5)){
    # Calculate the number of fish to sample based on empirical catch data
    # and, if sampling_sizeselect==4 sample uniformly from empirically based
    # size-frequency distribution. If sampling_sizeselect==5, sample from a
    # uniform size-frequency distribution.

    Lbin_sampling_sizeselect <- function(Lbin){
      ### Calculate number of females caught in each longline (n=1 longline)
      #   NOTE: the inclusion of sex ratio (sex_rat) in the calculation
      Ll_catch_f <- pmax(round(rnorm(n = 1, mean = Ll_cpue_mn * Ll_n, sd = Ll_cpue_sd) * (1-sex_rat)),0)
      ### Size frequency distribution
      #   Since size frequencies appear very normally distributed, generate
      # capture probabilities based on a normal distribution defined by
      # mean and sd from empirical data
      if(sampling_sizeselect==4){
        Ll_prob_f <- dnorm(Lbin,mean=Ll_TL_mean,sd=Ll_TL_sd)
      }
      #   Uniform size distribution
      if(sampling_sizeselect==5){
        Ll_prob_f <- rep(1,length(Lbin))
      }
      #   to generate size frequencies in each longline catch
      smp <- sample(x=Lbin, size=Ll_catch_f, prob=Ll_prob_f,replace=TRUE)
      tab_smp <- table(smp)
      as.numeric(rep(names(tab_smp),pmin(as.numeric(tab_smp),n_Lbinpsc)))
    }
  }



  # Create data frame of spawning indicator data of all size classes at each sampling event
  data_smp_raw    <- data.frame(cbind(Lbin,fm, sso$data_ov[,h_smp,"SI_type"]))
  data_smp_raw_SC <- data.frame(cbind(Lbin,fm, sso$data_ov[,h_smp,"SC"]))
  data_smp_raw_HO <- data.frame(cbind(Lbin,fm, sso$data_ov[,h_smp,"HO"]))
  data_smp_raw_PF <- data.frame(cbind(Lbin,fm, sso$data_ov[,h_smp,"PF"]))
  data_smp_raw_fb <- data.frame(cbind(Lbin,fm, sso$data_ov[,h_smp,"fb"]))

  # Select fish from each sampling event (time) then compile the data for the season
  # Initialize empty data frame (i.e. data sheet)
  data_smp <- data.frame("SamplingEvent"=integer(0),"date"=character(0),"time"=character(0),
                         "hs"=integer(0),"Lbin"=integer(0),"f"=integer(0), "fm"=integer(0),
                         "HO"=integer(0),"PF"=integer(0),
                         "fb"=integer(0))

  # Sample fish within the catch for each sampling event i
  for(i in seq_along(h_smp)){
    hs_i <- h_smp[i]

    # Choose lengths of females to sample
    Lbin_i <- Lbin_sampling_sizeselect(Lbin)

    # Number of observations per sampling event
    n_obs_i <- length(Lbin_i)

    if(n_obs_i>0){ # If any females were caught in a sampling event, "work them up" otherwise proceed
      # to the next sampling event
      # Sample females of appropriate lengths
      # Female maturity
      fm_i <- data_smp_raw[match(Lbin_i,data_smp_raw$Lbin),"fm"]
      # Female spawning capability
      SC_i <- local({data_smp_raw_SC[match(Lbin_i,data_smp_raw_SC$Lbin),paste("t",h_smp[i],sep=".")]})
      # Batch1 structure presence
      HO_i <- local({HO_i <- data_smp_raw_HO[match(Lbin_i,data_smp_raw_HO$Lbin),paste("t",h_smp[i],sep=".")]
      return(as.numeric(HO_i>0))})
      PF_i <- local({PF_i <- data_smp_raw_PF[match(Lbin_i,data_smp_raw_PF$Lbin),paste("t",h_smp[i],sep=".")]
      return(as.numeric(PF_i>0))})
      # Batch1 batch fecundity
      fb_i <- local({fb_i <- data_smp_raw_fb[match(Lbin_i,data_smp_raw_fb$Lbin),paste("t",h_smp[i],sep=".")]
      return(fb_i)})

      # Compile data from sampling event i
      data_i <- data.frame(SamplingEvent=rep(i,n_obs_i),
                           date=rep(as.Date(format(POSIX.t[hs_i],"%Y-%m-%d")),n_obs_i),
                           time=rep(format(POSIX.t[hs_i],"%H:%M"),n_obs_i),
                           hs=rep(hs_i,n_obs_i),
                           Lbin=Lbin_i,
                           f=1,
                           fm=fm_i,
                           SC=SC_i,
                           HO=HO_i,
                           PF=PF_i,
                           fb=fb_i)

      # Add data from sampling event i to sampling data
      # (NOTE: each row represents the sampling of a single female)
      data_smp <- rbind(data_smp,data_i)
    }
  }

  # my.set.seed(iter <- iter_i) # Reset random seed

  # Add column indicating if any spawning indicators were present (i.e. HO and/or PF)
  data_smp$SI <- pmin(rowSums(data_smp[,c("HO","PF")]),1)

  # Sampling warnings
  any_HO <- sum(data_smp$HO)>0 # Were any HO observed in samples?
  any_PF <- sum(data_smp$PF)>0 # Were any PF observed in samples?
  return(list(data_smp=data_smp))
}
