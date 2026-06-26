#' Simulate a spawning population of fish
#'
#' Simulate a spawning population of fish to generate fecundity data. NOTE FUNCTION UNDER CONSTRUCTION!!!
#' @param d_sim number of days in the simulation
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param Lmin minimum length of fish (same units as Linf)
#' @param Lbinw width of length bins (same units as Linf)
#' @param a50 Age at maturity using knife edge maturity function (analogous to age at 50 percent maturity from logistic function)
#' @param lw_a a parameter in length weight equation (W = aL^b)
#' @param lw_b b parameter in length weight equation (W = aL^b)
#' @param doy_SI_mn Mean appearance of spawning markers during the year (day of year)
#' @param doy_SI_sd StDev appearance of spawning markers during the year (day of year)
#' @param Ps_mn Mean (size-independent) spawning fraction
#' @param Ps_m spawning fraction~length parameter m (point of inflection)
#' @param Ps_s spawning fraction~length parameter s (scale parameter)
#' @param Ps_a spawning fraction~length parameter a (upper asymptote)
#' @param Ps_f spawning fraction~length parameter f ( lower asymptote (i.e. floor))
#' @param dd_Sn_mn Mean duration of spawning period (days)
#' @param h_HO_mn Mean time (hour of the day) that HO become evident
#' @param h_HO_sd StDev time (hour of the day) that HO become evident
#' @param h_spn_mn Mean time (hour of the day) that females spawn (POF become evident)
#' @param h_spn_sd StDev time (hour of the day) that females spawn (POF become evident)
#' @param dh_HO_mn mean duration of hydrated oocytes (in hours)
#' @param dh_PF_mn mean duration in hours of post ovulatory follicles (in hours)
#' @param dh_SI_mn mean duration of spawning indicators (in hours), used when spawning indicator is the union of two events (i.e. SI presence = presence of either HO or POF)
#' @param cf_SI spawning indicator correction factor: 24/dh_SI_mn
#' @param fb_a Batch fecundity-length intercept parameter
#' @param fb_b Batch fecundity-length coefficient parameter
#' @param fb_c Batch fecundity-length exponent parameter
#' @param fb_k Batch fecundity-length k parameter for negative binomial error distribution
#' @param fb_thresh Threshold length (TL mm), below which fecundity is fixed at the observed value
#' @param fb_min Fecundity (observed) for fish < threshold length
#' @param mon_Sn_beg Start of known spawning season: month
#' @param day_Sn_beg Start of known spawning season: day
#' @param mon_Sn_end End of known spawning season: month
#' @param day_Sn_end End of known spawning season: day
#' @param nm_fn_fb Name of function used to compute batch fecundity
#' @param first_batch_error Include an error distribution around the date that the first batch of the season is spawned?
#' @param diel_spawn_error  Include an error distribution around the time of day that batches are spawned?
#' @param POF_TempDep Should POF duration be calculated as a function of water temperature based on empirical data?
#' @param fb_size_error Should error be included around the batch fecundity-length function?
#' @param SI_size_error Should error be included around the spawning fraction-length function?
#' @param BI_size_independent Should batch interval be independent of body size?
#' @param covmat_lgs_msaf Covariance matrix for parameters in four parameter logistic model fit to spawning fraction-length data. Only used when SI_size_error=TRUE.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Simulate spawning with defaults
#' set.seed(23456)
#' out_spawn <- sim_spawn()
#'
#' # Set dh_HO and dh_PF to match values in \code{sim_spawn}
#' par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
#' plot_spawn(data_pop=  out_spawn$data_pop,
#'            time_unit = "month",
#'            dh_HO = 6,
#'            dh_PF = 8
#' )
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            time_unit = "hour",
#'            dh_HO = 6,
#'            dh_PF = 8,
#'            nm_trt_plot  = c("Bs","HO","PF"))
#'}
#'
sim_spawn <- function(d_sim=365,
                      Linf=911.36,
                      K=0.24,
                      t0=-0.33,
                      Lmin=200,
                      Lbinw=20,
                      a50=2,
                      lw_a=0.0000165,
                      lw_b=2.99,
                      doy_SI_mn=156,
                      doy_SI_sd=38.62560942,
                      Ps_mn=0.389112903,
                      Ps_m=580.1673777,
                      Ps_s=66.92958968,
                      Ps_a=0.606122263,
                      Ps_f=0.213780242,
                      dd_Sn_mn=78,
                      h_HO_mn=10.19994175,
                      h_HO_sd=1.770091189,
                      h_spn_mn=16.38637712,
                      h_spn_sd=4.863383191,
                      dh_HO_mn=6,
                      dh_PF_mn=8,
                      dh_SI_mn=14,
                      cf_SI=1.714285714,
                      fb_a=0.407975854,
                      fb_b=3.03E-08,
                      fb_c=4.774001527,
                      fb_k=2.618605045,
                      fb_thresh=375,
                      fb_min=58935,
                      mon_Sn_beg=4,
                      day_Sn_beg=1,
                      mon_Sn_end=10,
                      day_Sn_end=31,
                      nm_fn_fb="pow3",
                      first_batch_error=TRUE,
                      diel_spawn_error=TRUE,
                      POF_TempDep=FALSE,
                      fb_size_error=TRUE,
                      SI_size_error=TRUE,
                      covmat_lgs_msaf = covmat_SI_L,
                      BI_size_independent=FALSE
                      ) {
#### SETUP MODEL ####
# BASIC DEMOGRAPHIC PARAMETERS
  Lfmax <- floorf(Linf,Lbinw)         # Female length: maximum
  Lbin <- seq(Lmin,Lfmax, by=Lbinw)    # Female length bins

  n_Lbin <- length(Lbin)                                       # Number of female length bins
  ID <- sprintf("%02d", 1:n_Lbin)                              # Fish ID based on name of length bin (formerly size.class)
  a_f <- vb_age(L=Lbin,Linf=Linf,K=K,t0=t0)                    # Age of females as a VB function of L
  L50 <- vb_len(a=a50,Linf=Linf,K=K,t0=t0)                     # Female length: at maturity

  Lbin_imm <- Lbin[Lbin<L50]                          # Immature female length bins
  Lbin_mat <- Lbin[Lbin>=L50]                         # Mature female length bins
  fm <- c(rep(0,length(Lbin_imm)),rep(1,length(Lbin_mat))) # Female maturity status by length bin
  W.f <- lw_a*Lbin^lw_b                   # W in grams
  P.fm.L <- (Lbin>=L50)*1                 # Proportion of females that are mature by length bin
  P.fm.a <- P.fm.L                        # Proportion of females that are mature by age bin
  P.fm <- sum(P.fm.L)/n_Lbin              # Proportion of all females that are mature

  # Model structure
  h_sim <- 24*d_sim        # number of hours in simulation
  t_sim <- 1:h_sim   # time points in simulation (in hours)
  POSIX.t <- as.POSIXct(x=(0:(h_sim-1))*3600,origin="1970-01-01",tz="UTC") # time points in simulation (in POSIX date and hour)

  # Start and end of spawning season in POSIX-time (rounded to whole day units)
  # True
  POSIX_Sn_beg <- as.POSIXct(x=paste("1970",mon_Sn_beg,day_Sn_beg,sep="-"),tz="UTC") # Start
  POSIX_Sn_end <- as.POSIXct(x=paste("1970",mon_Sn_end,day_Sn_end,sep="-"),tz="UTC") # End

  # Start and end time of spawning season (midnight on first day) in hour of the year units
  t_Sn_beg <- (as.numeric(POSIX_Sn_beg)/3600)+1     # start
  t_Sn_end <- (as.numeric(POSIX_Sn_end)/3600)+1     # end

  # Initialize vectors
  ord.day <- floor((t_sim-1)/24) # Ordinal (Julian) day at each time point
  P.SC <- rep(x=0,times=h_sim)  # Proportion of females SC at each time point
  P.HO <- rep(x=0,times=h_sim)  # Proportion of females with HO at each time point
  P.PF <- rep(x=0,times=h_sim)  # Proportion of females with PF at each time point

  # Oocyte maturation timeline
  #   h_HO_mn=  h_HO_mn             # mean time (hour of the day) that HO become evident

  # mean duration in hours for..
  #      dh_HO_mn= dh_HO_mn    # HO
  #      dh_PF_mn= dh_PF_mn    # PF
  #      d_BI_mn= d_BI_mn.B0*Lbin^d_BI_mn.B1 # Batch interval (days between initiation of batches)

  # Deterministic batch fecundity function
  fn_fb <- get(nm_fn_fb)

  #-- Initialize data matrices to store ovary data.

  # Raw data matrices: females X time
  ov_empty <- matrix(data=0,
                     nrow=n_Lbin,
                     ncol=h_sim,
                     dimnames=list(females=ID,time=paste("t",1:h_sim,sep=".")))
  # Batch 1 data
  ov_SC <- ov_empty      # Batch1 spawning capability (1=yes, 0=no)
  ov_SI_age <- ov_empty  # Batch1 age
  ov_SI_type <- ov_empty # Batch1 structure type (1=HO, 2=PF)
  ov_HO <- ov_empty      # Batch1 HO age (in hours)
  ov_PF <- ov_empty      # Batch1 PF age (in hours)
  ov_fb <- ov_empty      # Batch1 batch fecundity (i.e. number of HO)
  ov_sp <- ov_empty      # Spawns (i.e. spawning events): record the moment a batch is spawned

  # Data events occurring in the population, where there is a row for each
  # female size group for each event

  # Hour at which HO appeared
  data_pop_HO=data.frame("POSIX.t"=as.POSIXct(0,origin="1970-01-01",tz="UTC"),
                      "Lbin"=as.integer(0),"HO"=as.integer(0))
  data_pop_HO=data_pop_HO[0,]    # Clears that first stupid row that I had to add to appease POSIXct

  # Hour at which a spawn occurred (i.e. POF appeared)
  data_pop_Bs=data.frame("POSIX.t"=as.POSIXct(0,origin="1970-01-01",tz="UTC"),
                      "Lbin"=as.integer(0),"Bs"=as.integer(0),
                      "BI"=as.integer(0), # Duration of batch interval following a spawning event
                      "fb"=as.numeric(0), # Batch fecundity
                      "d_PF"=as.integer(0)) # Actual post-ovulatory follicle duration
  data_pop_Bs=data_pop_Bs[0,]    # Clears that first stupid row that I had to add to appease POSIXct

  # Summary data vectors
  t_Pd_beg=rep(0,n_Lbin)     # First hour of spawning period
  t_Pd_end=rep(0,n_Lbin)     # Last hour of spawning period
  date_Pd_beg=rep(NA,n_Lbin)  # First date of spawning period
  date_Pd_end=rep(NA,n_Lbin)  # Last date of spawning period
  t_SC_beg=rep(0,n_Lbin)     # First hour of spawning capable period
  t_SC_end=rep(0,n_Lbin)     # Last hour of spawning capable period
  Ds =rep(0,n_Lbin)    # Duration of spawning period
  Ps=   rep(0,n_Lbin)  # Spawning fraction
  BI=   rep(0,n_Lbin)  # Batch interval
  nb=   rep(0,n_Lbin)  # Number of batches
  fb=   rep(0,n_Lbin)  # Batch fecundity
  fa=   rep(0,n_Lbin)  # Annual fecundity

  # Summary data frame
  data_pop_par=data.frame(cbind(a_f,Lbin,W.f,fm,date_Pd_beg,date_Pd_end,Ds,Ps,BI,nb,fb,fa))


##### SPAWNER MODEL ####

  #--Error density data frames (to keep out of for loop)
  D.t_Sn_beg <-    dtnorm2(mean=doy_SI_mn, sd=doy_SI_sd)               # time of first batch
  D.t.sp <- dtnorm2(mean=h_spn_mn, sd=h_spn_sd)                        # time of spawning (i.e. POF appearance)
  D.t.HO <-  data.frame(x=D.t.sp[,"x"]-round(dh_HO_mn),y=D.t.sp[,"y"]) # time of HO appearance

  #-- Calculate values for each mature female size-class at each time step
  # #@@@@@ BEGIN FEMALE SIZE CLASS LOOP @@@@@@@@@@@@@@@@@@@@@@@
  for(i in which(Lbin>=L50)) { # For each mature female size class calculate..
    # time at midnight nearest to spawning of the first batch (in simulation time)
    if(first_batch_error){t_sp1_i=sample(x=D.t_Sn_beg$x,size=1,replace=TRUE,prob=D.t_Sn_beg$y)*24 + 1
    }else{t_sp1_i=doy_SI_mn*24+1}
    # my.set.seed(iter=iter.i+1) # Reset random seed

    # time at midnight nearest to spawning of the last batch (in simulation time)
    t_Sn_end.i= t_sp1_i + dd_Sn_mn*24 + 1
    #@@@@@ BEGIN TIME STEP LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@
    j <- 0
    # While the time to initiate the next batch is less than the time of the end of the spawning season..
    while(t_sp1_i<t_Sn_end.i) {
    j <- j+1
      t_sp1_i=as.integer(t_sp1_i)                                                         # Convert t_sp1_i to integer
      if(diel_spawn_error){
        t_sp_ij <- t_sp1_i+sample(x=D.t.sp$x,size=1,replace=TRUE,prob=D.t.sp$y)  # time of simulation that a batch is spawned
      }else{
        t_sp_ij <- t_sp1_i+round(h_spn_mn)
      }

      t_HO <- t_sp_ij-dh_HO_mn            # time of simulation that a batch is initiated (i.e. HO appear)
      d.B1.HO <- dh_HO_mn                 # Duration of HO stage        (hours; fixed)
      d_PF <- if(POF_TempDep){
        # POF duration for current batch as a function of water temperature based on empirical water temp data
        D.Buoy$d.PF[match(floor((t_sp_ij-1)/24),D.Buoy$Date.j)] # hours; temperature dependent
      }else{dh_PF_mn}                                           # hours; fixed



      d.B1= d.B1.HO+d_PF               # Duration of HO+PF stages    (hours)

      # Batch fecundity of current batch
      # Mean batch fecundity based on body size
      if(nm_fn_fb=="pow3"){
        fb_mean <- fn_fb(x=Lbin[i], fb_a, fb_b, fb_c)
      }else{
        fb_mean <- fn_fb(x=Lbin[i], fb_a, fb_b)
      }

      fb_i <- fb_mean # NOTE: This used to be round(fb_mean), since you can't have a fraction
      # of an egg, but I removed the rounding because when you fit a model
      # to the simulated data, if it is rounded you get very slightly different
      # model parameters. This is probably not necessary, but I want to assure
      # that the perfect model is really perfect.
      # Add error if desired
      if(fb_size_error){
        fb_i= rnbinom(n=1,mu=fb_mean,size=fb_k)
      }

      # Duration of batch interval
      d.B1.BI= local({
        if(SI_size_error){
          ## Use bootstrap parameter values
          # Conduct Monte-Carlo draws from fitted parameters, incorporating covariance structure
          # Draw parameter values from multivariate normal distribution
          pars <- local({
            pars_mean <- c("m"=Ps_m,"s"=Ps_s,"a"=Ps_a,"f"=Ps_f)
            pars <- tmvtnorm::rtmvnorm(n=1,
                                       mean=  pars_mean,
                                       sigma= covmat_lgs_msaf,
                                       lower=c("m"=0,"s"=1,"a"=0,f=0),
                                       upper=c("m"=Inf,"s"=Inf,"a"=1,f=1))
            dimnames(pars)[[2]] <- names(pars_mean)
            as.data.frame(round(pars,4))
          })

        }else{pars=data.frame("m"=Ps_m,"s"=Ps_s,"a"=Ps_a,"f"=Ps_f)}
        P=lgs_msaf(m=pars$m, s=pars$s, a=pars$a, f=pars$f, x=Lbin[i]) # Calculate predicted probabilities at length
        if(BI_size_independent){  # Override P if batch interval is set to be size-independent
          P=Ps_mn}
        P=pmin(cf_SI*P,1) # Note that cf_SI is used to scale P based on empirical data, therefore
        # the durations that HO and POF are set to in the simulation do not
        # actually affect the setting of cf_SI
        I=round((1/P)*24) # Convert to hours and round to nearest hour
        return(I)
      })

      # Record structure age and type for each batch in the appropriate data matrices
      if(t_HO<h_sim){
        # Batch1 age
        b=local({a=c(t_HO:(t_HO+d.B1)); a[a<h_sim]})  # Vector indices restricted to simulation time
        ov_SI_age[i,b]= b-min(b); rm(b)

        # From initiation up to spawning...
        b=local({a=c(t_HO:(t_sp_ij-1)); a[a<h_sim]})    # Vector indices restricted to simulation time
        # Record HO presence
        ov_SI_type[i,b]=1                                  # Set structure type to 1 (HO)
        # Record HO age
        ov_HO[i,b]=b-min(b)+1
        # Record batch fecundity
        ov_fb[i,b]=fb_i                                  # Record batch fecundity from initiation up to spawning
        rm(b)

      }

      if(t_sp_ij<h_sim){
        # Record PF presence
        b=local({a=c(t_sp_ij:(t_sp_ij+d_PF-1)); a[a<h_sim]}) # Vector indices restricted to simulation time
        ov_SI_type[i,b]=2 ; rm(b)                            # Set structure type to 2 (PF) from spawning up to PFs disappear

        # Record PF age
        b=local({a=c(t_sp_ij:(t_sp_ij+d_PF-1)); a[a<h_sim]}) # Vector indices restricted to simulation time
        ov_PF[i,b]=b-min(b)+1 ; rm(b)                          # Record HO age from spawning up to PFs disappear

        # Record when a female spawns
        ov_sp[i,t_sp_ij]=1
      }

      # Add data to data_pop_HO as spawns occur
      data_pop_HO.row=nrow(data_pop_HO)+1                      # Determine the next row of data_pop_HO
      data_pop_HO[data_pop_HO.row,"POSIX.t"]=POSIX.t[t_HO]  # Add POSIX.t
      data_pop_HO[data_pop_HO.row,"Lbin"]=Lbin[i]                # Add Lbin
      data_pop_HO[data_pop_HO.row,"HO"]=1                      # Add HO

      # Add data to data_pop_Bs as spawns occur
      data_pop_Bs.row=nrow(data_pop_Bs)+1                      # Determine the next row of data_pop_Bs
      data_pop_Bs[data_pop_Bs.row,"POSIX.t"]=POSIX.t[t_sp_ij] # Add POSIX.t
      data_pop_Bs[data_pop_Bs.row,"Lbin"]=Lbin[i]                # Add Lbin
      data_pop_Bs[data_pop_Bs.row,"Bs"]=1                      # Add Bs
      data_pop_Bs[data_pop_Bs.row,"BI"]=d.B1.BI                # Add BI
      data_pop_Bs[data_pop_Bs.row,"fb"]=fb_i                   # Add fb
      data_pop_Bs[data_pop_Bs.row,"d_PF"]=d_PF           # Add d_PF

      # Recalculate t_sp1_i (a little clunky)
      # Determine the date that this recent spawn occurred on to restrict the model
      # from allowing females to spawn more than once on the same date
      t_sp1_i <- local({
        # Determine the calendar date of the current spawn
        a1=POSIX.t[t_sp_ij]
        a2=format(a1,"%Y-%m-%d")
        a3=as.POSIXct(paste(a2,"00:00",sep=" "),tz="UTC")
        # Calculate the calendar date of the next spawn
        b1=min(floorf(t_sp_ij+d.B1.BI,24)+1,h_sim) # Restrict date from going beyond the end of the simulation to avoid an error
        b2=POSIX.t[b1]    # Convert t_sp1_i to POSIX

        if(b2==a3){b2=b2+60*60*24} # If date of the next spawn is the same
        # date as the current spawn advance the date by 1 day
        b3=min(b2,POSIX.t[h_sim])     # Restrict date from going beyond the end of the simulation to avoid an error

        which(POSIX.t==b3) # Convert t_sp1_i back to h_sim format
      })
    } #@@@@@ END TIME STEP LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    data_pop_Bs$POSIX.t=as.POSIXct(data_pop_Bs$POSIX.t,origin="1970-01-01",tz="UTC")


    # First and last simulation hour of the spawning period
    t_Pd_beg[i]=min(which(ov_sp[i,]==1))
    t_Pd_end[i]=max(which(ov_sp[i,]==1))

    # First and last date of the spawning period
    date_Pd_beg[i]=as.Date(format(POSIX.t[t_Pd_beg[i]],"%Y-%m-%d"))
    date_Pd_end[i]=as.Date(format(POSIX.t[t_Pd_end[i]],"%Y-%m-%d"))

    # First and last simulation hour of the spawning capable period
    # Calculate first hour an HO is present, then determine the first hour (0:00h) of that day
    # (i.e. assume that the female has been capable of spawning the entire day that you first observe
    #  an HO. This is a simplifying assumption that should have minimal effect on calculations)
    t_SC_beg[i]=floorf(min(which(ov_HO[i,]>0)),24)+1
    # Calculate last hour a POF is present, then determine the last hour (23:00h) of that day
    # (i.e. assume that the female is capable of spawning the entire day that you last observe
    #  a POF. This is a simplifying assumption that should have minimal effect on calculations)
    t_SC_end[i]=ceilingf(max(which(ov_PF[i,]>0)),24)

    # Record spawning capability data
    ov_SC[i,c(t_SC_beg[i]:t_SC_end[i])]=1

    # Vector of batches by day of year
    Bs.d=as.numeric(tapply(X=ov_sp[i,],INDEX=ord.day,FUN=sum)) # Batches spawned per day
    # Endpoints of the true spawning period (NOTE: Estimates of spawning periods are based on HO or
    #   POF while these estimates are based on actual spawning events)
    Bs.d0=min(which(Bs.d==1))               # First day
    Bs.d1=max(which(Bs.d==1))               # Last day

    # NOTE: In this standard calculation of spawning fraction, there is a
    # slight upward bias.  This occurs because the spawning period begins and
    # ends with a spawning batch.  Thus there is one more batch per spawning
    # interval than we typically think of.  An unbiased calculation could be made
    # by ending the spawning period the day before the last batch was spawned:
    # sum(Bs.d[Bs.d0:(Bs.d1-1)])/length(Bs.d[Bs.d0:(Bs.d1-1)])
    # But for now, I'm going to keep the standard calculation, since that's what
    # everyone is familiar with.
    Ps[i]=sum(Bs.d[Bs.d0:Bs.d1])/(Bs.d1-Bs.d0+1)    # Spawning fraction
    Ds[i]=as.numeric(date_Pd_end[i]-date_Pd_beg[i])+1     # Spawning period duration
    nb[i]=sum(ov_sp[i,])                               # Number of batches
    fb[i]=exp(mean(log(data_pop_Bs[data_pop_Bs$Lbin==Lbin[i],"fb"]))) # (mean) batch fecundity logged and then back-transformed, since the error is lognormal
    fa[i]=sum(data_pop_Bs[data_pop_Bs$Lbin==Lbin[i],"fb"])  # Annual fecundity
  } #@@@@@ END FEMALE SIZE CLASS LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # Estimated batch fecundity for population
  mp.pop <- local({
    x <- data_pop_Bs$Lbin; y <- data_pop_Bs$fb
    fb_thresh <- min(x)

    if(nm_fn_fb=="lin"){
      mp.pop <- as.numeric(coef(lm(y~x)))
    }
    if(nm_fn_fb=="exp.lin"){
      mp.pop <- as.numeric(coef(lm(log(y)~x)))
    }
    if(nm_fn_fb=="pow.lin"){
      mp.pop <- as.numeric(coef(lm(log(y)~log(x))))
    }
    if(nm_fn_fb=="pow.nbinom"){
      if(fb_size_error){
        mp.pop <- as.numeric(fit_pow_nbinom(x,y)$par)
      }else{
        mp.pop <- as.numeric(fit.pow.norm(x,y)$par)
      }
    }
    if(nm_fn_fb=="pow3"){
      if(fb_size_error){
        mp.pop <- as.numeric(fit_pow3_nbinom(x,y)$par)
        fb_min <- round(pow3(x=fb_thresh,a=mp.pop[1],b=mp.pop[2],c=mp.pop[3]))
      }else{
        mp.pop <- as.numeric(fit_pow3_norm(x,y)$par)
        fb_min <- round(pow3(x=fb_thresh,a=mp.pop[1],b=mp.pop[2],c=mp.pop[3]))
      }
      mp.pop[4] <- fb_thresh
      mp.pop[5] <- fb_min
    }
    return(mp.pop)
  })

  # my.set.seed(iter=iter.i) # Reset random seed (2025-08-13 I think that this is being reset more than necessary within each iteration, but it's not a problem.)

  if(nm_fn_fb=="pow3"){
    fb <- c(rep(0,length(Lbin_imm)),fn_fb(Lbin_mat, mp.pop[1], mp.pop[2], mp.pop[3]))}else{
      fb <- c(rep(0,length(Lbin_imm)),fn_fb(Lbin_mat, mp.pop[1], mp.pop[2]))
    }

  # Sort data_pop_Bs by POSIX.t (i.e. put the spawn data in chronological order)
  data_pop_Bs <- data_pop_Bs[order(data_pop_Bs$POSIX.t),]
  # Add columns to data_pop_Bs
  data_pop_Bs$Date=as.Date(format(data_pop_Bs$POSIX.t,"%Y-%m-%d"))

  #-- Compile batch number data for the entire population during the simulation
  # First day of spawning periods
  data_pop_par[,"date_Pd_beg"]=as.Date(date_Pd_beg,origin=as.Date("1970-01-01"))
  # Last day of spawning periods
  data_pop_par[,"date_Pd_end"]=as.Date(date_Pd_end,origin=as.Date("1970-01-01"))
  # Spawning period duration
  data_pop_par[,"Ds"]=Ds
  # Spawning fraction
  data_pop_par[,"Ps"]=Ps
  # Batch interval
  data_pop_par[which(Ps>0),"BI"]=1/Ps[which(Ps>0)]
  # Number of batches
  data_pop_par[,"nb"]=nb
  # (mean) Batch fecundity
  data_pop_par[,"fb"]=fb
  # Annual fecundity
  data_pop_par[,"fa"]=fa

  # Build VERY LARGE data frame where there is a row for each female for each step in the simulation
  # Initialize data frame
  data_pop=data.frame("t"=rep(t_sim,each=n_Lbin),                                          # Simulation time
                   "date"=rep(as.Date(format(POSIX.t[t_sim],"%Y-%m-%d")),each=n_Lbin),  # date
                   "time"=rep(format(POSIX.t[t_sim],"%H:%M"),each=n_Lbin),              # time (24 hour)
                   "ID"=rep(ID,h_sim),                            # Size class
                   "Lbin"=rep(Lbin,h_sim),                                          # Female length
                   "f"=rep(1,h_sim),                                                # Female (0=no, 1=yes; used for counting females when summarizing data)
                   "fm"=rep(fm,h_sim))                                            # Female maturity (0=immature, 1=mature)
  # Add batches spawned
  data_pop <- add_col(df1=data_pop,
                    df2=reshape_lite(df=as.data.frame(ov_sp), pattern="t.",names=c("ID","t","Bs")),
                    nm_col1="ID",nm_col2="t",nm_x="Bs")
  # Add SC
  data_pop <- add_col(df1=data_pop,
                    df2=reshape_lite(df=as.data.frame(ov_SC),pattern="t.",names=c("ID","t","SC")),
                    nm_col1="ID",nm_col2="t",nm_x="SC")
  # Add HO
  data_pop <- add_col(df1=data_pop,
                    df2=reshape_lite(df=as.data.frame((ov_HO>0)*1),pattern="t.",names=c("ID","t","HO")),
                    nm_col1="ID",nm_col2="t",nm_x="HO")
  # Add PF
  data_pop <- add_col(df1=data_pop,
                    df2=reshape_lite(df=as.data.frame((ov_PF>0)*1),pattern="t.",names=c("ID","t","PF")),
                    nm_col1="ID",nm_col2="t",nm_x="PF")

  # Add column indicating if any spawning indicators were present (i.e. HO and/or PF)
  data_pop$SI <- pmin(rowSums(data_pop[,c("HO","PF")]),1)

  # Sort data
  o=order(data_pop$t,data_pop$ID)
  data_pop=data_pop[o,]


  #-- True prevalence of structures in population at each time step (used in plots)
  #-- Sum presence of each structure for each female for each time step
  for(i in 1:h_sim) {P.SC[i]=length(ov_SC[,i][ov_SC[,i]>0])/length(ov_SC[,i])}  # Proportion SC
  for(i in 1:h_sim) {P.HO[i]=length(ov_HO[,i][ov_HO[,i]>0])/length(ov_HO[,i])}  # Proportion with HO
  for(i in 1:h_sim) {P.PF[i]=length(ov_PF[,i][ov_PF[,i]>0])/length(ov_PF[,i])}  # Proportion with PF

  # True count count data for population at each hour of the simulation
  ov_HO.PA=(ov_HO>0)*(24/dh_HO_mn)  # Convert ov_HO to presence-absence scaled by HO duration
  ov_PF.PA=(ov_PF>0)*(24/dh_PF_mn)  # Convert ov_PF to presence-absence scaled by PF-duration

  data_ov <- array(c(ov_SI_type, ov_SC, ov_HO, ov_PF, ov_fb),
                   dim=c(dim(ov_empty),5),
                   dimnames = c(dimnames(ov_empty),list("trait"=c("SI_type","SC","HO","PF","fb")))
  )

  return(
    list(data_pop_par=data_pop_par,
         data_pop=data_pop,
         data_ov=data_ov
    )
  )
}
