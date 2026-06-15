#' Simulate a spawning population of fish
#'
#' Simulate a spawning population of fish to generate fecundity data. NOTE FUNCTION UNDER CONSTRUCTION!!!
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param Lmin minimum length of fish (same units as Linf)
#' @param Lbinw width of length bins (same units as Linf)
#' @param a50 Age at maturity using knife edge maturity function (analogous to age at 50 percent maturity from logistic function)
#' @param a_lw a parameter in length weight equation (W = aL^b)
#' @param b_lw b parameter in length weight equation (W = aL^b)
#' @param d_sim number of days in the simulation
#' @param j.mn.SI Mean appearance of spawning markers during the year (Julian day)
#' @param j.sd.SI StDev appearance of spawning markers during the year (Julian day)
#' @param Ps.mn Mean (size-independent) spawning fraction
#' @param Ps.m spawning fraction~length parameter m (point of inflection)
#' @param Ps.s spawning fraction~length parameter s (scale parameter)
#' @param Ps.a spawning fraction~length parameter a (upper asymptote)
#' @param Ps.f spawning fraction~length parameter f ( lower asymptote (i.e. floor))
#' @param d.mn.Sn Mean duration of spawning period (days)
#' @param h.mn.HO Mean time (hour of the day) that HO become evident
#' @param h.sd.HO StDev time (hour of the day) that HO become evident
#' @param h.mn.spn Mean time (hour of the day) that females spawn (POF become evident)
#' @param h.sd.spn StDev time (hour of the day) that females spawn (POF become evident)
#' @param d.mn.HO mean duration of hydrated oocytes (in hours)
#' @param d.mn.PF mean duration in hours of post ovulatory follicles (in hours)
#' @param SI.dur mean duration of spawning indicators (in hours), used when spawning indicator is the union of two events (i.e. SI presence = presence of either HO or POF)
#' @param SI_cf spawning indicator correction factor: 24/SI.dur
#' @param fb.a Batch fecundity-length intercept parameter
#' @param fb.b Batch fecundity-length coefficient parameter
#' @param fb.c Batch fecundity-length exponent parameter
#' @param fb.k Batch fecundity-length k parameter for negative binomial error distribution
#' @param fb.thresh Threshold length (TL mm), below which fecundity is fixed at the observed value
#' @param fb.min Fecundity (observed) for fish < threshold length
#' @param Sn0.m Start of known spawning season: month
#' @param Sn0.d Start of known spawning season: day
#' @param Sn1.m End of known spawning season: month
#' @param Sn1.d End of known spawning season: day
#' @param nm_fn_fb Name of function used to compute batch fecundity
#' @param first_batch_error Include an error distribution around the date that the first batch of the season is spawned?
#' @param diel_spawn_error  Include an error distribution around the time of day that batches are spawned?
#' @param POF_TempDep Should POF duration be calculated as a function of water temperature based on empirical data?
#' @param fb_size_error Should error be included around the batch fecundity-length function?
#' @param SI_size_error Should error be included around the spawning fraction-length function?
#' @param BI_size_independent Should batch interval be independent of body size?
#' @param covmat_lgs_msaf Covariance matrix for parameters in four parameter logistic model fit to spawning fraction-length data. Only used when SI_size_error=TRUE.
#'
#' @author Nikolai Klibansky
#' @keywords internal
#' @examples
#' \dontrun{
#' # Write example here

#'}
#'
sim_spawn <- function(Linf=911.36,
                      K=0.24,
                      t0=-0.33,
                      Lmin=200,
                      Lbinw=20,
                      a50=2,
                      a_lw=0.0000165,
                      b_lw=2.99,
                      d_sim=365,
                      j.mn.SI=156,
                      j.sd.SI=38.62560942,
                      Ps.mn=0.389112903,
                      Ps.m=580.1673777,
                      Ps.s=66.92958968,
                      Ps.a=0.606122263,
                      Ps.f=0.213780242,
                      d.mn.Sn=78,
                      h.mn.HO=10.19994175,
                      h.sd.HO=1.770091189,
                      h.mn.spn=16.38637712,
                      h.sd.spn=4.863383191,
                      d.mn.HO=6,
                      d.mn.PF=8,
                      SI.dur=14,
                      SI_cf=1.714285714,
                      fb.a=0.407975854,
                      fb.b=3.03E-08,
                      fb.c=4.774001527,
                      fb.k=2.618605045,
                      fb.thresh=375,
                      fb.min=58935,
                      Sn0.m=4,
                      Sn0.d=1,
                      Sn1.m=10,
                      Sn1.d=31,
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

  n_Lbin <- length(Lbin)                                         # Number of female length bins
  size.class <- paste("f",sprintf("%02d", 1:n_Lbin),sep=".")  # Names of length bins
  a_f <- vb_age(L=Lbin,Linf=Linf,K=K,t0=t0)                      # Age of females as a VB function of L
  L50 <- vb_len(a=a50,Linf=Linf,K=K,t0=t0)                     # Female length: at maturity

  Lbin_imm <- Lbin[Lbin<L50]                          # Immature female length bins
  Lbin_mat <- Lbin[Lbin>=L50]                         # Mature female length bins
  fm <- c(rep(0,length(Lbin_imm)),rep(1,length(Lbin_mat))) # Female maturity status by length bin
  W.f <- a_lw*Lbin^b_lw                   # W in grams
  P.fm.L <- (Lbin>=L50)*1                 # Proportion of females that are mature by length bin
  P.fm.a <- P.fm.L                        # Proportion of females that are mature by age bin
  P.fm <- sum(P.fm.L)/n_Lbin              # Proportion of all females that are mature

  # Model structure
  h_sim <- 24*d_sim        # number of hours in simulation

  # Start and end of spawning season in POSIX-time (rounded to whole day units)
  #   (NOTE: The endpoints of the spawning season here are effectively modeled as
  #    environmental traits.  In the spawner model, fish respond to these traits
  #    and begin or end spawning activity for the year.)
  # Expected (i.e. "known spawning season")
  POSIX.Sn0.kss <- as.POSIXct(x=paste("1970",Sn0.m,Sn0.d,sep="-"),tz="UTC")  # Start
  POSIX.Sn1.kss <- as.POSIXct(x=paste("1970",Sn1.m,Sn1.d,sep="-"),tz="UTC")  # End

  # True
  POSIX.Sn0 <- as.POSIXct(x=paste("1970",Sn0.m,Sn0.d,sep="-"),tz="UTC") # Start
  POSIX.Sn1 <- as.POSIXct(x=paste("1970",Sn1.m,Sn1.d,sep="-"),tz="UTC") # End

  # Start and end time of spawning season (midnight on first day) in hour of the year units
  t.Sn0 <- (as.numeric(POSIX.Sn0)/3600)+1     # start
  t.Sn1 <- (as.numeric(POSIX.Sn1)/3600)+1     # end

  # Initialize vectors
  t_sim <- 1:h_sim   # Time points in simulation (in hours)
  POSIX.t <- as.POSIXct(x=(0:(h_sim-1))*3600,origin="1970-01-01",tz="UTC") # Time points in simulation (in POSIX date and hour)

  ord.day <- floor((t_sim-1)/24) # Ordinal (Julian) day at each time point
  P.SC <- rep(x=0,times=h_sim)  # Proportion of females SC at each time point
  P.HO <- rep(x=0,times=h_sim)  # Proportion of females with HO at each time point
  P.PF <- rep(x=0,times=h_sim)  # Proportion of females with PF at each time point

  # Oocyte maturation timeline
  #   h.mn.HO=  h.mn.HO             # mean time (hour of the day) that HO become evident

  # mean duration in hours for..
  #      d.mn.HO= d.mn.HO    # HO
  #      d.mn.PF= d.mn.PF    # PF
  #      d.mn.BI= d.mn.BI.B0*Lbin^d.mn.BI.B1 # Batch interval (days between initiation of batches)

  # Deterministic batch fecundity function
  fn_fb <- get(nm_fn_fb)

  #-- Initialize data storage units:
  # Raw data matrices: females X time
  B1.blank <- matrix(data=0,nrow=n_Lbin,ncol=h_sim,
                  dimnames=list(females=size.class,time=paste("t",1:h_sim,sep=".")))
  # Batch 1 data
  B1.a <- B1.blank  # Batch1 age
  B1.o <- B1.blank  # Batch1 structure type (1=HO, 2=PF)
  B1.HO <- B1.blank # Batch1 HO age (in hours)
  B1.PF <- B1.blank # Batch1 PF age (in hours)
  B1.SC <- B1.blank # Batch1 spawning capability (1=yes, 0=no)
  B1.fb <- B1.blank # Batch1 batch fecundity (i.e. number of HO)
  Bs <- B1.blank    # Spawns (i.e. spawning events): record the moment a batch is spawned

  # Data events occurring in the population, where there is a row for each
  # female size group for each event

  # Hour at which HO appeared
  D.pop.HO=data.frame("POSIX.t"=as.POSIXct(0,origin="1970-01-01",tz="UTC"),
                      "Lbin"=as.integer(0),"HO"=as.integer(0))
  D.pop.HO=D.pop.HO[0,]    # Clears that first stupid row that I had to add to appease POSIXct

  # Hour at which a spawn occurred (i.e. POF appeared)
  D.pop.Bs=data.frame("POSIX.t"=as.POSIXct(0,origin="1970-01-01",tz="UTC"),
                      "Lbin"=as.integer(0),"Bs"=as.integer(0),
                      "BI"=as.integer(0), # Duration of batch interval following a spawning event
                      "fb"=as.numeric(0), # Batch fecundity
                      "d.B1.PF"=as.integer(0)) # Actual post-ovulatory follicle duration
  D.pop.Bs=D.pop.Bs[0,]    # Clears that first stupid row that I had to add to appease POSIXct

  # Summary data vectors
  t.Pd0=rep(0,n_Lbin)     # First hour of spawning period
  t.Pd1=rep(0,n_Lbin)     # Last hour of spawning period
  date.Pd0=rep(NA,n_Lbin)  # First date of spawning period
  date.Pd1=rep(NA,n_Lbin)  # Last date of spawning period
  t.SC0=rep(0,n_Lbin)     # First hour of spawning capable period
  t.SC1=rep(0,n_Lbin)     # Last hour of spawning capable period
  Ds =rep(0,n_Lbin)    # Duration of spawning period
  Ps=   rep(0,n_Lbin)  # Spawning fraction
  BI=   rep(0,n_Lbin)  # Batch interval
  nb=   rep(0,n_Lbin)  # Number of batches
  fb=   rep(0,n_Lbin)  # Batch fecundity
  fa=   rep(0,n_Lbin)  # Annual fecundity

  # Summary data frame
  D.pop.par=data.frame(cbind(a_f,Lbin,W.f,fm,date.Pd0,date.Pd1,Ds,Ps,BI,nb,fb,fa))


##### SPAWNER MODEL ####

  #--Error density data frames (to keep out of for loop)
  D.t.Sn0 <-    dtnorm2(mean=j.mn.SI, sd=j.sd.SI)                                # Time of first batch
  D.t.B1.spn <- dtnorm2(mean=h.mn.spn, sd=h.sd.spn)                              # Time of spawning (i.e. POF appearance)
  D.t.B1.HO <-  data.frame(x=D.t.B1.spn[,"x"]-round(d.mn.HO),y=D.t.B1.spn[,"y"]) # Time of HO appearance

  #-- Calculate values for each mature female size-class at each time step
  # #@@@@@ BEGIN FEMALE SIZE CLASS LOOP @@@@@@@@@@@@@@@@@@@@@@@
  for(i in which(Lbin>=L50)) { # For each mature female size class calculate..
    # Time at midnight nearest to spawning of the first batch (in simulation time)
    if(first_batch_error){t.B1=sample(x=D.t.Sn0$x,size=1,replace=TRUE,prob=D.t.Sn0$y)*24 + 1
    }else{t.B1=j.mn.SI*24+1}
    # my.set.seed(iter=iter.i+1) # Reset random seed

    # Time at midnight nearest to spawning of the last batch (in simulation time)
    t.Sn1.i= t.B1 + d.mn.Sn*24 + 1

    #@@@@@ BEGIN TIME STEP LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@
    while(t.B1<t.Sn1.i) {  # While the time to initiate the next batch is less
      # than the time of the end of the spawning season..
      t.B1=as.integer(t.B1)                                                         # Convert t.B1 to integer
      if(diel_spawn_error){
        t.B1.spn=t.B1+sample(x=D.t.B1.spn$x,size=1,replace=TRUE,prob=D.t.B1.spn$y)  # Time of simulation that a batch is spawned
      }else{
        t.B1.spn=t.B1+round(h.mn.spn)
      }



      t.B1.HO=t.B1.spn-d.mn.HO            # Time of simulation that a batch is initiated (i.e. HO appear)
      d.B1.HO= d.mn.HO                    # Duration of HO stage        (hours; fixed)
      d.B1.PF=if(POF_TempDep){
        # POF duration for current batch as a function of water temperature based on empirical water temp data
        D.Buoy$d.PF[match(floor((t.B1.spn-1)/24),D.Buoy$Date.j)] # hours; temperature dependent
      }else{d.mn.PF}                                           # hours; fixed



      d.B1= d.B1.HO+d.B1.PF               # Duration of HO+PF stages    (hours)

      # Batch fecundity of current batch
      # Mean batch fecundity based on body size
      if(nm_fn_fb=="pow3"){
        fb.mean <- fn_fb(x=Lbin[i], fb.a, fb.b, fb.c)
      }else{
        fb.mean <- fn_fb(x=Lbin[i], fb.a, fb.b)
      }

      fb.i <- fb.mean # NOTE: This used to be round(fb.mean), since you can't have a fraction
      # of an egg, but I removed the rounding because when you fit a model
      # to the simulated data, if it is rounded you get very slightly different
      # model parameters. This is probably not necessary, but I want to assure
      # that the perfect model is really perfect.
      # Add error if desired
      if(fb_size_error){
        fb.i= rnbinom(n=1,mu=fb.mean,size=fb.k)
      }



      # Duration of batch interval
      d.B1.BI= local({
        if(SI_size_error){
          ## Use bootstrap parameter values
          # Conduct Monte-Carlo draws from fitted parameters, incorporating covariance structure
          # Draw parameter values from multivariate normal distribution
          pars <- local({
            pars_mean <- c("m"=Ps.m,"s"=Ps.s,"a"=Ps.a,"f"=Ps.f)
            pars <- tmvtnorm::rtmvnorm(n=1,
                                       mean=  pars_mean,
                                       sigma= covmat_lgs_msaf,
                                       lower=c("m"=0,"s"=1,"a"=0,f=0),
                                       upper=c("m"=Inf,"s"=Inf,"a"=1,f=1))
            dimnames(pars)[[2]] <- names(pars_mean)
            as.data.frame(round(pars,4))
          })

        }else{pars=data.frame("m"=Ps.m,"s"=Ps.s,"a"=Ps.a,"f"=Ps.f)}
        P=lgs_msaf(m=pars$m, s=pars$s, a=pars$a, f=pars$f, x=Lbin[i]) # Calculate predicted probabilities at length
        if(BI_size_independent){  # Override P if batch interval is set to be size-independent
          P=Ps.mn}
        P=pmin(SI_cf*P,1) # Note that SI_cf is used to scale P based on empirical data, therefore
        # the durations that HO and POF are set to in the simulation do not
        # actually affect the setting of SI_cf
        I=round((1/P)*24) # Convert to hours and round to nearest hour
        return(I)
      })

      # Record structure age and type for each batch in the appropriate data matrices
      if(t.B1.HO<h_sim){
        # Batch1 age
        b=local({a=c(t.B1.HO:(t.B1.HO+d.B1)); a[a<h_sim]})  # Vector indices restricted to simulation time
        B1.a[i,b]= b-min(b); rm(b)

        # From initiation up to spawning...
        b=local({a=c(t.B1.HO:(t.B1.spn-1)); a[a<h_sim]})    # Vector indices restricted to simulation time
        # Record HO presence
        B1.o[i,b]=1                                  # Set structure type to 1 (HO)
        # Record HO age
        B1.HO[i,b]=b-min(b)+1
        # Record batch fecundity
        B1.fb[i,b]=fb.i                                  # Record batch fecundity from initiation up to spawning
        rm(b)

      }

      if(t.B1.spn<h_sim){
        # Record PF presence
        b=local({a=c(t.B1.spn:(t.B1.spn+d.B1.PF-1)); a[a<h_sim]}) # Vector indices restricted to simulation time
        B1.o[i,b]=2 ; rm(b)                                    # Set structure type to 2 (PF) from spawning up to PFs disappear

        # Record PF age
        b=local({a=c(t.B1.spn:(t.B1.spn+d.B1.PF-1)); a[a<h_sim]}) # Vector indices restricted to simulation time
        B1.PF[i,b]=b-min(b)+1 ; rm(b)                          # Record HO age from spawning up to PFs disappear

        # Record when a female spawns
        Bs[i,t.B1.spn]=1
      }

      # Add data to D.pop.HO as spawns occur
      D.pop.HO.row=nrow(D.pop.HO)+1                      # Determine the next row of D.pop.HO
      D.pop.HO[D.pop.HO.row,"POSIX.t"]=POSIX.t[t.B1.HO]  # Add POSIX.t
      D.pop.HO[D.pop.HO.row,"Lbin"]=Lbin[i]                # Add Lbin
      D.pop.HO[D.pop.HO.row,"HO"]=1                      # Add HO

      # Add data to D.pop.Bs as spawns occur
      D.pop.Bs.row=nrow(D.pop.Bs)+1                      # Determine the next row of D.pop.Bs
      D.pop.Bs[D.pop.Bs.row,"POSIX.t"]=POSIX.t[t.B1.spn] # Add POSIX.t
      D.pop.Bs[D.pop.Bs.row,"Lbin"]=Lbin[i]                # Add Lbin
      D.pop.Bs[D.pop.Bs.row,"Bs"]=1                      # Add Bs
      D.pop.Bs[D.pop.Bs.row,"BI"]=d.B1.BI                # Add BI
      D.pop.Bs[D.pop.Bs.row,"fb"]=fb.i                   # Add fb
      D.pop.Bs[D.pop.Bs.row,"d.B1.PF"]=d.B1.PF           # Add d.B1.PF

      # Recalculate t.B1 (a little clunky)
      # Determine the date that this recent spawn occurred on to restrict the model
      # from allowing females to spawn more than once on the same date
      t.B1=local({
        # Determine the calendar date of the current spawn
        a1=POSIX.t[t.B1.spn]
        a2=format(a1,"%Y-%m-%d")
        a3=as.POSIXct(paste(a2,"00:00",sep=" "),tz="UTC")
        # Calculate the calendar date of the next spawn
        b1=min(floorf(t.B1.spn+d.B1.BI,24)+1,h_sim) # Restrict date from going beyond the end of the simulation to avoid an error
        b2=POSIX.t[b1]    # Convert t.B1 to POSIX

        if(b2==a3){b2=b2+60*60*24} # If date of the next spawn is the same
        # date as the current spawn advance the date by 1 day
        b3=min(b2,POSIX.t[h_sim])     # Restrict date from going beyond the end of the simulation to avoid an error

        which(POSIX.t==b3) # Convert t.B1 back to h_sim format
      })
    } #@@@@@ END TIME STEP LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    D.pop.Bs$POSIX.t=as.POSIXct(D.pop.Bs$POSIX.t,origin="1970-01-01",tz="UTC")


    # First and last simulation hour of the spawning period
    t.Pd0[i]=min(which(Bs[i,]==1))
    t.Pd1[i]=max(which(Bs[i,]==1))

    # First and last date of the spawning period
    date.Pd0[i]=as.Date(format(POSIX.t[t.Pd0[i]],"%Y-%m-%d"))
    date.Pd1[i]=as.Date(format(POSIX.t[t.Pd1[i]],"%Y-%m-%d"))

    # First and last simulation hour of the spawning capable period
    # Calculate first hour an HO is present, then determine the first hour (0:00h) of that day
    # (i.e. assume that the female has been capable of spawning the entire day that you first observe
    #  an HO. This is a simplifying assumption that should have minimal effect on calculations)
    t.SC0[i]=floorf(min(which(B1.HO[i,]>0)),24)+1
    # Calculate last hour a POF is present, then determine the last hour (23:00h) of that day
    # (i.e. assume that the female is capable of spawning the entire day that you last observe
    #  a POF. This is a simplifying assumption that should have minimal effect on calculations)
    t.SC1[i]=ceilingf(max(which(B1.PF[i,]>0)),24)

    # Record spawning capability data
    B1.SC[i,c(t.SC0[i]:t.SC1[i])]=1

    # Vector of batches by Julian day
    Bs.d=as.numeric(tapply(X=Bs[i,],INDEX=ord.day,FUN=sum)) # Batches spawned per day
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
    Ds[i]=as.numeric(date.Pd1[i]-date.Pd0[i])+1     # Spawning period duration
    nb[i]=sum(Bs[i,])                               # Number of batches
    fb[i]=exp(mean(log(D.pop.Bs[D.pop.Bs$Lbin==Lbin[i],"fb"]))) # (mean) batch fecundity logged and then back-transformed, since the error is lognormal
    fa[i]=sum(D.pop.Bs[D.pop.Bs$Lbin==Lbin[i],"fb"])  # Annual fecundity
  } #@@@@@ END FEMALE SIZE CLASS LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # Estimated batch fecundity for population
  mp.pop <- local({
    x <- D.pop.Bs$Lbin; y <- D.pop.Bs$fb
    fb.thresh <- min(x)

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
        fb.min <- round(pow3(x=fb.thresh,a=mp.pop[1],b=mp.pop[2],c=mp.pop[3]))
      }else{
        mp.pop <- as.numeric(fit_pow3_norm(x,y)$par)
        fb.min <- round(pow3(x=fb.thresh,a=mp.pop[1],b=mp.pop[2],c=mp.pop[3]))
      }
      mp.pop[4] <- fb.thresh
      mp.pop[5] <- fb.min
    }
    return(mp.pop)
  })

  # my.set.seed(iter=iter.i) # Reset random seed (2025-08-13 I think that this is being reset more than necessary within each iteration, but it's not a problem.)

  if(nm_fn_fb=="pow3"){
    fb <- c(rep(0,length(Lbin_imm)),fn_fb(Lbin_mat, mp.pop[1], mp.pop[2], mp.pop[3]))}else{
      fb <- c(rep(0,length(Lbin_imm)),fn_fb(Lbin_mat, mp.pop[1], mp.pop[2]))
    }

  # Sort D.pop.Bs by POSIX.t (i.e. put the spawn data in chronological order)
  D.pop.Bs <- D.pop.Bs[order(D.pop.Bs$POSIX.t),]
  # Add columns to D.pop.Bs
  D.pop.Bs$Date=as.Date(format(D.pop.Bs$POSIX.t,"%Y-%m-%d"))

  #-- Compile batch number data for the entire population during the simulation
  # First day of spawning periods
  D.pop.par[,"date.Pd0"]=as.Date(date.Pd0,origin=as.Date("1970-01-01"))
  # Last day of spawning periods
  D.pop.par[,"date.Pd1"]=as.Date(date.Pd1,origin=as.Date("1970-01-01"))
  # Spawning period duration
  D.pop.par[,"Ds"]=Ds
  # Spawning fraction
  D.pop.par[,"Ps"]=Ps
  # Batch interval
  D.pop.par[which(Ps>0),"BI"]=1/Ps[which(Ps>0)]
  # Number of batches
  D.pop.par[,"nb"]=nb
  # (mean) Batch fecundity
  D.pop.par[,"fb"]=fb
  # Annual fecundity
  D.pop.par[,"fa"]=fa

  # Build VERY LARGE data frame where there is a row for each female for each step in the simulation
  # Initialize data frame
  D.pop=data.frame("t"=rep(t_sim,each=n_Lbin),                                          # Simulation time
                   "Date"=rep(as.Date(format(POSIX.t[t_sim],"%Y-%m-%d")),each=n_Lbin),  # Date
                   "Time"=rep(format(POSIX.t[t_sim],"%H:%M"),each=n_Lbin),              # Time (24 hour)
                   "size.class"=rep(size.class,h_sim),                            # Size class
                   "Lbin"=rep(Lbin,h_sim),                                          # Female length
                   "f"=rep(1,h_sim),                                                # Female (0=no, 1=yes; used for counting females when summarizing data)
                   "fm"=rep(fm,h_sim))                                            # Female maturity (0=immature, 1=mature)
  # Add batches spawned
  D.pop <- add_col(df1=D.pop,
                    df2=reshape_lite(df=as.data.frame(Bs), pattern="t.",names=c("size.class","t","Bs")),
                    nm_col1="size.class",nm_col2="t",nm_x="Bs")
  # Add SC
  D.pop <- add_col(df1=D.pop,
                    df2=reshape_lite(df=as.data.frame(B1.SC),pattern="t.",names=c("size.class","t","SC")),
                    nm_col1="size.class",nm_col2="t",nm_x="SC")
  # Add HO
  D.pop <- add_col(df1=D.pop,
                    df2=reshape_lite(df=as.data.frame((B1.HO>0)*1),pattern="t.",names=c("size.class","t","HO")),
                    nm_col1="size.class",nm_col2="t",nm_x="HO")
  # Add PF
  D.pop <- add_col(df1=D.pop,
                    df2=reshape_lite(df=as.data.frame((B1.PF>0)*1),pattern="t.",names=c("size.class","t","PF")),
                    nm_col1="size.class",nm_col2="t",nm_x="PF")

  # Add column indicating if any spawning indicators were present (i.e. HO and/or PF)
  D.pop$SI <- pmin(rowSums(D.pop[,c("HO","PF")]),1)

  # Sort data
  o=order(D.pop$t,D.pop$size.class)
  D.pop=D.pop[o,]


  # #-- True prevalence of structures in population at each time step (used in plots)
  # #-- Sum presence of each structure for each female for each time step
  # for(i in 1:h_sim) {P.SC[i]=length(B1.SC[,i][B1.SC[,i]>0])/length(B1.SC[,i])}  # Proportion SC
  # for(i in 1:h_sim) {P.HO[i]=length(B1.HO[,i][B1.HO[,i]>0])/length(B1.HO[,i])}  # Proportion with HO
  # for(i in 1:h_sim) {P.PF[i]=length(B1.PF[,i][B1.PF[,i]>0])/length(B1.PF[,i])}  # Proportion with PF
  #
  # # True count count data for population at each hour of the simulation
  # B1.HO.PA=(B1.HO>0)*(24/d.mn.HO)  # Convert B1.HO to presence-absence scaled by HO duration
  # B1.PF.PA=(B1.PF>0)*(24/d.mn.PF)  # Convert B1.PF to presence-absence scaled by PF-duration
  #
  # # Convert data to long-form binary data sets for model fitting
  #
  # # B1.SC
  # B1.SC.LB=LB(B1.SC,"SC",size.class=which(Lbin%in%Lbin_mat))
  # # B1.HO
  # B1.HO.LB=LB((B1.HO>0)*1,"HO",size.class=which(Lbin%in%Lbin_mat))
  # # B1.PF
  # B1.PF.LB=LB((B1.PF>0)*1,"PF",size.class=which(Lbin%in%Lbin_mat))

  return(
    list(D.pop.par=D.pop.par, D.pop=D.pop, Bs=Bs, B1.HO=B1.HO)
    )
}
