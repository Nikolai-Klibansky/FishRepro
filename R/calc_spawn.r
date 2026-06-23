#' Calculate number of batches produced by each spawning capable female
#'
#' Calculate number of batches produced by each spawning capable female
#' !!! FUNCTION UNDER CONSTRUCTION!!!
#' @param Date A vector of dates as a factor with numeric levels (e.g. "001", "365")
#' @param y A vector of spawning indicator (SI) presence/absence as a binary variable
#' @param nb.calc logical. Indicate if batch number calculations should be made
#' @param SC_pd Limits of main SC period. NEEDED ONLY for batch number calculations
#' @param dur_SI Estimated duration of spawning indicator in hours (defaults to 24 hours)
#' @param fd1 First date of the fiscal year (or calendar year is all dates occur within it). Defaults to the first date of 1970 (R standard origin date)
#' @param print.results logical. Indicate whether or not print out results. defaults to FALSE
#' @param time_unit unit of time to use for grouping observations for plotting proportions defaults to "day", but could also be "week" or "month"
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#'

    calc_spawn <- function(Date,
                   y,
                   nb.calc=FALSE,
                   SC_pd=NA,
                   dur_SI=24,
                   fd1="1970-01-01",
                   time_unit="day",
                   print.results=FALSE)              {

        if(length(unique(Date))<2){warning("Only one sampling event in dataset")}
        if(sum(y)==0){warning("No spawning indicators in dataset")}

      y.name=deparse(substitute(y))

      # Sort Date and y by Date
        od=order(Date)
        Date=Date[od]
        y=y[od]

      y=factor(y,levels=1:0)

      if(class(Date)=="Date") {Date=format(Date,"%j")}# If Date is given in R-date format, convert it to "001", "365"
      Date=as.numeric(as.character(Date))         # Convert Date into a numeric value
      oDate=order(Date); Date=Date[oDate]; y=y[oDate]         # Put Date and y in ascending order of Date

      if(class(SC_pd)=="Date") {SC_pd=format(SC_pd,"%j")}# If SC_pd is given in R-date format, convert it to "001", "365"
      SC_pd=as.numeric(as.character(SC_pd)) # Convert SC_pd into a numeric value
      fd1=as.Date(fd1)                      # Coerce fd1 into R Date format

      #   If you want to make batch number calculations..
        if(!FALSE%in%c(nb.calc==TRUE,!is.na(SC_pd))) { # if both nb.calc and SC_pd are supplied
            obs.SC_pd=Date>=SC_pd[1]&Date<=SC_pd[2]    # determine which observations are in the SC period
            Date= Date[obs.SC_pd]                      # and only use those Date..
            y= y[obs.SC_pd] }                          # ..and y values

      tab1=table(Date,y) # Make cross table where values are counts
        if(ncol(tab1)==1) {tab1=cbind(tab1,"1"=rep(0,nrow(tab1)))}
      tab2=t(tab1)                      # Transpose table 1
          row.names(tab2)=c("YES","NO") # Name rows
      tab3=prop.table(tab2,margin=2)    # Convert tab2 values to column proportions

      tab4=t(tab2)
      names(dimnames(tab4))=c("Date",y.name)

      dsmp=unique(Date)                 # Sampling dates, in day of year format

      dfy=as.Date(fd1)+dsmp-1                 # Sampling dates, in fiscal year (fy) equivalent, R date format
      success=as.vector(tab2["YES",])         # Number of successes (SI present) on each date
      trials=as.vector(colSums(tab2))         # Number of trials (females at the correct phase) on each date
      P.SI=pmin((success/trials)*24/dur_SI,1) # Spawning fraction for each date, adjusted by dur_SI and restricted to
                                              #   not exceed 1.
      D.nb=data.frame(dsmp,dfy,success,trials,P.SI) # Data frame
      names(D.nb)=c("Day","Date.fy","success","trials","P.SI") # Name columns

      # BATCH NUMBER CALCULATIONS
      if (nb.calc==TRUE) {
        ### INTERVAL METHOD ###
        #     NOTE that in the interval method, sampling dates are not restricted to the range
        #     of dates between the first and last evidence of spawning indicators. Thus, the area under
        #     the piecewise function will sometimes extend beyond this date range.
          nb= NULL     # Initialize null vector for batch measurements
            for (k in (1:(length(dsmp)-1))) {
              cf=ifelse(k==1,yes=1,no=0)  # Correction factor for the first interval so
                                          # that it includes the first date of spawning
              # Calculate area under P.SI function between two sampling dates
                A.tri=as.numeric((abs(P.SI[k+1]-P.SI[k])*(dsmp[k+1]-dsmp[k]))/2)  # Area of triangle
                A.rec=as.numeric(min(P.SI[k+1],P.SI[k])*(dsmp[k+1]-(dsmp[k]-cf))) # Area of rectangle
                nb[k]=A.tri+A.rec
                                            }

          #    If both nb.calc and SC_pd are supplied, batch number by the interval method
          # is calculated correctly, otherwise return NA.
            nb.I = ifelse(!FALSE%in%c(nb.calc==TRUE,!is.na(SC_pd)),
              yes=sum(nb),
              no=NA)

        ### STANDARD METHOD ###
          dSI=D.nb$Day[which(D.nb$P.SI!=0)] # Dates where P.SI!=0
            if(length(dSI>0)){
              dA=min(dSI)  # First date where SI were present
              dZ=max(dSI)  # Last date where SI were present

              dur=as.numeric(dZ-(dA-1))     #   Length of spawning season (with -1 so that the first
                                            # day of the spawning season is included in the spawning season.)
              dAtoZ=dsmp[dsmp>=dA&dsmp<=dZ] #   Sampling dates between the first and last dates in the spawning season
                                            # (Note this calculation includes dates within this range where no
                                            # females had spawning indicators.)

              # seasonal spawning fraction
                success.sn=sum(D.nb$success[match(dAtoZ,D.nb$Day)])
                trials.sn=sum(D.nb$trials[match(dAtoZ,D.nb$Day)])
                Ps =pmin((success.sn/trials.sn)*24/dur_SI,1)

              # 95% confidence limits for seasonal spawning fraction
                Ps.95CI.L=binom.test(x=c(success.sn,trials.sn-success.sn))$conf.int[1]
                Ps.95CI.U=binom.test(x=c(success.sn,trials.sn-success.sn))$conf.int[2]

              # seasonal spawning fraction with confidence limits in parentheses
              sig=2
              Ps.wCL=paste(signif(Ps,sig)," (",signif(Ps.95CI.L,sig),"-",
                signif(Ps.95CI.U,sig),")",sep="")

              # Number of samples that go into P.sn
                Ps.n=trials.sn

              #    If both nb.calc and SC_pd are supplied, calculate seasonal spawning
              # fraction for ALL dates in the main spawning capable period
                Ps.SC = ifelse(!FALSE%in%c(nb.calc==TRUE,!is.na(SC_pd)),
                  yes=sum(D.nb$success)/sum(D.nb$trials),
                  no=NA)

                nb.S=dur*Ps # Season length * seasonal spawning fraction

        # Date that SI first appeared
          Date.fy.dA=D.nb$Date.fy[which(D.nb$Day==dA)]
        # Date that SI last appeared
          Date.fy.dZ=D.nb$Date.fy[which(D.nb$Day==dZ)]
                             }

            if(length(dSI)==0){dSI=Date.fy.dA=Date.fy.dZ=dur=Ps=Ps.95CI.L=Ps.95CI.U=Ps.n=Ps.SC=nb.S=Ps.wCL=0}


        ### PRINT STUFF ###
        if (print.results==TRUE){
          cat("\n y=",y.name,"\n")
          cat("\nInterval method: Number of batches=",round(nb.I,2),"\n")
          cat("\nStandard method: Seasonal spawning fraction=",round(Ps,2))
          cat("\n                 Seasonal spawning fraction lower 95%CL=",round(Ps.95CI.L,2))
          cat("\n                 Seasonal spawning fraction upper 95%CL=",round(Ps.95CI.U,2))
          cat("\n                                                      n=",Ps.n)
          cat("\n                 Spawning season duration=",round(dur,2))
          cat("\n                 Number of batches=",round(nb.S,2),"\n")
          cat("\nSpawning fraction over entire SC period=",round(Ps.SC,2),"\n\n")
                            }
      }

      # Summarize output data by time_unit, is a value is given other than "day"
        if(time_unit=="day"){
          D.nb$day=format(D.nb$Date.fy,"%j")
          D.nb$day.day1=as.Date(D.nb$Date.fy)
        }

        if(time_unit=="week"){
          D.nb$week=format(D.nb$Date.fy,"%U")
          D.nb$week.day1=D.nb$Date.fy-as.numeric(format(D.nb$Date.fy,"%w"))
        }

        if(time_unit=="month"){
          D.nb$month=format(D.nb$Date.fy,"%m")
          D.nb$month.day1=as.Date(format(D.nb$Date.fy,"%Y-%m-01"))
        }
        D.nb.out<<-D.nb
        D.nb=cbind(unique(D.nb[,time_unit]),unique(D.nb[,paste(time_unit,"day1",sep=".")]),
          data.frame(apply(X=D.nb[,c("success","trials")],MARGIN=2,FUN=function(x)
          tapply(x,INDEX=D.nb[,time_unit],FUN=sum,simplify=FALSE))),
          stringsAsFactors=FALSE)
        names(D.nb)=c(time_unit,paste(time_unit,"day1",sep="."),c("success","trials"))
        D.nb$success=as.numeric(D.nb$success) # Simplify success column
        D.nb$trials=as.numeric(D.nb$trials)   # Simplify trials column
        D.nb$P.SI=pmin((D.nb$success/D.nb$trials)*24/dur_SI,1)

        if(nb.calc){
          out=list(y.name=y.name,Ps=Ps,
            Ps.95CI.L=Ps.95CI.L, Ps.95CI.U=Ps.95CI.U,Ps.wCL=Ps.wCL,
            Ps.n=Ps.n,
            Date.fy.dA=Date.fy.dA, Date.fy.dZ=Date.fy.dZ, dur=dur,
            nb.S=nb.S,nb.I=nb.I,Ps.SC=Ps.SC,
            data=D.nb)
        }else{out=D.nb}

        out
    }
