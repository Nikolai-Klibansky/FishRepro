#' Plot spawning activity over time
#'
#' Plot proportion of mature females spawning and with spawning indicators by time unit, using data produced by \code{sim_spawn} and \code{sim_sample}. Since spawning events occur instantaneously and evidence of these events (i.e. presence of HO or POF) is brief, plotting the data in a useful way can be challenging.
#' !!! FUNCTION UNDER CONSTRUCTION!!!
#' @param data_pop Population data returned by \code{sim_spawn}
#' @param data_smp Sample data returned by \code{sim_sample}
#' @param time_unit unit of time to use for grouping observations for plotting proportions defaults to "day", but could also be "week" or "month"
#' @param nm_spnr Name of column which represents numbers of spawners.
#' @param year Reference year for plotting data over time. Used to convert dates in data.
#' @param dh_Bs Duration of spawning events in hours. Spawning events are effectively instantaneous and are evident for 1 hr.
#' @param dh_SC Duration of spawning capable stage in hours. Even though spawning capable stage may last for months, it's best to set this to 24 hours, resulting in a correction factor of 1. (cf_SC = 24/dh_SC)
#' @param dh_HO Duration of hydrated oocytes (HO) in hours.
#' @param dh_PF Duration of post-ovulatory follicles (POF) in hours.
#' @param dh_SI Duration of spawning indicators in hours. Defaults to the sum of dh_HO and dh_POF
#' @param correct_pop Should population parameters, corrected based on dur, be plotted? If TRUE, plots of true proportions will be based on corrected data. If FALSE, the raw values will be plotted. Note that corrected values will tend to align better than uncorrected.
#' @param filter_data_h Should data_pop be filtered to only include dates when spawning was observed before plotting traits by hour?
#' @param nm_trt_plot Names of traits to plot.
#' @param args_matplot Additional arguments to pass to \code{\link[base]{matplot}}.
#' @param cols Plot colors associated with each plotted trait. Bs = spawns, SC = spawning capable, HO = hydrated oocytes, PF = spawning capable, SI = spawning indicators
#' @param types Plot types associated with each plotted trait. Bs = spawns, SC = spawning capable, HO = hydrated oocytes, PF = spawning capable, SI = spawning indicators
#' @param pchs Plotting characters associated with each plotted trait. Bs = spawns, SC = spawning capable, HO = hydrated oocytes, PF = spawning capable, SI = spawning indicators
#' @param ltys Line types used in plots to differentiate between population (pop) and sample (smp) data sets.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Simulate spawning with defaults
#' set.seed(23456)
#' out_spawn <- sim_spawn()
#'
#' # plot results with defaults
#' plot_spawn(data_pop=out_spawn$data_pop)
#'
#' # Simulate sampling with defaults
#' out_sample <- sim_sample(sim_spawn_out = out_spawn)
#'
#'
#' # plot results with defaults, adding sampling data
#' plot_spawn(data_pop=out_spawn$data_pop,
#'            data_smp = out_sample$data_smp
#' )
#' }
#'

plot_spawn <- function(data_pop,
                       data_smp=NULL,
                       time_unit="month",
                       nm_spnr="fm",
                       year="1970",
                       dh_Bs=1,
                       dh_SC=24,
                       dh_HO=12,
                       dh_PF=24,
                       dh_SI=dh_HO+dh_PF,
                       correct_pop=TRUE,
                       filter_data_h=TRUE,
                       nm_trt_plot  = c("Bs","SC","HO","PF","SI"),
                       args_matplot=list(),
                       cols =  c(
                         "Bs"=  rgb(0,0,0,0.75),          # black
                         "SC"=  rgb(1.00,0.55,0.00,0.75), # darkorange
                         "HO"=  rgb(0.33,0.10,0.55,0.75), # purple4
                         "PF"=  rgb(1.00,0.08,0.58,0.75), # deeppink
                         "SI"=  rgb(0,0.7,0,0.75)),       # darkgreen
                       types = c(
                         "Bs"= "o",
                         "SC"= "l",
                         "HO"= "l",
                         "PF"= "l",
                         "SI"= "l"
                       ),
                       pchs = c(
                         "Bs"= 16,
                         "SC"= NA,
                         "HO"= NA,
                         "PF"= NA,
                         "SI"= NA
                       ),
                       ltys = c("pop"="solid",
                                "smp"="dashed"
                       )
                       ){
  # Variable names
  cnm_dem <- c("f","fm") # Demographic
  cnm_trt_pop  <- c("Bs","SC","HO","PF","SI") # Spawning traits
  cnm_trt_pop_c <- paste(cnm_trt_pop,"_c",sep="") # Spawning traits corrected
  cnm_trt_smp <- cnm_trt_pop[!grepl("Bs",cnm_trt_pop)]
  cnm_trt_smp_c <- cnm_trt_pop_c[!grepl("Bs",cnm_trt_pop_c)]
  nm_trt_plot_c <- paste(nm_trt_plot,"_c",sep="")

  # Plotting parameters
  args_matplot_default <- list(lty=ltys[["pop"]], lwd=2, cex = .75,
                               ylim=c(0,1),xlab="", ylab="Proportion",xaxt="n")
  args_matplot <- modifyList(args_matplot_default, args_matplot)

  ntp_pop <- if(correct_pop){paste("P_",nm_trt_plot_c,sep="")}else{paste("P_",nm_trt_plot,sep="")}
  ntp_smp <- paste("P_",nm_trt_plot_c[!grepl("Bs",nm_trt_plot_c)],sep="")
  nm_tr_smp <- gsub("P_([A-Z]{2}).*","\\1",ntp_smp) # used in plotting to match formatting

  # Check for data_smp
  if(is.null(data_smp)){
    warning("data_smp is NULL and will not be analyzed")
  }else{
    # Subset data_smp to only include "spawners"
    data_smp <- data_smp[data_smp[,nm_spnr]==1,]
    }

  # Summarize data by time unit
  tud1 <- local({
    a <- paste(time_unit,"day1",sep=".")
    if(time_unit=="day"){a <- time_unit}
    a
  })
  jan1 <- as.Date(paste(year,"01","01",sep="-")) # First date of the year

  ## Summarize data and compute corrected counts
  if(time_unit=="hour"){
  # by hour
  data_pop_h <- local({
    a <- data_pop
    if(filter_data_h){
    a <- a[a$date%in%unique(a$date[a$SC==1]),]
    }
    a <- data.frame(apply(X=a[,c(cnm_dem,cnm_trt_pop)],MARGIN=2,
                          FUN=function(x) tapply(x,INDEX=a[,"time"],FUN=sum)))
    a <- cbind(data.frame("hour"=as.numeric(gsub("([0-9]+):([0-9]+)","\\1",row.names(a))),
                          "time"=row.names(a)),a)
    a$Bs_c <- (24/dh_Bs)*a$Bs
    a$SC_c <- (24/dh_SC)*a$SC
    a$HO_c <- (24/dh_HO)*a$HO
    a$PF_c <- (24/dh_PF)*a$PF
    a$SI_c <- (24/dh_SI)*a$SI
    a
  })
  ## Calculate proportions with traits
  data_pop_h[,paste("P_",cnm_trt_pop,sep="")] <- data_pop_h[,cnm_trt_pop]/data_pop_h[,nm_spnr]     # Observed
  data_pop_h[,paste("P_",cnm_trt_pop_c,sep="")] <- data_pop_h[,cnm_trt_pop_c]/data_pop_h[,nm_spnr] # Corrected

  if(!is.null(data_smp)){
  data_smp_h <- local({
    a <- data_smp
    if(filter_data_h){
      a <- a[a$date%in%unique(a$date[a$SC==1]),]
    }
    a <- data.frame(apply(X=a[,c(cnm_dem,cnm_trt_smp)],MARGIN=2,
                          FUN=function(x) tapply(x,INDEX=a[,"time"],FUN=sum)))
    a <- cbind(data.frame("hour"=as.numeric(gsub("([0-9]+):([0-9]+)","\\1",row.names(a))),
                          "time"=row.names(a)),a)
    a$SC_c <- (24/dh_SC)*a$SC
    a$HO_c <- (24/dh_HO)*a$HO
    a$PF_c <- (24/dh_PF)*a$PF
    a$SI_c <- (24/dh_SI)*a$SI
    a
  })

  ## Calculate proportions with traits
  data_smp_h[,paste("P_",cnm_trt_smp,sep="")] <- data_smp_h[,cnm_trt_smp]/data_smp_h[,nm_spnr]     # Observed
  data_smp_h[,paste("P_",cnm_trt_smp_c,sep="")] <- data_smp_h[,cnm_trt_smp_c]/data_smp_h[,nm_spnr] # Corrected
  }


  # diel plot
    nm_tr_pop <- gsub("P_([A-Za-z]{2}).*","\\1",ntp_pop)
    do.call(matplot,c(list(x=data_pop_h[,"hour"],y=data_pop_h[,ntp_pop],
                           type=types[nm_tr_pop],
                           col=cols[nm_tr_pop],
                           pch=pchs[nm_tr_pop]),
                      args_matplot))
    axis(side = 1, at = pretty(data_pop_h[,"hour"]))
    title(xlab = "hour")
    out <- list(data_pop_t=data_pop_h,data_smp_t=NULL)

    if(!is.null(data_smp)){
    matpoints(x=data_smp_h[,"hour"],y=data_smp_h[,ntp_smp],
              type=types[nm_tr_smp],
              col=cols[nm_tr_smp],
              pch=pchs[nm_tr_smp],
              lwd=2,
              lty=2)
    out$data_smp_t=data_smp_h
    }

  out <- list(data_pop_t=data_pop_h, data_smp_t=data_smp_h)
  } # end if(time_unit=="hour")


  if(time_unit!="hour"){
  # by day
  data_pop_d <- local({
    a <- data.frame(apply(X=data_pop[,c(cnm_dem,cnm_trt_pop)],MARGIN=2,
                          FUN=function(x) tapply(x,INDEX=data_pop[,"date"],FUN=sum)))
    dates <- as.Date(row.names(a))
    a <- cbind(data.frame("date"=dates,
                          "day"=format(dates,"%j"),
                          "week"=format(dates,"%U"),
                          "month"=format(dates,"%m")
                          ),
               a)

    a$Bs_c <- (24/dh_Bs)*a$Bs
    a$SC_c <- (24/dh_SC)*a$SC
    a$HO_c <- (24/dh_HO)*a$HO
    a$PF_c <- (24/dh_PF)*a$PF
    a$SI_c <- (24/dh_SI)*a$SI
    a
  })
  if(!is.null(data_smp)){
    data_smp_d <- local({
      a <- data.frame(apply(X=data_smp[,c(cnm_dem,cnm_trt_smp)],MARGIN=2,
                            FUN=function(x) tapply(x,INDEX=data_smp[,"date"],FUN=sum)))
      dates <- as.Date(row.names(a))
      a <- cbind(data.frame("date"=dates,
                            "day"=format(dates,"%j"),
                            "week"=format(dates,"%U"),
                            "month"=format(dates,"%m")
      ),
      a)
      a$SC_c <- (24/dh_SC)*a$SC
      a$HO_c <- (24/dh_HO)*a$HO
      a$PF_c <- (24/dh_PF)*a$PF
      a$SI_c <- (24/dh_SI)*a$SI
      a
    })
  }

  if(time_unit=="day"){
    data_pop_t <- data_pop_d
    if(!is.null(data_smp)){data_smp_t <- data_smp_d}
  }

  if(time_unit=="week"){
    data_pop_t <- local({
      a <- data.frame(apply(X=data_pop_d[,c(cnm_dem,cnm_trt_pop,cnm_trt_pop_c)],MARGIN=2,
                            FUN=function(x) tapply(x,INDEX=data_pop_d[,"week"],FUN=sum)))
      cbind(data.frame("week"=row.names(a),
                       "week.day1"=as.Date(jan1-as.numeric(format(jan1,"%w"))+as.numeric(row.names(a))*7)),a)
    })

    if(!is.null(data_smp)){
    data_smp_t <- local({
      a <- data.frame(apply(X=data_smp_d[,c(cnm_dem,cnm_trt_smp,cnm_trt_smp_c)],MARGIN=2,
                            FUN=function(x) tapply(x,INDEX=data_smp_d[,"week"],FUN=sum)))
      cbind(data.frame("week"=row.names(a),
                       "week.day1"=as.Date(jan1-as.numeric(format(jan1,"%w"))+as.numeric(row.names(a))*7)),a)
    })
    }
  }

  if(time_unit=="month"){
    data_pop_t <- local({
      a <- data.frame(apply(X=data_pop_d[,c(cnm_dem,cnm_trt_pop,cnm_trt_pop_c)],MARGIN=2,
                            FUN=function(x) tapply(x,INDEX=data_pop_d[,"month"],FUN=sum)))
      cbind("month"=row.names(a),
            "month.day1"=as.Date(paste(year,row.names(a),"01",sep="-")),a)
    })

    if(!is.null(data_smp)){
    data_smp_t <- local({
      a <- data.frame(apply(X=data_smp_d[,c(cnm_dem,cnm_trt_smp,cnm_trt_smp_c)],MARGIN=2,
                            FUN=function(x) tapply(x,INDEX=data_smp_d[,"month"],FUN=sum)))
      cbind("month"=row.names(a),
            "month.day1"=as.Date(paste(year,row.names(a),"01",sep="-")),a)
    })
    }
  }

  # names(data_pop_t)[1:2] <- c(time_unit,tud1)

  ## Calculate proportions with traits by time unit (day, week, or month)
  data_pop_t[,paste("P_",cnm_trt_pop,sep="")] <- data_pop_t[,cnm_trt_pop]/data_pop_t[,nm_spnr]     # Observed
  data_pop_t[,paste("P_",cnm_trt_pop_c,sep="")] <- data_pop_t[,cnm_trt_pop_c]/data_pop_t[,nm_spnr] # Corrected

  if(!is.null(data_smp)){
  data_smp_t[,paste("P_",cnm_trt_smp,sep="")] <- data_smp_t[,cnm_trt_smp]/data_smp_t[,nm_spnr]     # Observed
  data_smp_t[,paste("P_",cnm_trt_smp_c,sep="")] <- data_smp_t[,cnm_trt_smp_c]/data_smp_t[,nm_spnr] # Corrected
  }

# year plot
  nm_tr_pop <- gsub("P_([A-Za-z]{2}).*","\\1",ntp_pop)
  do.call(matplot,c(list(x=data_pop_t[,tud1],
                         y=data_pop_t[,ntp_pop],
                         type=types[nm_tr_pop],
                         col=cols[nm_tr_pop],
                         pch=pchs[nm_tr_pop]),
                         args_matplot))
  axis.Date(
    side = 1,
    at = pretty(data_pop_t[,tud1]),
    format = "%b %d"
  )
  out <- list(data_pop_t=data_pop_t, data_smp_t=NULL)
  if(!is.null(data_smp)){
  matpoints(x=data_smp_t[,tud1],y=data_smp_t[,ntp_smp],
            type=types[nm_tr_smp],
            col=cols[nm_tr_smp],
            pch=pchs[nm_tr_smp],
            lwd=2,
            lty=2)
  out$data_smp_t <- data_smp_t
  }
  } # end if(time_unit!="hour")

  # Add legend to plots
  leg_text <- c(paste(ntp_pop,"pop"))
  if(!is.null(data_smp)){
    leg_text <- c(leg_text,paste(ntp_smp,"smp"))
  }

  legend("topright",
         legend=leg_text,
         col=cols[gsub("P_([A-Za-z]+).*","\\1",leg_text)],
         cex=cex,
         lty=ltys[gsub("(.*) ([A-Za-z]+)","\\2",leg_text)],
         pch=pchs[gsub("P_([A-Za-z]+).*","\\1",leg_text)],
         lwd=2,
         bty="n")


  invisible(out)
}
