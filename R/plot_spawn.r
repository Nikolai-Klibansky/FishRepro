#' Plot spawning activity over time
#'
#' Plot proportion of mature females spawning and with spawning indicators by time unit, using data produced by \code{sim_spawn} and \code{sim_sample}. Since spawning events occur instantaneously and evidence of these events (i.e. presence of HO or POF) is brief, plotting the data in a useful way can be challenging.
#' !!! FUNCTION UNDER CONSTRUCTION!!!
#' @param data_pop Population data returned by \code{sim_spawn}
#' @param data_smp Sample data returned by \code{sim_sample}
#' @param time_unit unit of time to use for grouping observations for plotting proportions defaults to "day", but could also be "week" or "month"
#' @param nm_spnr Name of column which represents numbers of spawners.
#' @param year Reference year for plotting data over time. Used to convert dates in data.
#' @param dur_Bs Duration of spawning events in hours. Spawning events are effectively instantaneous and are evident for 1 hr.
#' @param dur_SC Duration of spawning capable stage in hours. Even though spawning capable stage may last for months, it's best to set this to 24 hours, resulting in a correction factor of 1. (cf_SC = 24/dur_SC)
#' @param dur_HO Duration of hydrated oocytes (HO) in hours.
#' @param dur_PF Duration of post-ovulatory follicles (POF) in hours.
#' @param dur_SI Duration of spawning indicators in hours. Defaults to the sum of dur_HO and dur_POF
#' @param correct_pop Should population parameters, corrected based on dur, be plotted? If TRUE, plots of true proportions will be based on corrected data. If FALSE, the raw values will be plotted. Note that corrected values will tend to align better than uncorrected.
#' @param args_matplot Additional arguments to pass to \code{\link[base]{matplot}}.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#'

plot_spawn <- function(data_pop,
                       data_smp=NULL,
                       time_unit="month",
                       nm_spnr="fm",
                       year="1970",
                       dur_Bs=1,
                       dur_SC=24,
                       dur_HO=12,
                       dur_PF=24,
                       dur_SI=dur_HO+dur_PF,
                       correct_pop=TRUE,
                       args_matplot=list()){
  # Variable names
  nm_dem <- c("f","fm") # Demographic
  nm_spn  <- c("Bs","SC","HO","PF","SI") # Spawning
  nm_spn_c <- c("Bs_c","SC_c","HO_c","PF_c","SI_c") # Spawning corrected

  # Modify data_pop
  data_pop <- as.data.frame(apply(X=data_pop[,c(nm_dem,nm_spn)],MARGIN=2,FUN=function(x) tapply(x,INDEX=data_pop[,"Date"],FUN=sum)))
  data_pop <- cbind(data.frame("Date"=as.Date(row.names(data_pop))),data_pop)

  # Modify data_smp
  if(is.null(data_smp)){
    message("data_smp is NULL. Setting data_smp equal to data_pop")
    data_smp <- data_pop}

  # Subset data_smp to only include "spawners"
  data_smp <- data_smp[data_smp[,nm_spnr]==1,]

  # Summarize data by time unit
  tud1 <- paste(time_unit,"day1",sep=".")
  jan1 <- as.Date(paste(year,"01","01",sep="-")) # First date of the year

  # Summarize data by day and compute corrected counts
  D_day <- local({
    a <- data.frame(apply(X=data_pop[,c(nm_dem,nm_spn)],MARGIN=2,FUN=function(x) tapply(x,INDEX=data_pop[,"Date"],FUN=sum)))
    a <- cbind(data.frame("day"=format(as.Date(row.names(a)),"%j"),
                           "Date"=as.Date(row.names(a))),a)
    a$Bs_c <- (24/dur_Bs)*a$Bs
    a$SC_c <- (24/dur_SC)*a$SC
    a$HO_c <- (24/dur_HO)*a$HO
    a$PF_c <- (24/dur_PF)*a$PF
    a$SI_c <- (24/dur_SI)*a$SI
    a
  })



  if(time_unit=="day"){
    # Spawn data
    data_pop_tu <- D_day
  }

  if(time_unit=="week"){
    # Spawn data
    D_day$week <- format(D_day$Date,"%U")
    data_pop_tu <- local({
      a <- data.frame(apply(X=D_day[,c(nm_dem,nm_spn,nm_spn_c)],MARGIN=2,FUN=function(x) tapply(x,INDEX=D_day[,"week"],FUN=sum)))
      cbind(data.frame("week"=row.names(a),
                       "week.day1"=as.Date(jan1-as.numeric(format(jan1,"%w"))+as.numeric(row.names(a))*7)),a)
    })
  }

  if(time_unit=="month"){
    # Spawn data
    D_day$month <- format(D_day$Date,"%m")
    data_pop_tu <- local({
      a <- data.frame(apply(X=D_day[,c(nm_dem,nm_spn,nm_spn_c)],MARGIN=2,FUN=function(x) tapply(x,INDEX=D_day[,"month"],FUN=sum)))
      cbind("month"=row.names(a),
            "month.day1"=as.Date(paste(year,row.names(a),"01",sep="-")),a)
    })
    }

  names(data_pop_tu)[1:2] <- c(time_unit,tud1)

  ## Calculate proportions with traits by time unit
  # Observed
  data_pop_tu$P_Bs <- data_pop_tu$Bs/data_pop_tu[,nm_spnr]
  data_pop_tu$P_SC <- data_pop_tu$SC/data_pop_tu[,nm_spnr]
  data_pop_tu$P_HO <- data_pop_tu$HO/data_pop_tu[,nm_spnr]
  data_pop_tu$P_PF <- data_pop_tu$PF/data_pop_tu[,nm_spnr]
  data_pop_tu$P_SI <- data_pop_tu$SI/data_pop_tu[,nm_spnr]
  # Corrected
  data_pop_tu$P_Bs_c <- data_pop_tu$Bs_c/data_pop_tu[,nm_spnr]
  data_pop_tu$P_SC_c <- data_pop_tu$SC_c/data_pop_tu[,nm_spnr]
  data_pop_tu$P_HO_c <- data_pop_tu$HO_c/data_pop_tu[,nm_spnr]
  data_pop_tu$P_PF_c <- data_pop_tu$PF_c/data_pop_tu[,nm_spnr]
  data_pop_tu$P_SI_c <- data_pop_tu$SI_c/data_pop_tu[,nm_spnr]

  # Estimated proportion with spawning indicators
  # # SC
  # D_PSC <- calc_spawn(Date=data_smp$Date, y=data_smp$SC, dur_SI=dur_SC, time_unit=time_unit)
  #
  # # HO
  # D_PHO <- calc_spawn(Date=data_smp$Date, y=data_smp$HO, dur_SI=dur_HO, time_unit=time_unit)
  #
  # # PF
  # D_PPF <- calc_spawn(Date=data_smp$Date, y=data_smp$PF,dur_SI=dur_PF, time_unit=time_unit)
  #
  # # SI
  # D_PSI <- calc_spawn(Date=data_smp$Date, y=data_smp$SI,dur_SI=dur_SI, time_unit=time_unit)

  # PLOT
  Tr <- 0.75 # Transparency
  cols <- list(
    "spawn"=rgb(0,0,0,Tr),          # black
    "SC"=   rgb(1.00,0.55,0.00,Tr), # darkorange
    "HO"=   rgb(0.33,0.10,0.55,Tr), # purple4
    "POF"=  rgb(1.00,0.08,0.58,Tr), # deeppink
    "SI"=   rgb(0,0.7,0,Tr))        # darkgreen
  types <- list(
    "spawn"="o",
    "SC"=   "l",
    "HO"=   "l",
    "POF"=  "l",
    "SI"=   "l"
  )
  cex <- .75

  ltys <- list("true"="solid",
               "est"="dashed"
  )
  # pch.SC <- 16
  # pch.HO <- 16
  # pch.PF <- 16
  # pch.SI <- 16

  pch <- 16

  # True proportions spawning
  args_matplot_default <- list(col=unlist(cols), lty=ltys$true, lwd=1, pch=pch, type=types,
                              ylim=c(0,1),xlab="Time", ylab="Proportion")
  args_matplot <- modifyList(args_matplot_default, args_matplot)

  nm_plot <- if(correct_pop){paste("P_",nm_spn_c,sep="")}else{paste("P_",nm_spn,sep="")}
  do.call(matplot,c(list(x=data_pop_tu[,tud1],y=data_pop_tu[,nm_plot]),args_matplot))

  # Estimated proportion with spawning indicators
  # points(D_PSC[,tud1],D_PSC$P_SI, lty=lty.SC, lwd=2, col=cols[3], type="l")
  # points(D_PHO[,tud1],D_PHO$P_SI, lty=lty.HO, lwd=2, col=cols[1], type="l")
  # points(D_PPF[,tud1],D_PPF$P_SI, lty=lty.PF, lwd=2, col=cols[2], type="l")
  # points(D_PSI[,tud1],D_PSI$P_SI, lty=lty.SI, lwd=2, col=cols[4], type="l")

  legend("topright",legend=c("True spawns","True SC","True HO","True POF","True SI","Est. SC","Est. HO","Est. POF", "Est. SI"),
         bty="n",col=rep(unlist(cols),2),cex=cex,lty=rep(unlist(ltys),each=5),
         pch=c(pch,rep(NA,8)))

  return(data_pop_tu)
}
