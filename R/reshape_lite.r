#' Quickly reshape a large data frame from wide to long format
#'
#' Quickly reshape a large data frame from wide to long format
#' @param df wide form data frame
#' @param pattern pattern to delete in column names of `df`
#' @param names names of column names in output data set
#' @author Nikolai Klibansky
#' @note This function doesn't do as much as \code{\link[stats]{reshape}} but it's way faster for reshaping certain large data frames from wide to long format.
#' @export
#' @examples
#' # Simulate a data set representing hydrated oocyte presence (by femele size class, by hour)
#' d_sim <- 365                      # number of days in simulation
#' h_sim <- 24*d_sim                # number of hours in simulation
#' Lbin <- seq(100,1000, by=50)     # Female length bins
#' n_Lbin <- length(Lbin)           # Number of female length bins
#' L50 <- Lbin[ceiling(length(Lbin)*0.2)] # Female length: at maturity
#' Lbin_imm <- Lbin[Lbin<L50]       # Immature female length bins
#' Lbin_mat <- Lbin[Lbin>=L50]      # Mature female length bins
#'
#' fm <- c(rep(0,length(Lbin_imm)),rep(1,length(Lbin_mat))) # Female maturity status by length bin
#' size.class <- paste("f",sprintf("%02d", 1:n_Lbin),sep=".")  # Names of length bins
#' dur_HO <- 10 # duration of hydrated oocytes in hours
#'
#' # Initialize empty matrix of simulated HO presence (by size class, by hour)
#' data_HO <- matrix(data=0,nrow=n_Lbin,ncol=h_sim,
#'                   dimnames=list(females=size.class,time=paste("t",1:h_sim,sep=".")))
#'
#' # Fill in matrix with simulated values
#' HO_h1 <- 24*sort(sample(x=ceiling(d_sim*0.2):floor(d_sim*.8),size=round(d_sim*0.2)))
#' HO_h2 <- HO_h1+dur_HO
#' HO_h <- lapply(1:length(HO_h1),FUN=function(i){HO_h1[i]:HO_h2[i]})
#' for(i in seq_along(HO_h)){
#'   xi <- HO_h[[i]]
#'   data_HO[,xi] <- fm
#' }
#'
#' data_HO_long <- reshape_lite(df=as.data.frame(data_HO), pattern="t.",names=c("size.class","t","HO"))
#'

reshape_lite <- function(df,pattern,names){
  x <- rep(rownames(df),)
  x2 <- rep(as.numeric(sub(pattern=pattern,replacement="",x=names(df))),each=nrow(df))
  x3 <- as.numeric(unlist(df))
  df_out <- data.frame(x,x2,x3)
  if(!missing(names)){names(df_out)=names}
  return(df_out)
}
