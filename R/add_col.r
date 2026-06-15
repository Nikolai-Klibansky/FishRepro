#' Add a column `x` to a data frame df1 from data frame df2, using two index columns to match rows
#'
#' Add a column `x` to a data frame df1 from data frame df2, using two index columns to match rows
#' @param df1 data frame 1. The main data frame you want to add a new column to.
#' @param df2 data frame 2. The data frame you want to add the
#' @param nm_col1 name of indexing column 1
#' @param nm_col2 name of indexing column 2
#' @param nm_x name of variable of interest `x`
#' @author Nikolai Klibansky
#' @note You can use \code{\link[base]{merge}} to do this instead, but this is faster. I'm not sure how widely useful this function is, but it is used by other functions in this package.
#' @export
#' @examples
#' # Here's a slightly obnoxious but hopefully informative example which first generates a couple of data sets
#' # Scroll down to the end of the example to see the usage of add_col.
#'
#' # Simulate a data set representing hydrated oocyte presence (by femele size class, by hour)
#' d_sim <- 365                     # number of days in simulation
#' h_sim <- 24*d_sim                # number of hours in simulation
#' t_sim <- 1:h_sim                 # Time points in simulation (in hours)
#' POSIX.t <- as.POSIXct(x=(0:(h_sim-1))*3600,origin="1970-01-01",tz="UTC") # Time points in simulation (in POSIX date and hour)
#'
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
#' # data frame for collecting population data
#' df_pop <- data.frame("t"=rep(t_sim,each=n_Lbin),                                      # Simulation time
#'                      "Date"=rep(as.Date(format(POSIX.t[t_sim],"%Y-%m-%d")),each=n_Lbin),  # Date
#'                      "Time"=rep(format(POSIX.t[t_sim],"%H:%M"),each=n_Lbin),              # Time (24 hour)
#'                      "size.class"=rep(size.class,h_sim),                                  # Size class
#'                      "Lbin"=rep(Lbin,h_sim),                                              # Female length
#'                      "f"=rep(1,h_sim),                                                    # Female (0=no, 1=yes; used for counting females when summarizing data)
#'                      "fm"=rep(fm,h_sim))
#'
#' # Add HO to df_pop
#' df_pop <- add_col(df1=df_pop,
#'                   df2=reshape_lite(df=as.data.frame(data_HO),pattern="t.",names=c("size.class","t","HO")),
#'                   nm_col1="size.class",nm_col2="t",nm_x="HO")
#'

add_col <- function(df1,df2,nm_col1,nm_col2,nm_x){
  df1_ix <- paste(df1[,nm_col1],df1[,nm_col2])   # Create indexing vector for first data frame
  df2_ix <- paste(df2[,nm_col1],df2[,nm_col2])   # Create indexing vector for second data frame

  df1[,nm_x] <- df2[match(df1_ix,df2_ix),nm_x]  # Add values of x where df1 and df2 match
  return(df1)
}
