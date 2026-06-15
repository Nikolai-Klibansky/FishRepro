## code to prepare `data_SI_L` dataset goes here

# SI~L (spawning indicator~length)
data_SI_L <- read.csv(file=paste("data_SI_L.csv"))

usethis::use_data(data_SI_L, overwrite = TRUE)
