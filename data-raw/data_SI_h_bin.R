## code to prepare `data_SI_h_bin` dataset goes here

# SI~h (spawning indicator~hour of the day) where each row represents a single fish
data_SI_h_bin <- read.csv(file="data_SI_h_bin.csv")
usethis::use_data(data_SI_h_bin, overwrite = TRUE)
