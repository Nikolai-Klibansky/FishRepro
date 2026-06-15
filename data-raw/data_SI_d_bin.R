## code to prepare `data_SI_d_bin` dataset goes here

# SI~d (spawning indicator~day of year) where each row represents a single fish
data_SI_d_bin <- read.csv(file="data_SI_d_bin.csv")

usethis::use_data(data_SI_d_bin, overwrite = TRUE)
