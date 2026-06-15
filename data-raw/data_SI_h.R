## code to prepare `data_SI_h` dataset goes here

# SI~h (spawning indicator~hour of the day)
data_SI_h <- read.csv(file="data_SI_h.csv")

usethis::use_data(data_SI_h, overwrite = TRUE)
