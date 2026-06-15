## code to prepare `data_SI_d` dataset goes here

# SI~d (spawning indicator~day of the year)
data_SI_d <- read.csv(file="data_SI_d.csv")

usethis::use_data(data_SI_d, overwrite = TRUE)
