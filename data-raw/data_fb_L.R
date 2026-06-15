## code to prepare `data_fb_L` dataset goes here

# fb~L (batch fecundity~length)
data_fb_L <- read.csv(file="data_fb_L.csv")
data_fb_L$fb <- as.integer(data_fb_L$fb)

usethis::use_data(data_fb_L, overwrite = TRUE)
