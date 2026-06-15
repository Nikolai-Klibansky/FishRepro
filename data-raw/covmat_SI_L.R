## code to prepare `covmat_SI_L` dataset goes here

# covariate matrix resulting from fitting 4-parameter logistic model to SI~L data
covmat_SI_L <- as.matrix(read.csv(file=paste("covmat_SI_L.csv"),row.names = 1))

usethis::use_data(covmat_SI_L, overwrite = TRUE)
