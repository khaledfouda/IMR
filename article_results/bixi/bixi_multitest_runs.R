library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(dplyr)
library(magrittr)
source("./article_results/bixi/data_bixi.R")
source("./other_models/SoftImpute_cv.R")
source("./article_results/bixi/fit_bktr.R")
source("./article_results/bixi/helpers_bixi.R")
#-------------------------------------------------------------------------------
#' we will create multiple test set with the same size but with different seeds
#' First: we will run BKTR on all of them and save them to a file.
#' Start with 50 different test sets.
#------------------------------------------------------------
# step 1: generate the data
num_testsets = 50
timestamp = "Jan"
original_missing_pct = 0.1163287
test_pct <- 0.2
total_miss = test_pct + original_missing_pct
seed = 2025
for(prefix in 1:num_testsets){
  seed = seed + 1
  preprocess_bixi_data(total_miss, timestamp, seed, prefix,F)
}
#----------------------------------------------------------------
# step 2: train BKTR on all of them

for(prefix in 14:num_testsets){
  o = fit_BKTR_to_Bixi(total_miss, timestamp, prefix, seed)
}
#-------------------------------------------------------------
# step 3: train our model + soft-impute >>



dat <- prepare_bixi_data(total_miss, timestamp, seed, val_prop = 0.2, prefix="test",
                  spatial = "original", temporal="original",
                  spatial_jitter = TRUE,
                  temporal_jitter = TRUE,
                  kappa_max = 1e3,
                  tau_max =1e-1)




o = prepare_bixi_data(test_pct, "25Sep")

as.matrix(o$modd$obs_mask) * 1 -> m1
as.matrix(o$test_mask) * 1 -> m2

# m1 = 0 when test or missing
# m2 = 0 when train or missing

sum(m1 == 0) / length(m1)
sum(m2 == 0) / length(m2)

(-sum(m2 == 1) + sum(m1 == 0)) / length(m1)




mean( m2 * m1) # m2 = 1 and m1 = 1 # both test and observed? this is wrong!!
mean(())


sum(o$test_mask@x)/length(o$test_mask) #-
sum(o$modd$obs_mask@x)/length(o$test_mask)

is.na(train_df$y) -> mask1
length(test_df$y)
