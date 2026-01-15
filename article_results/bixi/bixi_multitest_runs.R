library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
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
date = "Jan14"
test_pct = 0.5

o = prepare_bixi_data(test_pct, date)





