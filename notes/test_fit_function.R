library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
require(tidyverse)
source("./notes/generate_simu_dat.R")


dat <-
generate_simulated_data(300, 200, 3, 5, 3, 0.7,
                        prepare_for_fitting = T,mv_coeffs = T,seed = 2025)




Y = dat$fit_data$train
X = dat$fit_data$X$Q
Z = NULL
r = 3
lambda_M = .001
lambda_beta = .0001
lambda_gamma = .001
intercept_row = T
intercept_col = F
L_a = NULL
lambda_a = 0
L_b = NULL
lambda_b = 0
maxit = 15
thresh = 1e-5
trace = T
warm_start = NULL
ls_initial = T

fit <- IMR::imr.fit(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,

                r=6, lambda_M = 3.23, lambda_beta=.000, lambda_gamma=0,
                trace=T, ls_initial = T,intercept_row = T, intercept_col = T)
quick_camc_simu_res(dat, fit)

fitm <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,
                                 .0001, .0001, T, T, trace=T)

fit2 <- IMR::imr.cv.M(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$X$Q,
                      dat$fit_data$Z$Q, dat$fit_data$Y_full, 0, 0, T, T)

quick_camc_simu_res(dat, fit2$fit)


microbenchmark::microbenchmark(
a = partial_crossprod(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
b =  partial_crossprod_cpp(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
times = 500
)
