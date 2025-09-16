library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
devtools::uninstall(); devtools::install()
require(tidyverse)
source("./notes/generate_simu_dat.R")


dat <-
generate_simulated_data(300, 200, 3, 5, 3, 0.7,
                        prepare_for_fitting = T,mv_coeffs = T,seed = 2025)




Y = dat$fit_data$train

Y = dat$fit_data$Y_full
X = dat$fit_data$X$Q
val_prop = 0.2
hpar = IMR::get_imr_default_hparams()
Z = NULL
r = 3
lambda_M = .001
lambda_beta = NULL#.0001
lambda_gamma = NULL#.001
intercept_row = T
intercept_col = F
L_a = NULL
lambda_a = 0
L_b = NULL
lambda_b = 0
maxit = 15
thresh = 1e-5
trace = T
seed=2023
verbose=3
error_function <- IMR::error_metric$rmse
warm_start = NULL
ls_initial = T

fit <- IMR::imr.fit(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,

                r=6, lambda_M = 3.23, lambda_beta=.000, lambda_gamma=0,
                trace=T, ls_initial = T,intercept_row = T, intercept_col = T)
quick_camc_simu_res(dat, fit)

fitm <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,
                                 .0001, .0001, T, T, trace=T)

fit2 <- IMR::imr.cv_M(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$X$Q,
                      dat$fit_data$Z$Q, dat$fit_data$Y_full, 0, 0, T, T)

quick_camc_simu_res(dat, fit2$fit)


IMR::get_lambda_lasso_max(dat$fit_data$train, NULL,
                          dat$fit_data$Z$Q, dat$fit_data$valid,
                          dat$fit_data$Y_full, T, T, verbose = T
                          )



future::plan(future::multisession, workers = 4)
fit3 <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,verbose=3)

f <- IMR::imr.cv_M

worker <- function(args) {
  on.exit(if (!is.null(p)) p(), add = TRUE)
  do.call(f, c(args, list(y_train = y_train,
                          y_valid = y_valid,
                          X = X,
                          Z = Z,
                          Y_full = Y,
                          intercept_row = intercept_row,
                          intercept_col = intercept_col,
                          hpar = hpar,
                          error_function = error_function,
                          thresh = thresh,
                          maxit = maxit,
                          trace = inner_trace,
                          seed = seed)))
}

hpar$M$lambda_max = 30

IMR::imr.cv_M(
  lambda_beta = 0,
  lambda_gamma = 0,
  y_train = y_train,
  y_valid = y_valid,
  X = X,
  Z = Z,
  Y_full = Y,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  hpar = hpar,
  error_function = error_function,
  thresh = thresh,
  maxit = maxit,
  trace = T,
  seed = seed
)


IMR::imr.fit(
  lambda_beta = 0,
  lambda_gamma = 0,
  Y = y_train,
  r = 2,
  X = X,
  Z = Z,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  thresh = thresh,
  maxit = 300,
  trace = T,
  ls_initial = F
) -> o

dims <- dim(y_train)
nr <- dims[1]
nc <- dims[2]
init <- opt_svd(IMR:::naive_MC(y_train), 2, nr, nc, FALSE, FALSE)

mat <- as.matrix(y_train)


microbenchmark::microbenchmark(
a = partial_crossprod(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
b =  partial_crossprod_cpp(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
times = 500
)
