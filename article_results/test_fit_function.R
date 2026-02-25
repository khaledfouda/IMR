library(devtools)
clean_dll(); Rcpp::compileAttributes(); document();
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

n = 800; m = 900;
dat <-
generate_simulated_data(n, m, 3, 5, 5, 0.7,sparsity_beta = .50, sparsity_gamma = .50,
                        prepare_for_fitting = T,mv_coeffs = T,seed = 2025)

#======
d1 <-  matrix(rbinom(n*n,n,.2), n, n);d1 <- (d1 + t(d1)) / 2
d2 <-  matrix(rbinom(m*m,m,.2), m, m);d2 <- (d2 + t(d2)) / 2
S <- generate_similarity(d1, invert = T, jitter = 1);S
S2 <- generate_similarity(d2, invert = T, jitter = 1);S2
data <- IMR::prepare_data(dat$Y, dat$X, dat$Z,val_prop = 0.2,
                          similarity_rows = S, similarity_cols = S2);data
hp <- IMR::imr_hparameters();hp
convergence <- IMR::imr_convergence(400, 1e-4, FALSE,ls_initial = T); convergence

fit <- IMR::imr_fit(data, rank = 5, lambda_m = 5,
                    lambda_beta = 0.2, lambda_gamma = 0.2,
                    intercept_row = F, intercept_col = T,
                    shared_beta = F, shared_gamma = T, warm_start = NULL,
                    convergence = convergence);fit

rec <- IMR::reconstruct(fit, data)
dat$test <- (dat$theta * (1 - dat$mask)) %>% IMR::as.Incomplete()
recp <- IMR::reconstruct_partial(fit, data, dat$test, TRUE)
IMR::evaluate(recp@x, dat$test@x, "all") %>% as_tibble()


get_lambda_m_max(data, intercept_row = T, intercept_col = T,
                 lambda_beta = 0, lambda_gamma = 0,
                 shared_beta = F, shared_gamma = F, convergence = convergence)


IMR::get_lambda_lasso_max(data, target = "gamma", rank = 5, lambda_m = 0.2,
                     intercept_row = TRUE, intercept_col = TRUE,
                     shared_effects = TRUE, convergence = convergence, verbose=2)

# parameters for lambda_lasso_max
target = "beta"
rank = 5
lambda_m = 0.1
intercept_row = FALSE
intercept_col = FALSE
shared_beta = shared_gamma = shared_effects = FALSE
bisection_iter = 15
verbose = 2

#====


bench::mark(
  a = svd::propack.svd(as.matrix(data$Y),
                   neig = 1, opts = list(kmax = svd_maxit)
  )$d[1],
  a2 = IMR::svd_opt(data$Y,1)$d[1],
  iterations = 1000
)

#
# Y = dat$fit_data$train
#
# Y = dat$fit_data$Y_full
# X = dat$fit_data$X$Q
# val_prop = 0.2
# hpar = IMR::get_imr_default_hparams()
# Z = NULL
# r = 3
# lambda_m = .001
# lambda_beta = NULL#.0001
# lambda_gamma = NULL#.001
# intercept_row = T
# intercept_col = F
# L_a = NULL
# lambda_a = 0
# L_b = NULL
# lambda_b = 0
# maxit = 15
# thresh = 1e-5
# trace = T
# seed=2023
# verbose=3
# error_function <- IMR::error_metric$rmse
# warm_start = NULL
# ls_initial = T

fit <- IMR::imr.fit(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,

                r=6, lambda_m = 3.23, lambda_beta=.000, lambda_gamma=0,
                trace=T, ls_initial = T,intercept_row = T, intercept_col = T)
quick_camc_simu_res(dat, fit)

fitm <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,
                                 .0001, .0001, T, T, trace=T)

fit2 <- IMR::imr.cv_M(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$X$Q,
                      dat$fit_data$Z$Q, dat$fit_data$Y_full, 0, 0, T, T)

quick_camc_simu_res(dat, fit2$fit)

mfit <- IMR::imr.cv_M(
  y_train     = dat$fit_data$train,
  y_valid     = dat$fit_data$valid,
  X           = NULL,
  Z           = dat$fit_data$Z$Q,
  Y_full      = dat$fit_data$Y_full,
  lambda_beta = 0,
  lambda_gamma = 0,
  intercept_row = T,
  intercept_col = T,
  trace       = T
)


fit <- IMR::imr.fit(
  Y = dat$fit_data$train,
  X = NULL,
  Z = dat$fit_data$Z$Q,
  r = 7,
  lambda_m = 62.34,
  lambda_beta = 0,
  lambda_gamma = 0,
  intercept_row = T,
  intercept_col = T,
  #warm_start = old_fit,
  trace = T,
  thresh = 1e-6,
  maxit = 300
)



IMR::get_lambda_lasso_max(dat$fit_data$train, NULL,
                          dat$fit_data$Z$Q, dat$fit_data$valid,
                          dat$fit_data$Y_full, T, T, verbose = T
                          )


future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()
hpar$M$lambda_factor <- 1
hpar$M$n.lambda <- 20
hpar$M$rank.init <- 2
hpar$M$rank.step <- 2
hpar$M$rank.max  <- 15
hpar$M$lambda_max <- 80
hpar$beta$n.lambda <- 10
hpar$gamma$n.lambda <- 10
fit32 <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                    Z = dat$fit_data$Z$Q,intercept_row = F,
                    hpar = hpar, seed = 2025,
                    intercept_col = F, verbose=2)
quick_camc_simu_res(dat, fit32$fit)



fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                    trace=T,tol = 5, n.lambda=80, rank.init = 2,
                    rank.step = 1
                    #lambda0_fun = IMR::get_lambda_m_max
                    )
quick_camc_simu_res(dat, fitsi$fit)




fitmmci <- MCCI.cv(dat$Y, dat$X, dat$mask, numCores = 9)
quick_camc_simu_res(dat, fitmmci$fit, T)

#-------------------------------------------------------
#---- testing IMC no bad local minima 2022
library(reticulate)
reticulate::py_run_string("import sys; sys.path.insert(0, r'~/Research/IMR/other_models/')")
np <- reticulate::import("numpy")
sp    <- import("scipy.sparse")
reticulate::source_python("./other_models/GNIMC.py")
GNIMC <- reticulate::py$GNIMC

pyopt <- list(
  X = dat$Y,
  omega = dat$mask,#sp$csr_matrix(dat$mask, dtype = np$float64),
  A = dat$X,
  B = dat$Z,
  verbose = TRUE,
  init_option = py$INIT_WITH_SVD
)
out <- do.call(GNIMC, c(list(rank=10L), pyopt))
estimate <- out[[1]]

xnorm <- np$linalg$norm(dat$theta[dat$Y==0])
np$linalg$norm(estimate[dat$Y==0] - dat$theta[dat$Y==0])/xnorm

IMR::error_metric$rmse(estimate[dat$Y==0], dat$theta[dat$Y==0])

#---------------------------------------------------------
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
init <- svd_opt(IMR:::naive_MC(y_train), 2, nr, nc, FALSE, FALSE)

mat <- as.matrix(y_train)


microbenchmark::microbenchmark(
a = partial_crossprod(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
b =  partial_crossprod_cpp(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
times = 500
)
