library(devtools)
clean_dll()
Rcpp::compileAttributes()
document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article_results/simulation/generate_simu_dat.R")

sim_res <- function(dat, fit, name = "", ortho = TRUE,
                    error_metric = IMR:::error_metric$rmse,
                    digits = 5) {
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  refit <- IMR::reconstruct(fit, dat)
  test.obs <- dat$Y == 0

  out <- data.frame(model = name)

  out$theta_rmse <- error_metric(refit$estimates, dat$theta)
  out$test_rmse <- error_metric(refit$estimates[test.obs], dat$theta[test.obs])
  out$rank <- qr(refit$estimates)$rank
  true_singular <- IMR::opt_svd(dat$theta, tol = 1e-6)$d
  num_singular <- length(true_singular)
  if (length(fit$d) > num_singular) {
    estim_singular <- fit$d[1:num_singular]
  } else if (length(fit$d) < num_singular) {
    estim_singular <- c(fit$d, rep(0, num_singular - length(fit$d)))
  } else {
    estim_singular <- fit$d
  }
  out$eigen_rmse <- IMR:::error_metric$rmse(estim_singular, true_singular)
  out <- out %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
  return(out)
}


#-------------------------------------------------

n <- 800
m <- 900
r <- 6
seed <- 2025

dat <-
  generate_simulated_data(n, m, r, 0, 0, 0.8,
    sparsity_beta = .5, sparsity_gamma = 0.0,
    structured_error_A = T,
    structured_error_B = T,
    prepare_for_fitting = T, mv_coeffs = T, seed = 2025
  )

dat$similarity_rows %<>% solve()
dat$similarity_cols %<>% solve()

inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
data <- inp.dat$model

# the following are input parameters to the cv function >>
lambda_beta <- lambda_gamma <- 0
intercept_row <- intercept_col <- FALSE
error_function <- IMR:::error_metric$rmse
n_streaks <- 2; thresh = 1e-6; maxit=500; trace=TRUE; old_fit=NULL; ls_initial=FALSE;seed=2025

# set-up the hpar >>
hpar <- get_imr_default_hparams(dat$similarity_rows, data$similarity_cols, 0, 0)
hpar$laplace$lambda_step_sizes <- c(5, 1)
hpar$laplace$alpha_step_sizes <- c(0.1)
hpar$laplace$alpha_min <- hpar$laplace$alpha_max <- 0.5
hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1

# initialize parallel workers
initialize_parallel_workers(9)

#--- done -------- go run from the function >>
# delete this part later



#--------------


#---- we now test the cv lambda_laplace function: >>>>
bench::bench_time(fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                                      trace = T,
                                      tol = 2, seed = seed,
                                      maxit = 600,
                                      test_error = IMR:::error_metric$rmse,
                                      n.lambda = 20)) -> time.si

round(lubridate::time_length(time.si, "minute"), 2)

bench::bench_time(res <-
                    IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
                                         intercept_col = F,
                             num_cores = 9)) -> time.lap
round(lubridate::time_length(time.lap, "minute"), 2)
sim_res(dat, res$best_fit, "laplace")
sim_res(dat, fitsi$fit, "SoftImpute")


#--------------------------------------------------------
res <- data.frame()
for (b in 1:1) {
  timesim <- system.time(fitsim <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar, trace = T))
  timesi <- system.time(fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
    trace = T,
    tol = 2, seed = seed,
    n.lambda = 20,
    test_error = IMR:::error_metric$rmse,
    #rank.init = hpar$M$rank.min, maxit = 600,
    #rank.step = 1, rank.limit = hpar$M$rank.max
    # lambda0_fun = IMR::get_lambda_M_max
  ))


  timer <- system.time(fit.imrS <- IMR::imr.cv(inp.dat$model,
    intercept_row = F,
    seed = seed, ls_initial = FALSE, hpar = hpar,
    maxit = 600,
    intercept_col = F, verbose = 0
  ))

  res <- bind_rows(res, rbind(
    sim_res(dat, fitsi$fit, "SoftImpute") %>%
      dplyr::mutate(valerr = fitsi$error, time = timesi[[3]]),
    sim_res(dat, fit.imrS$fit, "IMR") %>%
      dplyr::mutate(valerr = fit.imrS$error, time = timer[[3]]),
    sim_res(dat, fitsim$fit, "IMR-Similarity") %>%
      dplyr::mutate(valerr = fitsim$cols$best_error, time = timesim[[3]])
  ))
  print(res)
}



fitsi$fit$d
fit.imrS$fit$d
fitsim$fit$d
svd(dat$theta)$d[1:r]

fit.imrS$lambda_M
fitsim$cols$best_parameter
fitsim$rows$best_parameter




# goal: compare cv_laplace - cv_M - Soft-impute > no  covariates but with structured error
# we want to follow what happens inside
lambda_betaa <- 0
lambda_gammaa <- 0
intercept_row <- FALSE
intercept_col <- FALSE
trace <- T
thresh <- 1e-6
maxit <- 300
ls_initial <- FALSE
