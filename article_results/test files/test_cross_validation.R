library(devtools)
clean_dll()
Rcpp::compileAttributes()
document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article_results/simulation/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")

sim_res <- function(dat, fit, name = "",
                    error_metric = IMR:::error_metric$rmse,
                    digits = 5) {
  # dat should include: Y, theta, beta, X, Xr, Z, Zr
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  refit <- IMR::reconstruct(fit, dat)
  test.obs <- dat$Y == 0

  out <- data.frame(model = name)

  out$theta_rmse <- error_metric(refit$estimates, dat$theta)
  out$test_rmse <- error_metric(refit$estimates[test.obs], dat$theta[test.obs])
  out$beta_rmse <- error_metric(refit$beta, dat$beta)
  out$M_rmse <- error_metric(refit$M, dat$M)

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
m <- 600
r <- 6
seed <- 2025

dat <-
  generate_simulated_data(n, m, r, 3, 5, 0.8,
    sparsity_beta = 0, sparsity_gamma = 0,
    intercept = FALSE,
    structured_error_A = T,
    structured_error_B = T,
    prepare_for_fitting = T, mv_coeffs = T, seed = 2025
  )

dat$similarity_rows %<>% solve()
dat$similarity_cols %<>% solve()


# dat$similarity_rows <- diag(1, nrow(dat$theta))
# dat$similarity_cols <- diag(1, ncol(dat$theta))

inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
data <- inp.dat
dat$Xr <- data$Xr;


# the following are input parameters to the cv function >>
lambda_beta <- lambda_gamma <- 0
intercept_row <- intercept_col <- FALSE
error_function <- IMR:::error_metric$rmse
n_streaks <- 2; thresh = 1e-6; maxit=500; trace=2;
old_fit=NULL; ls_initial=TRUE;seed=2025
num_cores = 0; warm_start = NULL;

# set-up the hpar >>
hpar <- get_imr_default_hparams(data$similarity_rows, data$similarity_cols, 0, 0)
hpar$laplace$lambda_step_sizes <- c(.2)
hpar$laplace$alpha_step_sizes <- c(0.1)
hpar$laplace$lambda_max <- 10
hpar$laplace$alpha_min <- 0.5; hpar$laplace$alpha_max <-0.5
hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1
hpar$beta$step_sizes <- c(0.3, 0.1, 0.05)
hpar$beta$n_streaks <- 2
# hpar$beta$max <- hpar$beta$value <- 500
# hpar$gamma$max <- hpar$gamma$value <- 500

# initialize parallel workers
IMR::initialize_parallel_workers(9)

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

data_nocov <- data
data_nocov$X <- data_nocov$Xq <- data_nocov$Xr <- NULL


#-  1. IMR without intercepts or covariates
res1 <-
  IMR:::imr.cv_laplace(data_nocov, trace=2, hpar=hpar, intercept_row = F,
                       intercept_col = F, ls_initial = F,
                       seed = 2025, warm_start = NULL,
                       num_cores = 0)

#-  2. IMR with intercept
res2 <-
  IMR:::imr.cv_laplace(data_nocov, trace=2, hpar=hpar, intercept_row = T,
                       intercept_col = T, ls_initial = F,
                       seed = 2025, warm_start = NULL,
                       num_cores = 0)

#-  3.  IMR with covariates
res3 <-
  IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
                       intercept_col = F, ls_initial = F,
                       seed = 2025, warm_start = NULL,
                       num_cores = 0)


#-  4. IMR with covariates but initialize covariates with IMR with covariates
warm <- res3$best_fit
resid <- data$Y
resid@x <- resid@x - IMR:::partial_crossprod(dat$X, warm$beta, Y2@i, Y2@p)
init <- opt_svd(naive_MC(resid), 2, n, m, FALSE, FALSE)
warm$u <- init$u; warm$v = init$v; warm$d = init$d;
res4 <-
  IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
                       intercept_col = F, ls_initial = F,
                       seed = 2025, warm_start = warm,
                       num_cores = 0)


#- 5. IMR with covariates and intercept
# res5 <-
#   IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
#                        intercept_col = F, ls_initial = F,
#                        seed = 2025, warm_start = NULL,
#                        num_cores = 0)

#- 6. IMR with covariates but initialize covariates to TRUE value
warm <- res3$best_fit
warm$beta <- dat$beta
resid <- data$Y
resid@x <- resid@x - IMR:::partial_crossprod(dat$X, warm$beta, Y2@i, Y2@p)
init <- opt_svd(naive_MC(resid), 2, n, m, FALSE, FALSE)
warm$u <- init$u; warm$v = init$v; warm$d = init$d;
res6 <-
  IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
                       intercept_col = F, ls_initial = F,
                       seed = 2025, warm_start = warm,
                       num_cores = 0)


#- 7. IMR with covariates but initialize AB  to IMR without covariates
warm <- res3$best_fit
warm$beta <- matrix(0, nrow(warm$beta), ncol(warm$beta))
res7 <-
  IMR:::imr.cv_laplace(data, trace=2, hpar=hpar, intercept_row = F,
                       intercept_col = F, ls_initial = F,
                       seed = 2025, warm_start = warm,
                       num_cores = 0)



errormet <- IMR:::error_metric$rmse
require(gt)

rbind(
  sim_res(dat, res1$best_fit, "None",errormet),
  sim_res(dat, res2$best_fit, "intercept",errormet ),
  sim_res(dat, res3$best_fit, "Covariates",errormet),
  sim_res(dat, res4$best_fit, "Covariates (beta from Covariates)",errormet),
  #sim_res(dat, res5$best_fit, "Covariates (but with random init)",errormet),
  sim_res(dat, res6$best_fit, "Covariates(beta is true beta)",errormet),
  sim_res(dat, res7$best_fit, "Covariates(AB^T is from None)",errormet),
  sim_res(dat, fitsi$fit, "SoftImpute",errormet)
) %>%
  arrange(test_rmse) %>%
  dplyr::select(-eigen_rmse) %>%
  gt(rowname_col = "model") |>
  fmt_number(
    columns = c(theta_rmse, test_rmse, beta_rmse, M_rmse),
    decimals = 3
  ) |>
  cols_label(
    theta_rmse = html("RMSE(&theta;)"),
    test_rmse  = "Test RMSE",
    beta_rmse  = html("RMSE(&beta;)"),
    M_rmse     = html("RMSE(M)"),
    rank       = "Rank (true=6)"
   # eigen_rmse = "Eigen RMSE"
  ) |>
  tab_header(
    title = "Model comparison based on RMSE metrics"
  ) |>
  opt_table_outline() |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  )

# inp.dat <- IMR::prepare_data(dat$Y, NULL, NULL, dat$similarity_rows, dat$similarity_cols)
# data <- inp.dat$model

#
# Y2 <- IMR::as.Incomplete(dat$Y)
# Y2@x <- Y2@x - IMR:::partial_crossprod(dat$X, res1$best_fit$beta, Y2@i, Y2@p)
#
# data2 <- IMR::prepare_data(as.matrix(Y2), dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
# warm <- res1$best_fit
# resid <- data$Y
# resid@x <- resid@x - IMR:::partial_crossprod(dat$X, warm$beta, Y2@i, Y2@p)
# init <- opt_svd(naive_MC(resid), 2, n, m, FALSE, FALSE)
# warm$u <- init$u; warm$v = init$v; warm$d = init$d;
#
# bench::bench_time(res3 <-
#                     IMR:::imr.cv_laplace(data, trace=2, hpar=hpar2, intercept_row = F,
#                                          intercept_col = F, ls_initial = T,
#                                          seed = 2025, warm_start = warm,
#
#                              num_cores = 0)) -> time.lap
# round(lubridate::time_length(time.lap, "minute"), 2)
# res3$best_fit$beta <- res1$best_fit$beta
#
# hpar2 <- hpar;
# hpar2$beta$value <- 0;
#
# bench::bench_time(res1 <-
#                     IMR:::imr.cv_laplace(data, trace=2, hpar=hpar2, intercept_row = F,
#                                          intercept_col = F, ls_initial = T,
#                                          seed = 2025,
#
#                                          num_cores = 0)) -> time.lap
#
# bench::bench_time(res2 <-
#                     IMR:::imr.cv_laplace(data, trace=2, hpar=hpar2, intercept_row = T,
#                                          intercept_col = T, ls_initial = T,
#                                          seed = 2025,
#
#                                          num_cores = 0)) -> time.lap
#
#
#
#
# bench::bench_time(res1 <-
#                     IMR:::imr.cv(data, trace=2, hpar=hpar, intercept_row = F,
#                                          intercept_col = F,
#                                  warm_start = NULL,
#                                  ls_initial = T,
#                                          num_cores = 0)) -> time.lap1
# round(lubridate::time_length(time.lap1, "minute"), 2)


errormet <- IMR:::error_metric$rmse

rbind(
sim_res(dat, res1$best_fit, "laplace+cov",errormet),
sim_res(dat, res3$best_fit, "laplace+cov+int",errormet ),
sim_res(dat, res$best_fit, "laplace",errormet),
sim_res(dat, fitsi$fit, "SoftImpute",errormet)
) %>%
  arrange(test_rmse)

res1$best_fit$beta0 %>% summary()
colMeans(dat$Y) %>% summary()
res1$best_fit$gamma0 %>% summary()
rowMeans(dat$Y) %>% summary()

hpar$laplace$lambda_min <-
  hpar$laplace$lambda_max <- res$best_fit$lambda_laplace
hpar$laplace$alpha_min <-
  hpar$laplace$alpha_max <- res$best_fit$alpha
hpar$rank$min <- 2; hpar$rank$max <- 10;# res$best_fit$r
hpar$beta$max <- NULL; hpar$beta$value <- 1; hpar$beta$min <- 0;



#-- detele this later

#--- end delete


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
