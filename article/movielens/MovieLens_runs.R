require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")

data <- load_movielens1m()
seed = 2025

model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0.2)
print(model_data)

model_combn <- data.frame(
  row_covariates = c(F,T,T),
  col_covariates = c(F,F,T),
  intercepts = c(T,T,T)
)
#-----


convergence <- imr_convergence(maxit=600, thresh=1e-5)
convergence2 <- imr_convergence(maxit=5000, thresh=1e-7)


grid <- imr_tune_grid(beta = c(0, 0.4, 60),
                      gamma = c(0, 0.4, 60),
                      laplace = c(0, 20, 80, 3),
                      rank = c(5, 20, 1, 3))

for(model_id in seq_along(model_combn)[c(1,3)]){

row_covariates <- model_combn$row_covariates[model_id]
col_covariates <- model_combn$col_covariates[model_id]

model_data <- update(model_data,
                     row_intercept = TRUE,
                     col_intercept = TRUE,
                     row_covariates = row_covariates,
                     col_covariates = col_covariates); model_data

# notes: max for model IMR-I: laplace = 120 but make it 45-16
#       max for model IMR-IX: beta = 2,  laplace = 135
#       max for model IMR-IXZ: beta: 1.9, laplace = 109, gamma = 4
# grid <- imr_set_grid_limits(model_data, grid,default_rank = 10,
#                             bisection_iter = 5,
#                             convergence=convergence, verbose=2)


fitimr <- IMR::imr_tune(model_data, grid, final_fit = FALSE,
                                          fast_laplace = FALSE,
                                          laplace_log_scale = FALSE,
                                          convergence=convergence, n_cores=7,
                                          seed = seed, verbose=1)

saveRDS(fitimr, paste0("article/movielens/data/saved_models/March_IMR_I",
               ifelse(row_covariates, "X",""),
               ifelse(col_covariates, "Z", ""),
               "_tune.rds"))

start <- Sys.time()
fitimr_fit <- IMR::imr_fit(model_data,
                       rank = fitimr$params$rank,
                       lambda_m = fitimr$params$lambda_laplace,
                       lambda_beta = fitimr$params$lambda_beta,
                       lambda_gamma = fitimr$params$lambda_gamma,
                        convergence=convergence2)
time <- Sys.time() - start
fitimr_fit$time_secs <- as.numeric(time, "secs")
saveRDS(fitimr_fit, paste0("article/movielens/data/saved_models/March_IMR_I",
                       ifelse(row_covariates, "X",""),
                       ifelse(col_covariates, "Z", ""),
                       "_fit_1e7.rds"))

print(fitimr_fit)
print(summary(fitimr_fit))
datrec <- reconstruct(fitimr_fit, model_data)

print(prepare_output_movielens(
  "IMR-...",
  time = fitimr_fit$time_secs,
  X = data$X,
  Z = data$Z,
  beta.estim = datrec$beta,
  gamma.estim = datrec$gamma,
  estim.test = datrec$estimates[data$test.idx],
  estim.train = as.Incomplete(datrec$estimates * data$obs_mask)@x,
  obs.test = data$test.truths,
  test_error = IMR:::error_metrics$rmse,
  obs.train = data$Y[data$Y != 0],
  M.estim = datrec$M,
  rank.M = fitimr$params$rank
))

}


fitimr <- readRDS("article/movielens/data/saved_models/March_IMR_I_tune.rds")
