library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article_results/simulation/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")



#' 1. Does randomization alter the outcome? also test the initial values >>
#' Data: row and column covariates, one setting, no similarity, shared information; same dataset
#' Models: IMR-XZ-Random, IMR-XZ-LS, LS
n = m = 800
r = 6;
p=q=2; missp = 0.8

dat <-
  generate_simulated_data(n, m, r, p, q, missp,
                          sparsity_beta = 0, sparsity_gamma = 0,
                          intercept = FALSE,
                          structured_error_A = F, SNR = 1,
                          structured_error_B = F,
                          prepare_for_fitting = T, mv_coeffs = F, seed = 2025
  )
mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols,seed = 2025)
dat$Xr <- mdata$Xr;
dat$Zr <- mdata$Zr;
mdata$X <- mdata$Z <- NULL
#-----------------------------------------------------------------------
hparam <- get_imr_default_hparams(mdata$similarity_rows, mdata$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(1, .1)
hparam$laplace$min <- 0;
hparam$laplace$max <- 17
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
#-------------------------------------------------------------------------------------
results1_1 <- list(IMRXZR = list(), IMRXZLS = list(), LS = list())
results1_1$DATA <- list(dat=dat, data=mdata, hparam=hparam)
#--------------------------------------------------------------------------------------
initialize_parallel_workers(9)
# 20 replications
B = 20
for(i in 1:B){
  seed = 2025 + i
  set.seed(seed)
  print(i)
#--------------------------------------------------------------------------------------
  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols,seed = 2025)
  bench::bench_time(results1_1$IMRXZR[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                trace=1, hpar=hparam, intercept_row = T,
                                                intercept_col = T, ls_initial = F,
                                                seed = seed, warm_start = NULL, maxit=600,
                                                shared_information = T,
                                                num_cores = 0)) -> timee
results1_1$IMRXZR[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)

mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols,seed = 2025)
bench::bench_time(results1_1$IMRXZLS[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                             trace=1, hpar=hparam, intercept_row = T,
                                                             intercept_col = T, ls_initial = T,
                                                             maxit=600,
                                                             seed = seed, warm_start = NULL,
                                                             shared_information = T,
                                                             num_cores = 0)) -> timee

results1_1$IMRXZLS[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)
  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols,seed = 2025)

bench::bench_time(results1_1$LS[[i]] <- IMR::imr.fit_no_low_rank(
  Y = mdata$y_train, X = mdata$Xq, Z = mdata$Zq,
  lambda_beta = 0, lambda_gamma = 0, intercept_row = F,
  intercept_col = F, shared_information = T, maxit = 600,
  trace = F
)) -> timee
sdd = IMR::opt_svd(IMR::naive_MC(results1_1$LS[[i]]$resid), r,n,m,F,F)
results1_1$LS[[i]][names(sdd)] <- sdd
results1_1$LS[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)

if(i %% 4 == 0 || i == 1){
  print(rbind(sim_res(dat, results1_1$IMRXZR[[i]]$best_fit, "XZ",T),
              sim_res(dat, results1_1$IMRXZLS[[i]]$best_fit, "None",T),
              sim_res(dat, results1_1$LS[[i]], "LS",T)))

}
 saveRDS(results1_1, "article_results/simulation/data/december_2025_results1_1.rds")
}
#-------------------------------------------------------------------------------
#sim_res(dat, results1$IMRXZLS[[20]]$best_fit, "Covariates (beta from Covariates)",errormet,T)
#---------------------------------------------------------------------------------
#' scanrio 1:
#' Data:  with X and Z. no correlation
#' Models: 1. IMR, 2. IMR-XZ
#-------------------------------------------------------------------------------------

# 20 replications
B = 20
results2_1 <- list(IMR = list(), IMRXZ = list(), data=list())
for(i in 1:B){
seed = 2025 + i
set.seed(seed)
print(i)
#--------------------------------------------------------------------------------------
dat <-
  generate_simulated_data(n, m, r, q, q, missp,
                          sparsity_beta = 0, sparsity_gamma = 0,
                          intercept = FALSE,
                          structured_error_A = F, SNR = 1,
                          structured_error_B = F,
                          prepare_for_fitting = T, mv_coeffs = F, seed = seed
  )
dat$Xr <- dat$fit_data$X$R;
dat$Zr <- dat$fit_data$Z$R;
mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
mdata$X <- mdata$Z <- NULL

data_nocov <- mdata
data_nocov$X <- data_nocov$Xq <- data_nocov$Xr <- NULL
data_nocov$Z <- data_nocov$Zq <- data_nocov$Zr <- NULL

#-----------------------------------------------------------------------
hparam <- get_imr_default_hparams(mdata$similarity_rows, mdata$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(1, .1)
hparam$laplace$max <- 17
hparam$laplace$min <- 0;
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
#-------------------------------------------------------------------------------------
results2_1$data[[i]] <- dat
#--------------------------------------------------------------------------------------
  bench::bench_time(results2_1$IMR[[i]] <- IMR:::imr.cv_laplace(data_nocov,
                                                                 trace=1, hpar=hparam,
                                                                intercept_row = T,
                                                                 intercept_col = T, ls_initial = T,
                                                                 seed = seed, warm_start = NULL,
                                                                  maxit=600,
                                                                 shared_information = TRUE,
                                                                 num_cores = 9)) -> timee
  results2_1$IMR[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)

  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  mdata$X <- mdata$Z <- NULL
  bench::bench_time(results2_1$IMRXZ[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                                  trace=1, hpar=hparam,
                                                                  intercept_row = T,
                                                                  intercept_col = T, ls_initial = T,
                                                                  maxit=600,
                                                                  seed = seed, warm_start = NULL,
                                                                  shared_information = TRUE,
                                                                  num_cores = 9)) -> timee

  results2_1$IMRXZ[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)
  if(i %% 4 == 0 || i == 1){
    print(rbind(sim_res(dat, results2_1$IMRXZ[[i]]$best_fit, "XZ"),
                sim_res(dat, results2_1$IMR[[i]]$best_fit, "None")))

  }
  saveRDS(results2_1, "article_results/simulation/data/december_2025_results2_1.rds")
}


#---------------------------------------------------------------------------------
#' scanrio 2:
#' Data:  with X and Z. row correlation. no column correlation
#' Models: 1. IMR, 2. IMR-Laplace,  IMR-XZ, IMR-XZ-Laplace
#-------------------------------------------------------------------------------------

# 20 replications
B = 20
results3 <- list(IMR = list(), IMRXZ = list(), data=list())
for(i in 1:B){
  seed = 2025 + i
  set.seed(seed)
  print(i)
  #--------------------------------------------------------------------------------------
  dat <-
    generate_simulated_data(n, m, r, q, q, missp,
                            sparsity_beta = 0, sparsity_gamma = 0,
                            intercept = FALSE,
                            structured_error_A = T, SNR = 1,
                            structured_error_B = F,
                            prepare_for_fitting = T, mv_coeffs = F, seed = seed
    )
  dat$Xr <- dat$fit_data$X$R;
  dat$Zr <- dat$fit_data$Z$R;
  dat$similarity_rows %<>% solve()

  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  mdata$X <- mdata$Z <- NULL

  data_nocov <- mdata
  data_nocov$X <- data_nocov$Xq <- data_nocov$Xr <- NULL
  data_nocov$Z <- data_nocov$Zq <- data_nocov$Zr <- NULL

  #-----------------------------------------------------------------------
  hparam <- IMR::get_imr_default_hpar(mdata$similarity_rows, mdata$similarity_cols, 0, 0)
  hparam$laplace$step_sizes <- c(1, .1)
  hparam$laplace$max <- 17
  hparam$laplace$min <- 0;
  hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
  #-------------------------------------------------------------------------------------
  results3$data[[i]] <- dat
  #--------------------------------------------------------------------------------------
  bench::bench_time(results3$IMRL[[i]] <- IMR:::imr.cv_laplace(data_nocov,
                                                              trace=1, hpar=hparam, intercept_row = F,
                                                              intercept_col = F, ls_initial = T,
                                                              seed = seed, warm_start = NULL,
                                                              maxit=600,
                                                              shared_information = TRUE,
                                                              num_cores = 0)) -> timee
  results3$IMRL[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  mdata$X <- mdata$Z <- NULL
  bench::bench_time(results3$IMRXZL[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                                trace=1, hpar=hparam, intercept_row = F,
                                                                intercept_col = F, ls_initial = T,
                                                                maxit=600,
                                                                seed = seed, warm_start = NULL,
                                                                shared_information = TRUE,
                                                                num_cores = 0)) -> timee

  results3$IMRXZL[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  data_nocov <- IMR::prepare_data(dat$Y, NULL, NULL, NULL, NULL)
  bench::bench_time(results3$IMR[[i]] <- IMR:::imr.cv_laplace(data_nocov,
                                                               trace=1, hpar=hparam, intercept_row = F,
                                                               intercept_col = F, ls_initial = T,
                                                               seed = seed, warm_start = NULL,
                                                               maxit=600,
                                                               shared_information = TRUE,
                                                               num_cores = 0)) -> timee
  results3$IMR[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, NULL, NULL)
  mdata$X <- mdata$Z <- NULL
  bench::bench_time(results3$IMRXZ[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                                 trace=1, hpar=hparam, intercept_row = F,
                                                                 intercept_col = F, ls_initial = T,
                                                                 maxit=600,
                                                                 seed = seed, warm_start = NULL,
                                                                 shared_information = TRUE,
                                                                 num_cores = 0)) -> timee

  results3$IMRXZ[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  if(i %% 4 == 0 || i == 1){
    print(rbind(sim_res(dat, results3$IMRXZ[[i]]$best_fit, "XZ"),
                sim_res(dat, results3$IMR[[i]]$best_fit, "None"),
                sim_res(dat, results3$IMRXZL[[i]]$best_fit, "XZL"),
                sim_res(dat, results3$IMRL[[i]]$best_fit, "L")))

  }


  saveRDS(results3, "article_results/simulation/data/december_2025_results3.rds")
}




#---------------------------------------------------------------------------------
#' scanrio 3:
#' Data:  with X and Z. row and column correlation
#' Models: 1. IMR, 2. IMR-Laplace,  IMR-XZ, IMR-XZ-Laplace
#-------------------------------------------------------------------------------------


# 20 replications
B = 20
results4 <- list(IMR = list(), IMRXZ = list(), data=list())
for(i in 1:B){
  seed = 2025 + i
  set.seed(seed)
  print(i)
  #--------------------------------------------------------------------------------------
  dat <-
    generate_simulated_data(n, m, r, q, q, missp,
                            sparsity_beta = 0, sparsity_gamma = 0,
                            intercept = FALSE,
                            structured_error_A = T, SNR = 1,
                            structured_error_B = T,
                            prepare_for_fitting = T, mv_coeffs = F, seed = seed
    )
  dat$Xr <- dat$fit_data$X$R;
  dat$Zr <- dat$fit_data$Z$R;
  dat$similarity_rows %<>% solve()
  dat$similarity_cols %<>% solve()

  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  mdata$X <- mdata$Z <- NULL

  data_nocov <- mdata
  data_nocov$X <- data_nocov$Xq <- data_nocov$Xr <- NULL
  data_nocov$Z <- data_nocov$Zq <- data_nocov$Zr <- NULL

  #-----------------------------------------------------------------------
  hparam <- get_imr_default_hparams(mdata$similarity_rows, mdata$similarity_cols, 0, 0)
  hparam$laplace$step_sizes <- c(1, .1)
  hparam$laplace$max <- 17
  hparam$laplace$min <- 0;
  hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
  #-------------------------------------------------------------------------------------
  results4$data[[i]] <- dat
  #--------------------------------------------------------------------------------------
  bench::bench_time(results4$IMRL[[i]] <- IMR:::imr.cv_laplace(data_nocov,
                                                               trace=1, hpar=hparam, intercept_row = F,
                                                               intercept_col = F, ls_initial = T,
                                                               seed = seed, warm_start = NULL,
                                                               maxit=600,
                                                               shared_information = TRUE,
                                                               num_cores = 0)) -> timee
  results4$IMRL[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  mdata$X <- mdata$Z <- NULL
  bench::bench_time(results4$IMRXZL[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                                 trace=1, hpar=hparam, intercept_row = F,
                                                                 intercept_col = F, ls_initial = T,
                                                                 maxit=600,
                                                                 seed = seed, warm_start = NULL,
                                                                 shared_information = TRUE,
                                                                 num_cores = 0)) -> timee

  results4$IMRXZL[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  data_nocov <- IMR::prepare_data(dat$Y, NULL, NULL, NULL, NULL)
  bench::bench_time(results4$IMR[[i]] <- IMR:::imr.cv_laplace(data_nocov,
                                                              trace=1, hpar=hparam, intercept_row = F,
                                                              intercept_col = F, ls_initial = T,
                                                              seed = seed, warm_start = NULL,
                                                              maxit=600,
                                                              shared_information = TRUE,
                                                              num_cores = 0)) -> timee
  results4$IMR[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)


  mdata <- IMR::prepare_data(dat$Y, dat$X, dat$Z, NULL, NULL)
  mdata$X <- mdata$Z <- NULL
  bench::bench_time(results4$IMRXZ[[i]] <- IMR:::imr.cv_laplace(mdata,
                                                                trace=1, hpar=hparam, intercept_row = F,
                                                                intercept_col = F, ls_initial = T,
                                                                maxit=600,
                                                                seed = seed, warm_start = NULL,
                                                                shared_information = TRUE,
                                                                num_cores = 0)) -> timee

  results4$IMRXZ[[i]]$time <- round(lubridate::time_length(timee[2], "seconds"), 2)

  if(i %% 4 == 0 || i == 1){
    print(rbind(sim_res(dat, results4$IMRXZ[[i]]$best_fit, "XZ"),
                sim_res(dat, results4$IMR[[i]]$best_fit, "None"),
                sim_res(dat, results4$IMRXZL[[i]]$best_fit, "XZL"),
                sim_res(dat, results4$IMRL[[i]]$best_fit, "L")))

  }

  saveRDS(results4, "article_results/simulation/data/december_2025_results4.rds")
}

