library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article/simulation/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")
source("./other_models/MCCI.R")
#==============================================
sim1_res <- function(dat, fit, name="",
                     error_metric=IMR:::error_metric$rel.rmse,
                     coeffs_transformed = TRUE){
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  has.beta <- "beta" %in% names(fit) & ! is.null(fit$beta)
  has.gamma <- "gamma" %in% names(fit) & ! is.null(fit$gamma)
  has.M    <- "M" %in% names(fit) | all(c("u","d","v") %in% names(fit))
  has.intercept <- "beta0" %in% names(fit) & !is.null(fit$beta0)
  estimates <- 0
  out <- data.frame(model=name, beta=NA, gamma=NA, M=NA, theta=NA, test=NA, rank=NA)
  # check M
  if(all(c("u", "d", "v") %in% names(fit)))
    fit$M <- fit$u %*% (fit$d * t(fit$v))
  if(has.M) {
    estimates <- fit$M
    out$M <- error_metric(fit$M, dat$M)
  }
  if(has.beta){
    if(coeffs_transformed)
      fit$beta <- solve(dat$Xr) %*% fit$beta

    estimates <- estimates + dat$X %*% fit$beta
    out$beta <- error_metric(fit$beta, dat$beta)
  }
  if(has.gamma){
    if(coeffs_transformed)
      fit$gamma <- fit$gamma %*% solve(t(dat$Zr))
    estimates <- estimates + fit$gamma %*% t(dat$Z)
    out$gamma <- error_metric(fit$gamma, dat$gamma)
  }
  if(has.intercept){
    estimates <- estimates + fit$beta0 %*% matrix(1,1,ncol(dat$Y))
  }


  stopifnot(all(estimates!=0))
  out$theta <- error_metric(estimates, dat$theta)
  test.obs <- dat$Y == 0
  out$test_rel <- IMR:::error_metric$rel.rmse(estimates[test.obs], dat$theta[test.obs])
  out$test <- error_metric(estimates[test.obs], dat$theta[test.obs])
  train.obs <- dat$Y != 0
  out$train <- error_metric(estimates[train.obs], dat$theta[train.obs])
  out$rank <- qr(estimates)$rank
  return(out)
}

#===========================================================================================
# setting 2)
n = m = 400
p = 5;
q = 5;
r = 4;
missing_pct = seq(.7, .98, .2)
models <- c("IMR", "SImpute")
all_res <- res <- data.frame()
for(b in 1:500){
  seed = 2025 + b
  set.seed(seed)
for(pct in missing_pct){

dat <-
  generate_simulated_data(n, m, r, p, q, pct,
                          sparsity_beta = 0,
                          sparsity_gamma = 0,
                          intercept = FALSE,
                          structured_error_A = F, SNR = 1,
                          structured_error_B = F,
                          prepare_for_fitting = F, mv_coeffs = T, seed = seed)


mdat <- IMR::prepare_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2)

fitsi <- simpute.cv(y_full = mdat$Y,
                    y_train = mdat$y_train,
                    y_valid = mdat$y_valid,
                    trace = FALSE,
                    print.best = FALSE,
                    tol = 2,
                    maxit = 1000,
                    thresh = 1e-6,
                    n.lambda = 20,
                    test_error = IMR:::error_metric$rmse,
                    seed = seed)

hparam <- IMR::get_imr_default_hparams()
hparam$laplace$step_sizes = c(1)#c(5,1,0.1)
hparam$laplace$min = 10
hparam$laplace$max = 30

hparam$rank$max = 10
hparam$beta$length = 10
hparam$gamma$length = 10
hparam$beta$max <- hparam$gamma$max <- 2
#hparam$beta$max = 0

# first tune with a single laplace value
hparam$laplace$min = hparam$laplace$max = 18.3; hparam$laplace$step_sizes = c(1)
fitimr = IMR:::imr.cv_2(mdat, intercept_row = FALSE, intercept_col = FALSE,
                        hpar = hparam, thresh = 1e-6, maxit = 1000,
                        trace = 3, ls_initial = TRUE, shared_information = FALSE,
                        seed = seed, num_cores = 7,
)
# then we refit to re-tune lammbda_laplace



bench::mark(

  fitimr = IMR:::imr.cv_2(mdat, intercept_row = FALSE, intercept_col = FALSE,
                          hpar = hparam, thresh = 1e-6, maxit = 1000,
                          trace = 3, ls_initial = TRUE, shared_information = FALSE,
                          seed = seed, num_cores = 5,
  ),
  # fitimr = IMR:::imr.cv_3(mdat, intercept_row = FALSE, intercept_col = FALSE,
  #                         hpar = hparam, thresh = 1e-6, maxit = 1000,
  #                         trace = 3, ls_initial = TRUE, shared_information = FALSE,
  #                         seed = seed, num_cores = 7,
  # ),
iterations = 1,
check = T

) -> bb



bb$result[[1]]$fit$params

print(bb)

dat$Xr <- mdat$Xr
dat$Zr <- mdat$Zr
errorm <- IMR:::error_metric$rmse

sim1_res(dat, bb$result[[1]]$fit, "IMR", errorm)
sim1_res(dat, fitimr$fit, "IMR", errorm)


bb2$result[[1]]$hpar$beta

fitimr$fit$params

fitimr$fit$beta

b1 <- bb$result[[1]]$fit$beta
b2 <- bb$result[[1]]$fit2$beta
b3 <- bb2$result[[1]]$fit$beta

mask = (b1 != 0) * 1



fit4 <- IMR:::imr.fit_debias(
  Y = mdat$Y,
  X = mdat$Xq,
  Z = mdat$Zq,
  r = 5,
  lambda_M = 14.1,
  mask_beta = mask,
  #lambda_gamma = 0,
  intercept_row = F,
  intercept_col = F,
  shared_information = F,
  warm_start = bb$result[[1]]$fit,
  trace = T,
  thresh = 1e-6,
  maxit = 1000,
  ls_initial = F
)


b1[1:5,1:5]
b2[1:5,1:5]
b3[1:5,1:5]


# hparam$beta$step_sizes <- hparam$gamma$step_sizes <-
#   step_size <- function(min_val, max_val, n=20) {
#     step = (max_val - min_val) / (n - 1L)
#     print(step)
#     return(step)
#   }
#
# IMR::initialize_parallel_workers(9)
# fitimro <- IMR::imr.cv(mdat, intercept_row = FALSE, intercept_col = FALSE,
#                       hpar = hparam, thresh = 1e-5, maxit = 1000,
#                       trace = 2, ls_initial = TRUE, shared_information = FALSE,
#                       seed = 2025, num_cores = 9,
#                       )


fitmcci <- MCCI.cv(Y = dat$Y, X = dat$X, W = dat$mask, n_folds = 5,numCores = 9,
                   seed = seed,
                   test_error = IMR:::error_metric$rmse,
                   lambda_1_grid = c(0),#seq(0, 1, length = 10),
                   lambda_2_grid = seq(2.9, 0.1, length = 18),
                   alpha_grid = c(1),#seq(0.992, 1, length = 10),
                   n1n2_optimized = TRUE,
                   return_diagn = FALSE)

dat$Xr <- mdat$Xr
errorm <- IMR:::error_metric$rmse
sim1_res(dat, fitsi$fit, "SI", errorm) %>% rbind(
  sim1_res(dat, fitmcci$fit, "MCCI", errorm, coeffs_transformed = FALSE),
  sim1_res(dat, fitimr$fit, "IMR", errorm)
) -> res
errorm <- IMR:::error_metric$rel.rmse
res %<>% rbind(sim1_res(dat, fitsi$fit, "SI", errorm) %>% rbind(
  sim1_res(dat, fitmcci$fit, "MCCI", errorm, coeffs_transformed = FALSE),
  sim1_res(dat, fitimr$fit, "IMR", errorm)
))
res$dim = d
all_res %<>% rbind(res)
print(d)
}
  print(b)
  print(res)
saveRDS(all_res, "./article/simulation/data/sim2_res.rds")
}




errorm <- IMR:::error_metric$rmse
sim1_res(dat, fitsi$fit, "SI", errorm) %>% rbind(
sim1_res(dat, fitmcci$fit, "MCCI", errorm, coeffs_transformed = FALSE),
sim1_res(dat, fitimr2$fit, "IMR2", errorm),
sim1_res(dat, fitimr2$fit2, "IMR22", errorm),
sim1_res(dat, fitimro$best_fit, "IMRo", errorm)
#sim1_res(dat, fitimr$fit, "IMR", errorm)
) %>%
  arrange(test)

fitimr$fit$params

fitimr$best_fit$r

fitmcci$best_parameters
fitmcci$best_score
fitimr$best_fit$
#======== delete start ====
hparam$beta$length = 10;
hparam$gamma$length = 1; hparam$gamma$max = 0;


data <- mdat; intercept_row = FALSE; intercept_col = FALSE;
hpar = hparam; error_function = IMR:::error_metric$rmse;
thresh = 1e-4; maxit = 300; trace = 3; ls_initial = TRUE;
shared_information = FALSE; num_cores = 9; warm_start=NULL; seed = 2025;

lambda_beta = 0.2; lambda_gamma = 0.2; lambda_laplace = 1

#==== delete end =====
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

