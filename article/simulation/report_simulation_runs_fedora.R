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
#==================================================================
#' we now work on second part of the simulation
#' Results needed: Rmse > beta, gamma, M, test
#'                 obtain an average hyperparameters (r, 3 x lambda)
#' Variables: Fix n,m,p,q; change observation rate
#' Models : IMC, SoftImpute
#' Second step: use the fit function only to get speed vs model and observation rate
increase_sparsity <- function(dat, step=0.05){
  current_sparsity <- mean(dat$mask == 0)
  target_sparsity <- step + current_sparsity
  print(paste("Target sparsity is ", target_sparsity))
  stopifnot(target_sparsity < 1)
  extra_nonzero_frac <- step / (1 - current_sparsity)
  nonzero_idx <- which(dat$mask == 1)
  to_zero_ind <- sample(nonzero_idx, extra_nonzero_frac*length(nonzero_idx),replace = F)

  dat$Y[to_zero_ind] <- 0
  #dat$Y %<>% IMR::as.Incomplete()
  #-- we now recreate the train/test splits
  dat$mask <- as.matrix(dat$Y != 0)
  dat$sparsity <- target_sparsity
  return(dat)
}


# --=--=-=-=-=



n = m = 1000
p = 5;
q = 5;
r = 5;
missing_pct = seq(.7, .98, .05)
models <- c("IMR", "SImpute")
all_res <- res <- data.frame()
# b = 1; pct=1
for(b in 1:500){
  seed = 2025 + b
  start1 = Sys.time()
  set.seed(seed)
  dat <-
    generate_simulated_data(n, m, r, p, q, .7,
                            sparsity_beta = 0,
                            sparsity_gamma = 0,
                            intercept = FALSE,
                            structured_error_A = F, SNR = 1,
                            structured_error_B = F,
                            prepare_for_fitting = F, mv_coeffs = T, seed = seed)
  for(pct in 1:length(missing_pct)){
    start2 = Sys.time()
    if(pct > 1)
      dat <- increase_sparsity(dat, .05)

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


    hparam$rank$max = 10
    hparam$beta$length = 5
    hparam$gamma$length = 5
    hparam$beta$max <- hparam$gamma$max <- 0.5
    #hparam$beta$max = 0

    # first tune with a single laplace value
    hparam$laplace$min = hparam$laplace$max = 18.9; hparam$laplace$step_sizes = c(1)
    fitimr = IMR:::imr.cv_2(mdat, intercept_row = FALSE, intercept_col = FALSE,
                            hpar = hparam, thresh = 1e-6, maxit = 1000,
                            trace = 1, ls_initial = TRUE, shared_information = FALSE,
                            seed = seed, num_cores = 10)

    # then we refit to re-tune lammbda_laplace
    hparam$beta$value = fitimr$fit$params$lambda_beta
    hparam$gamma$value = fitimr$fit$params$lambda_gamma
    hparam$laplace$step_sizes = c(5,1,0.1)
    hparam$laplace$min = 10
    hparam$laplace$max = 25

    fitimr = IMR::imr.cv_laplace(mdat, intercept_row = FALSE, intercept_col = FALSE,
                                 hpar = hparam, thresh = 1e-6, maxit = 1000,
                                 trace = 1, ls_initial = TRUE, shared_information = FALSE,
                                 seed = seed, num_cores = 10,final_fit = TRUE,
                                 warm_start = fitimr$fit)



    fitimr$best_fit$lambda_laplace


    dat$Xr <- mdat$Xr
    dat$Zr <- mdat$Zr
    errorm <- IMR:::error_metric$rmse
    sim1_res(dat, fitsi$fit, "SI", errorm) %>%
      rbind(sim1_res(dat, fitimr$best_fit, "IMR", errorm)) %>%
      mutate(metric = "RMSE")->
      res

    errorm <- IMR:::error_metric$rel.rmse

    res %<>%
      rbind(sim1_res(dat, fitsi$fit, "SI", errorm) %>%
              rbind(sim1_res(dat, fitimr$best_fit, "IMR", errorm)) %>%
              mutate(metric = "Rel.RMSE")
      )
    res$dim = n
    res$p = p
    res$q = q
    res$miss_pct = mean(dat$mask == 0)
    res$r = r
    res$lambda_beta = hparam$beta$value
    res$lambda_gamma = hparam$gamma$value
    res$lambda_laplace = fitimr$best_fit$lambda_laplace

    all_res %<>% rbind(res)
    print(paste(pct, " in ", round(Sys.time() - start2)))
    print(res)
  }
  print(paste(b, " in ", round(Sys.time() - start1)))
  saveRDS(all_res, "./article/simulation/data/sim2_res.rds")
}

########################################################################
# simulation 2 part 2
results2 <- readRDS("article/simulation/data/sim2_res.rds")

results2 %>%
  filter(model == "IMR", metric == "RMSE") %>%
  dplyr::select(lambda_beta, lambda_gamma, lambda_laplace, miss_pct) %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  filter(miss_pct != .95) %>%
  #group_by(miss_pct) %>%
  summarise_all(mean) %>%
  ungroup() -> best_hparams

results2 %>%
  filter(model == "SI", metric == "RMSE") %>%
  dplyr::select(rank, miss_pct) %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  filter(miss_pct != .95) %>%
  summarise_all(mean) %>%
  round() -> simpute_rank

n = m = 1000
p = 5;
q = 5;
r = 5;
missing_pct = seq(.7, .98, .05)
models <- c("IMR", "SImpute")
all_res <- res <- data.frame()
#simpute_ranks <- c(13, 12, 12, 12, 11, 10)
# b = 1; pct=0.9
for(b in 1:500){
  seed = 2025 + b
  start1 = Sys.time()
  set.seed(seed)
  dat <-
    generate_simulated_data(n, m, r, p, q, .7,
                            sparsity_beta = 0,
                            sparsity_gamma = 0,
                            intercept = FALSE,
                            structured_error_A = F, SNR = 1,
                            structured_error_B = F,
                            prepare_for_fitting = F, mv_coeffs = T, seed = seed)
  for(pct in 1:length(missing_pct)){
    start2 = Sys.time()
    if(pct > 1)
      dat <- increase_sparsity(dat, .05)

    mdat <- IMR::prepare_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2)

    start = Sys.time()
    fitsi <- softImpute::softImpute(dat$Y,
                                    rank.max = simpute_rank$rank,
                                    lambda = best_hparams$lambda_laplace,
                                    thresh = 1e-6,
                                    maxit = 1000,
                                    trace.it = FALSE,final.svd = TRUE, type = "als")
    time.si = as.numeric(Sys.time() -  start, units = "secs")
    start = Sys.time()
    fitimr <- IMR::imr.fit(Y = mdat$Y,
                           X = mdat$Xq,
                           Z = mdat$Zq,
                           intercept_row = FALSE,
                           intercept_col = FALSE,
                           r = r,
                           lambda_M = best_hparams$lambda_laplace,
                           lambda_beta = best_hparams$lambda_beta,
                           lambda_gamma = best_hparams$lambda_gamma,
                           maxit = 1000,
                           thresh = 1e-6,
                           trace=FALSE,
                           shared_information = FALSE,
                           ls_initial = TRUE
    );
    time.imr = as.numeric(Sys.time() -  start, units = "secs")

    dat$Xr <- mdat$Xr
    dat$Zr <- mdat$Zr
    errorm <- IMR:::error_metric$rel.rmse
    mutate(sim1_res(dat, fitsi, "SI", errorm),time=time.si) %>%
      rbind(mutate(sim1_res(dat, fitimr, "IMR", errorm), time=time.imr)) %>%
      mutate(metric = "RMSE")->
      res

    errorm <- IMR:::error_metric$rel.rmse

    res %<>%
      rbind(mutate(sim1_res(dat, fitsi, "SI", errorm),time=time.si) %>%
              rbind(mutate(sim1_res(dat, fitimr, "IMR", errorm), time=time.imr)) %>%
              mutate(metric = "Rel.RMSE")
      )
    res$dim = n
    res$p = p
    res$q = q
    res$miss_pct = mean(dat$mask == 0)
    res$r = r
    res$lambda_beta = best_hparams$lambda_beta
    res$lambda_gamma = best_hparams$lambda_gamma
    res$lambda_laplace = best_hparams$lambda_laplace

    all_res %<>% rbind(res)
    print(paste(pct, " in ", round(Sys.time() - start2)))
    print(res)
  }
  print(paste(b, " in ", round(Sys.time() - start1)))
  saveRDS(all_res, "./article/simulation/data/sim2_3_res.rds")
}



#======================== end ============================================
