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
                     error_metric=IMR:::error_metrics$rrmse,
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
  out$test_rel <- IMR:::error_metrics$rrmse(estimates[test.obs], dat$theta[test.obs])
  out$test <- error_metric(estimates[test.obs], dat$theta[test.obs])
  train.obs <- dat$Y != 0
  out$train <- error_metric(estimates[train.obs], dat$theta[train.obs])
  out$rank <- qr(estimates)$rank
  return(out)
}

#===========================================================================================
# setting 1)
#===========================================================================================
dims <- c(400, 600, 800, 1000)
p = 10;
q = 0;
r = 4;
missing_pct = 0.8
sparsity_beta = 0.5
models <- c("IMR", "SImpute", "MCCI")
all_res <- res <- data.frame()
convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6, trace=FALSE, ls_initial = TRUE)
grid <- IMR::imr_tune_grid(rank = c(2, 10, 1, 2), beta = 0, laplace = c(0,40,40,2))

for(b in 1:500){
  seed = 2025 + b
  set.seed(seed)
  for(d in dims){

    n = m = d
    dat <-
      generate_simulated_data(n, m, r, p, 0, .8,
                              sparsity_beta = 0,
                              intercept = FALSE,
                              structured_error_A = F, SNR = 1,
                              structured_error_B = F,
                              prepare_for_fitting = F, mv_coeffs = T, seed = seed)


    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, seed = seed, val_prop = 0.2);

    fitsi <- simpute.cv(y_full = mdat$Y,
                        y_train = mdat$y_train,
                        y_valid = mdat$y_valid,
                        trace = FALSE,
                        print.best = FALSE,
                        tol = grid$laplace$streaks,
                        n.lambda = grid$laplace$length,
                        maxit = convergence$maxit,
                        thresh = convergence$thresh,
                        test_error = IMR:::error_metrics$rmse,
                        seed = seed)




    fitimr <- IMR::imr_tune(mdat, grid, convergence=convergence, fast_laplace = FALSE,
                            seed = seed, n_cores = 7, verbose = 1)


    fitmcci <- MCCI.cv(Y = dat$Y, X = dat$X, W = dat$mask, n_folds = 5,numCores = 9,
                       seed = seed,
                       test_error = IMR:::error_metrics$rmse,
                       lambda_1_grid = c(0),#seq(0, 1, length = 10),
                       lambda_2_grid = seq(2.9, 0.1, length = 18),
                       alpha_grid = c(1),#seq(0.992, 1, length = 10),
                       n1n2_optimized = TRUE,
                       return_diagn = FALSE)

    dat$Xr <- mdat$Xr
    errorm <- IMR:::error_metrics$rmse
    sim1_res(dat, fitsi$fit, "SI", errorm) %>% rbind(
      sim1_res(dat, fitmcci$fit, "MCCI", errorm, coeffs_transformed = FALSE),
      sim1_res(dat, fitimr$fit$coefficients, "IMR", errorm)
    ) -> res
    errorm <- IMR:::error_metrics$rrmse
    res %<>% rbind(sim1_res(dat, fitsi$fit, "SI", errorm) %>% rbind(
      sim1_res(dat, fitmcci$fit, "MCCI", errorm, coeffs_transformed = FALSE),
      sim1_res(dat, fitimr$fit$coefficients, "IMR", errorm)
    ))
    res$dim = d
    all_res %<>% rbind(res)
    print(d)
  }
  print(b)
  print(res)
  saveRDS(all_res, "./article/simulation/data/sim1_res.rds")
}




#===================================================
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


convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6, trace=FALSE, ls_initial = TRUE)
grid <- IMR::imr_tune_grid(rank = c(2, 10, 1, 2), beta=c(0), gamma=c(0), laplace=c(0,120,60,2));

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

    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2);

    fitsi <- simpute.cv(y_full = mdat$Y,
                        y_train = mdat$y_train,
                        y_valid = mdat$y_valid,
                        trace = FALSE,
                        print.best = FALSE,
                        tol = grid$laplace$streaks,
                        maxit = convergence$maxit,
                        thresh = convergence$thresh,
                        n.lambda = grid$laplace$length,
                        test_error = IMR:::error_metrics$rmse,
                        seed = seed)





    fitimr <- IMR::imr_tune(mdat, grid, convergence=convergence, seed = seed, fast_laplace = FALSE,
                            n_cores = 7, verbose = 1)

    #print(fitimr$fit)
    #summary(fitimr$fit)



    dat$Xr <- mdat$Xr
    dat$Zr <- mdat$Zr
    errorm <- IMR:::error_metrics$rmse
    sim1_res(dat, fitsi$fit, "SI", errorm) %>%
      rbind(sim1_res(dat, fitimr$fit$coefficients, "IMR", errorm)) %>%
      mutate(metric = "RMSE")->
      res

    errorm <- IMR:::error_metrics$rrmse

    res %<>%
      rbind(sim1_res(dat, fitsi$fit, "SI", errorm) %>%
              rbind(sim1_res(dat, fitimr$fit$coefficients, "IMR", errorm)) %>%
              mutate(metric = "Rel.RMSE")
      )
    res$dim = n
    res$p = p
    res$q = q
    res$miss_pct = mean(dat$mask == 0)
    res$r = r
    res$lambda_beta = fitimr$fit$meta$lambdas["beta"]
    res$lambda_gamma = fitimr$fit$meta$lambdas["gamma"]
    res$lambda_laplace = fitimr$fit$meta$lambdas["M"]
    res$rank_m = fitimr$fit$meta$rank

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
  dplyr::select(lambda_beta, lambda_gamma, lambda_laplace, rank_m, miss_pct) %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  filter(miss_pct != .95) %>%
  #group_by(miss_pct) %>%
  summarise_all(mean) %>%
  mutate(rank_m = round(rank_m)) %>%
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
convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6, trace=FALSE, ls_initial = TRUE)

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

    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2)

    start = Sys.time()
    fitsi <- softImpute::softImpute(dat$Y,
                                    rank.max = simpute_rank$rank,
                                    lambda = best_hparams$lambda_laplace,
                                    thresh = convergence$thresh,
                                    maxit = convergence$maxit,
                                    trace.it = FALSE,final.svd = TRUE, type = "als")
    time.si = as.numeric(Sys.time() -  start, units = "secs")
    start = Sys.time()

    fitimr <- IMR::imr_fit(mdat, rank = best_hparams$rank_m,
                           lambda_m = best_hparams$lambda_laplace,
                           lambda_beta = best_hparams$lambda_beta,
                           lambda_gamma = best_hparams$lambda_gamma,
                           convergence=convergence)
    time.imr = as.numeric(Sys.time() -  start, units = "secs")

    dat$Xr <- mdat$Xr
    dat$Zr <- mdat$Zr
    errorm <- IMR:::error_metrics$rrmse
    mutate(sim1_res(dat, fitsi, "SI", errorm),time=time.si) %>%
      rbind(mutate(sim1_res(dat, fitimr$coefficients, "IMR", errorm), time=time.imr)) %>%
      mutate(metric = "RMSE")->
      res

    errorm <- IMR:::error_metrics$rrmse

    res %<>%
      rbind(mutate(sim1_res(dat, fitsi, "SI", errorm),time=time.si) %>%
              rbind(mutate(sim1_res(dat, fitimr$coefficients, "IMR", errorm), time=time.imr)) %>%
              mutate(metric = "Rel.RMSE")
      )
    res$dim = n
    res$p = p
    res$q = q
    res$miss_pct = mean(dat$mask == 0)
    res$r = r
    res$rank_m = best_hparams$rank_m
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
