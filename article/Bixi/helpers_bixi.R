BKTR_Bixi_Wrapper <- function(
    dat,
    timestamp = format(Sys.Date(), "%d%b"),
    seed = NULL,
    miss = .8,
    return_fit = FALSE,
    prefix = "",
    train_prefix = "",
    file_dir = "./article/bixi/data/splits/",
    test_error = IMR:::error_metric$rmse,
    ...
){
  if(!is.null(seed)) set.seed(seed)
  file_prefix <- paste0(
    "./article_results/bixi/data/bktrfit_",
    round(100*miss),
    "percent_",
    timestamp,
    "_"
  )

    file_prefix <- paste0(
      "./article/Bixi/data/bktr_fits/bktrfit_",
      round(100*miss),
      "percent_",
      timestamp,
      "_",
      prefix,
      "_train",
      prefix,
      "_train",
      train_prefix,
      "_.rds")

  bktr_fit <- readRDS(file_prefix)


  if(train_prefix != "")
    train_prefix <- paste0(prefix, "_train", train_prefix)
  train <- mutate_bixi_file(NULL, "train", miss, timestamp, train_prefix,
                               out_dir = file_dir, read=TRUE)
  test <- mutate_bixi_file(NULL, "test", miss, timestamp, prefix,
                            out_dir = file_dir, read=TRUE)


    train  %<>% rename(location=column, time=row)
    test   %<>% rename(location=column, time=row)

  # obtain test estimates
  bktr_fit$fit$imputed_y_estimates |>
    as.data.frame() |>
    merge(test, by = c("location", "time")) |>
    dplyr::select(location, time, y_est, y) ->
    test.estimates

  # obtain train estimates
  bktr_fit$fit$imputed_y_estimates |>
    as.data.frame() |>
    merge(filter(train, ! is.na(y)), by = c("location", "time")) |>
    dplyr::select(location, time, y_est, y) ->
    train.estimates

  results <- list(
    model = "BKTR",
    lambda_m = NA,
    lambda_beta = NA
  )

  results <- c(results,
               prepare_output_bixi(X = dat$X,
                                   Z = dat$Z,
                                   estim.test = test.estimates$y_est,
                                   estim.train = train.estimates$y_est,
                                   obs.test = test.estimates$y,
                                   obs.train = train.estimates$y,
                                   test_error = test_error))

  # results$total_num_fits = 1000
  results$time0 = bktr_fit$time
  results$time1 = bktr_fit$time1
  results$time2 = bktr_fit$time2
  # results$time_per_fit = results$time  / results$total_num_fits
  results$cov_summaries = bktr_fit$fit$beta_covariates_summary
  results$sparsity = 0
  if(return_fit) return(list(fit=bktr_fit, results=results))
  results
}



generate_similarity_bixi <- function(miss      = 0.8,
                                     timestamp = "25Sep",
                                     prefix    = "",
                                     train_prefix = "",
                                     file_dir = "./article/bixi/data/splits/",
                                     spatial = "simulated",
                                     temporal = "simulated",
                                     spatial_jitter = TRUE,
                                     temporal_jitter = TRUE,
                                     matern_range = function(x){x <- x[upper.tri(x)]; median(x[x>0])},
                                     matern.cor.target = 0.5,
                                     jitter_kappa_max = 1e3,
                                     jitter_tau_max = 1e-2,
                                     RBF_ell_t = 1.3,
                                     RBF_ell_s = 1.3,
                                     matern_scale = NULL,
                                     return_distance){
  library(BKTR)
  bkdat <- BixiData$new()
  stopifnot(temporal %in% c("Matern", "original", "RBF", "none", "simulated"))
  stopifnot(spatial %in% c("Matern", "original", "none", "RBF", "simulated"))

  if(train_prefix != "")
    train_prefix <- paste0(prefix, "_train", train_prefix)
  train.df <- mutate_bixi_file(NULL, "train", miss, timestamp, train_prefix,
                               out_dir = file_dir, read=TRUE)


  train.df %<>% rename(location=column, time = row)

  bkdat$temporal_positions_df %<>%
    filter(time %in% train.df$time)


  p_lgth <- KernelParameter$new(value = 7, is_fixed = TRUE)

  if(temporal == "simulated"){
    se_lgth <- KernelParameter$new(value = 6.427, is_fixed = TRUE)
    per_lgth <- KernelParameter$new(value = 1.039, is_fixed = TRUE)
  }else{
    se_lgth <- KernelParameter$new(value = 6.448, is_fixed = TRUE)
    per_lgth <- KernelParameter$new(value = 0.941, is_fixed = TRUE)
  }
  temporal_kernel <- KernelSE$new(lengthscale = se_lgth) *
    KernelPeriodic$new(lengthscale = per_lgth, period_length = p_lgth)
  #temporal_kernel = BKTR::KernelSE$new()
  temporal_kernel$set_positions(bkdat$temporal_positions_df)

  bkdat$spatial_positions_df %<>%
    filter(location %in% train.df$location)
  if(spatial == "simulated"){
    sp_lgth <- KernelParameter$new(value = 0.018, is_fixed = TRUE)
  }else
    sp_lgth <- KernelParameter$new(value = 21.128, is_fixed = TRUE)

  spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 5,lengthscale = sp_lgth)
  spatial_kernel$set_positions(bkdat$spatial_positions_df)

  distance = list()
  if(return_distance)
    distance = list(
      spatial = as.matrix(spatial_kernel$distance_matrix),
      temporal = as.matrix(temporal_kernel$distance_matrix)
    )

  choose_jitter <- function(K, jitter_kappa_max = 1e4, jitter_tau_max = 1e-2, tau0 = 1e-12) {
    s <- mean(diag(K))
    tau2 <- tau0 * s
    n <- nrow(K)
    repeat {
      Kt <- K + diag(tau2, n)
        kap <- kappa(Kt)
        #print(kap)
        if (is.finite(kap) && kap <= jitter_kappa_max) return(tau2)

      tau2 <- tau2 * 10
      if (tau2 > jitter_tau_max * s) return(tau2)
    }
  }

  if(temporal == "Matern"){
    d = temporal_kernel$distance_matrix %>% as.matrix()
    range0 <- fields::Matern.cor.to.range(matern_range(d), nu=3/2)
    temporal_kernel <- fields::Matern(d, smoothness = 3/2, range = range0)

  }else if (temporal == "none"){
    temporal_kernel = NULL
    #n = train.df$time %>% unique() %>% length()
    #temporal_kernel = diag(1, n, n)
  }else if (temporal == "RBF"){
    d = temporal_kernel$distance_matrix %>% as.matrix()
    temporal_kernel <- exp(-(d^2) / (2 * RBF_ell_t^2))

  }else if (temporal %in% c("original", "simulated")) {

    temporal_kernel <- temporal_kernel$kernel_gen() %>% as.matrix()
  }
  if(temporal_jitter & temporal != "none") temporal_kernel <- temporal_kernel +
    diag(choose_jitter(temporal_kernel, jitter_kappa_max, jitter_tau_max), nrow(temporal_kernel))
  if(temporal != "none")
    temporal_kernel <- IMR:::inv(temporal_kernel)

  if(spatial == "Matern"){

    d = spatial_kernel$distance_matrix %>% as.matrix()
    if(is.null(matern_scale))
      matern_scale <- fields::Matern.cor.to.range(matern_range(d), nu=5/2, matern.cor.target = matern.cor.target)
    spatial_kernel <- fields::Matern(d, smoothness = 5/2, range = matern_scale)

  }else if (spatial == "none"){
    #n = train.df$location %>% unique() %>% length()
    #spatial_kernel = diag(1, n, n)
    spatial_kernel = NULL
    }else if (spatial %in% c("original", "simulated")){
    spatial_kernel = spatial_kernel$kernel_gen() %>% as.matrix()
  }else if (spatial == "RBF"){
    d = spatial_kernel$distance_matrix %>% as.matrix()
    spatial_kernel <- exp(-(d^2) / (2 * RBF_ell_s^2))
  }

  if(spatial_jitter & spatial != "none") spatial_kernel <- spatial_kernel +
    diag(choose_jitter(spatial_kernel, jitter_kappa_max, jitter_tau_max), nrow(spatial_kernel))
  if(spatial != "none")
    spatial_kernel <- IMR:::inv(spatial_kernel)

  list(spatial=spatial_kernel, temporal=temporal_kernel, distance=distance)
}


prepare_output_bixi <- function(
    X,
    Z,
    estim.test,
    estim.train,
    obs.test,
    obs.train,
    beta.estim  = NA,
    gamma.estim = NA,
    M.estim     = NA,
    test_error  = IMR:::error_metrics$rmse,
    digits = 5
) {
  # Core metrics
  results <- list(
    error.test  = test_error(estim.test, obs.test),
    corr.test   = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    rank_M      = tryCatch(
      qr(M.estim)$rank,
      error = function(e) NA
    ),
    rank_beta   = tryCatch(
      qr(beta.estim)$rank,
      error = function(e) NA
    ),
    sparsity_beta    = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    ),
    rank_gamma   = tryCatch(
      qr(gamma.estim)$rank,
      error = function(e) NA
    ),
    sparsity_gamma    = tryCatch(
      sum(gamma.estim == 0) / length(gamma.estim),
      error = function(e) NA
    )
  )

  # Covariate coefficient summaries
  results$beta_summaries <- tryCatch({
    apply(beta.estim, 1, summary) |>
      as.data.frame() |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        prop_non_zero = apply(beta.estim, 1, function(x)
          sum(x != 0) / length(x)
        )
      ) |>
      `rownames<-`(colnames(X))
  }, error = function(e) NA)

  results$gamma_summaries <- tryCatch({
    apply(gamma.estim, 2, summary) |>
      as.data.frame() |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        prop_non_zero = apply(gamma.estim, 2, function(x)
          sum(x != 0) / length(x)
        )
      ) |>
      `rownames<-`(colnames(Z))
  }, error = function(e) NA)

  results
}


output_wrapper_bixi <- function(fit, dat, shared_information = FALSE,
                                test_error  = IMR:::error_metrics$rmse){



  out <- IMR:::reconstruct(fit, dat$modd,trace = TRUE)

  return(list(
    rec = out,
    res = prepare_output_bixi(
      X           = dat$X,
      Z           = dat$Z,
      estim.test  = out$estimates[as.matrix(dat$test_mask==1)],
      estim.train = out$estimates[as.matrix(dat$modd$obs_mask == 1)],
      obs.test    = dat$test@x,
      obs.train   = dat$modd$Y@x,
      beta.estim  = out$beta,
      gamma.estim = out$gamma,
      M.estim     = out$M,
      test_error  = test_error
    )
  ))
}
