BKTR_Bixi_Wrapper <- function(
    dat,
    timestamp = format(Sys.Date(), "%d%b"),
    seed = NULL,
    miss = .8,
    return_fit = FALSE,
    ...
){
  if(!is.null(seed)) set.seed(seed)

  bktr_fit <- readRDS(paste0('./article_results/bixi/data/bktr_',
                             paste0(round(miss*100)),
                             '_fit.rds'))
  file_prefix <- paste0(
    "./article_results/bixi/data/splits/",
    round(100*miss),
    "percent_",
    timestamp,
    "_"
  )
  train <- readRDS(paste0(file_prefix, "train.rds"))
  test  <- readRDS(paste0(file_prefix, "test.rds"))


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
    lambda_M = NA,
    lambda_beta = NA
  )

  results <- c(results,
               prepare_output_bixi(X = dat$X,
                                   Z = dat$Z,
                                   estim.test = test.estimates$y_est,
                                   estim.train = train.estimates$y_est,
                                   obs.test = test.estimates$y,
                                   obs.train = train.estimates$y))

  # results$total_num_fits = 1000
  results$time = bktr_fit$time
  # results$time_per_fit = results$time  / results$total_num_fits
  results$cov_summaries = bktr_fit$fit$beta_covariates_summary
  results$sparsity = 0
  if(return_fit) return(list(fit=bktr_fit, results=results))
  results
}


generate_similarity_bixi <- function(miss      = 0.8,
                                     timestamp = format(Sys.Date(), "%d%b")){
  require(BKTR)
  bkdat <- BixiData$new()

  file_prefix <- paste0(
    "./article_results/bixi/data/splits/",
    round(100*miss),
    "percent_",
    timestamp,
    "_"
  )
  train.df <- readRDS(paste0(file_prefix, "train.rds"))

  train.df %<>% rename(location=column, time = row)
  #train.df <- setkey(as.data.table(train.df), location, time)

  bkdat$temporal_positions_df %<>%
    filter(time %in% train.df$time)
  bkdat$spatial_positions_df %<>%
    filter(location %in% train.df$location)

  spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 3)
  temporal_kernel = BKTR::KernelSE$new()
  spatial_kernel$set_positions(bkdat$spatial_positions_df)
  temporal_kernel$set_positions(bkdat$temporal_positions_df)

  symmetrize <- function(x) (x + t(x))/2

  spatial_kernel$kernel_gen()  %>%
    as.matrix() -> #%>%
    solve() %>%
    symmetrize() ->
    spatial_kernel
  temporal_kernel$kernel_gen() %>% as.matrix() %>%
    solve() %>% symmetrize() -> temporal_kernel
  list(spatial=spatial_kernel, temporal=temporal_kernel)
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
    test_error  = IMR:::error_metric$rmse,
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


output_wrapper_bixi <- function(fit, dat, shared_information = FALSE){


  out <- IMR:::reconstruct(fit, dat$modd,shared_information =  shared_information)

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
      M.estim     = out$M
    )
  ))
}
