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
    select(location, time, y_est, y) ->
    test.estimates

  # obtain train estimates
  bktr_fit$fit$imputed_y_estimates |>
    as.data.frame() |>
    merge(filter(train, ! is.na(y)), by = c("location", "time")) |>
    select(location, time, y_est, y) ->
    train.estimates

  results <- list(
    model = "BKTR",
    lambda_M = NA,
    lambda_beta = NA
  )

  results <- c(results,
               prepare_output_bixi(NA, dat$X, test.estimates$y_est,
                                   train.estimates$y_est, test.estimates$y,
                                   train.estimates$y))

  results$total_num_fits = 1000
  results$time = bktr_fit$time
  results$time_per_fit = results$time  / results$total_num_fits
  results$cov_summaries = bktr_fit$fit$beta_covariates_summary
  results$sparsity = 0
  if(return_fit) return(list(fit=bktr_fit, results=results))
  results
}



prepare_output_bixi <- function(
    time,
    X,
    estim.test,
    estim.train,
    obs.test,
    obs.train,
    time_per_fit = NA,
    total_num_fits = NA,
    beta.estim  = NA,
    M.estim     = NA,
    test_error  = IMR::error_metric$rmse
) {
  # Core metrics
  results <- list(
    time           = time,
    time_per_fit   = time_per_fit,
    total_num_fits = total_num_fits,
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
    sparsity    = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    )
  )


  # Covariate coefficient summaries
  results$cov_summaries <- tryCatch({
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

  results
}


