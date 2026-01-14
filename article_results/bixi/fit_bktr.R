



#-----------------------------------------------------------------------
#' fit BKTR model to the data
fit_BKTR_to_Bixi <- function(miss      = 0.8,
                             timestamp = format(Sys.Date(), "%d%b"),
                             seed = 2025){
  require(BKTR)
  require(data.table)
  set.seed(seed)

  TSR$set_params(seed = seed, fp_type = "float32", fp_device = 'cpu')

    bixi.dat <- BixiData$new()
  file_prefix <- paste0(
    "./article_results/bixi/data/splits/",
    round(100*miss),
    "percent_",
    timestamp,
    "_"
  )
  train.df <- readRDS(paste0(file_prefix, "train.rds"))




  train.df %<>% rename(location=column, time = row)

  train.df <- setkey(as.data.table(train.df), location, time)

  bixi.dat$temporal_positions_df %<>%
    filter(time %in% train.df$time)
  bixi.dat$spatial_positions_df %<>%
    filter(location %in% train.df$location)

  start_time <- Sys.time()

  p_lgth <- KernelParameter$new(value = 7, is_fixed = TRUE)
  k_local_periodic <- KernelSE$new() * KernelPeriodic$new(period_length = p_lgth)

  bktr_fit <- BKTRRegressor$new(
    formula = nb_departure ~ 1 + mean_temp_c + area_park + total_precip_mm,
    data_df = train.df,
    spatial_positions_df = bixi.dat$spatial_positions_df,
    temporal_positions_df = bixi.dat$temporal_positions_df,
    rank = 8,
    spatial_kernel = KernelMatern$new(smoothness_factor = 5),
    temporal_kernel = k_local_periodic,
    burn_in_iter = 1000,
    sampling_iter = 500)
  bktr_fit$mcmc_sampling()


  time <- round(as.numeric(
    difftime(Sys.time(), start_time, units = "secs")))
  return_obj <- list(fit = bktr_fit, time = time)

  saveRDS(return_obj, paste0("./article_results/bixi/data/bktr_",round(100*miss),"_fit.rds"))
  #bktr_fit$beta_covariates_summary %>% knitr::kable() %>% print()
  #bktr_fit |> summary() %>% print()

  return(return_obj)
}

#-----------------------------------------------------------------------
