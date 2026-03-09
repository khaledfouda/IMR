



#-----------------------------------------------------------------------
#' fit BKTR model to the data
fit_BKTR_to_Bixi <- function(miss      = 0.8,
                             timestamp = format(Sys.Date(), "%d%b"),
                             prefix    = "",
                             train_prefix = "",
                             file_dir = "./article/bixi/data/splits/",
                             seed = 2025,
                             burn_in_iter = 1000,
                             sampling_iter = 500){
  require(BKTR)
  require(data.table)
  set.seed(seed)

  TSR$set_params(seed = seed, fp_type = "float32", fp_device = 'cpu')

    bixi.dat <- BixiData$new()

  if(train_prefix != "")
    train_prefix <- paste0(prefix, "_train", train_prefix)
  train.df <- mutate_bixi_file(NULL, "train", miss, timestamp, train_prefix,
                               out_dir = file_dir, read=TRUE)

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
    formula = y ~ 1 + x_mean_temp_c + z_area_park + x_total_precip_mm,
    data_df = train.df,
    spatial_positions_df = bixi.dat$spatial_positions_df,
    temporal_positions_df = bixi.dat$temporal_positions_df,
    rank = 8,
    spatial_kernel = KernelMatern$new(smoothness_factor = 5),
    temporal_kernel = k_local_periodic,
    burn_in_iter = burn_in_iter,
    sampling_iter = sampling_iter)

  ttime <- bench::bench_time(bktr_fit$mcmc_sampling())


  time <- round(as.numeric(
    difftime(Sys.time(), start_time, units = "secs")))
  return_obj <- list(fit = bktr_fit, time = time,
                     time1 = lubridate::time_length(ttime[1], "seconds"),
                     time2 = lubridate::time_length(ttime[2], "seconds"))


  file_prefix <- paste0(
    "./article/bixi/data/bktr_fits/bktrfit_",
    round(100*miss),
    "percent_",
    timestamp,
    "_"
  )
  if(train_prefix != "")
    prefix <- paste0(prefix, "_train", train_prefix)
  if(prefix != "")
    file_prefix <- paste0(file_prefix, prefix, "_")

  saveRDS(return_obj, paste0(file_prefix, ".rds"))
  #bktr_fit$beta_covariates_summary %>% knitr::kable() %>% print()
  #bktr_fit |> summary() %>% print()
  return(return_obj)
}

#-----------------------------------------------------------------------

