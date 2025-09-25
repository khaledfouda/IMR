# required libraries:
# magrittr, tidyverse

preprocess_bixi_data <- function(miss_pct = 0.85,
                                 timestamp = format(Sys.Date(), "%d%b"),
                                 seed = 2025) {
  require(BKTR)
  require(corrr)
  # Set seed for reproducibility ---------------------------
  set.seed(seed)

  # Load raw data ------------------------------------------
  bixi_data  <- BixiData$new()
  data_df    <- bixi_data$data_df

  # Remove rows/locations with all-missing departures --------
  # remove empty rows and columns of response matrix ---------
  data_df %<>%
    group_by(time) %>%
    filter(!all(is.na(nb_departure))) %>%
    ungroup() %>%
    group_by(location) %>%
    filter(!all(is.na(nb_departure))) %>%
    ungroup()


  # Select covariates & reshape for matrix input -----------
  # For now, I will keep all covariates.
  z_vars    <- c("mean_temp_c", "total_precip_mm", "holiday", "max_temp_f", "humidity")


  data_df %<>%
    rename(column = location, row = time, y = nb_departure) %>%
    arrange(row, column) %>%
    rename_with(~ paste0("z_", .x), .cols = setdiff(names(.),
              c(z_vars, "row", "column", "y"))) %>%
    rename_with(~ paste0("x_", .x), .cols=any_of(z_vars))



  # check correlations with the response variable ------
  data_df %>%
    as.data.frame() %>%
    dplyr::select(-row, -column) %>%
    as.matrix() %>%
    replace_na(0) %>%
    corrr::correlate() %>%
    corrr::stretch() %>%
    as.data.frame() %>%
    filter(x == 'y') %>%
    arrange(desc(abs(r))) %>%
    mutate(r = round(r, 2)) %>%
    print()

  #----------------------------------------------------------------------------------

  # Initialize train/test -----------------------------------


  # Determine dimensions & thresholds -----------------------
  num_rows    <- length(unique(data_df$row))
  num_columns <- length(unique(data_df$column))
  min_obs     <- 13
  # but consider that we remove 95% or missing_rate of each column
  #min_obs <- min_obs / ( (1 - miss_pct) * 0.8)
  # the 0.8 is for the CV with 5 folds
  message(paste0("Response data matrix dimension is:", num_rows,
                 "x",num_columns))

  data_df %>%
    group_by(column) %>%
    summarize(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    filter(na.sum < min_obs) %>%
    select(column) -> cols_to_remove
  message("Removing the following columns for containing less observations",
          " than required ", cols_to_remove)
  # Identify columns with enough data ----------------------
  test_df <-
    train_df <-
    data_df %>%
    group_by(column) %>%
    filter(sum(!is.na(y)) >= min_obs) %>%
    ungroup()

  train_df %>%
    group_by(column) %>%
    summarise(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    head(10) %>% print()
  #
  ## Sample half of them for masking [seed] ---------------
  # actually no, use all of eligible columns + those
  set.seed(seed)
  #====================================================================
  # we will leave the columns with very low number of observations
  # we will choose those with min_obs*90
  train <- train_df; test <- test_df
  train_df <-
    train_df |>
    mutate(
      row_id = row_number(),
      orig_y = y
    ) %>%
    group_by(column) %>%
    mutate(num_observed = sum(!is.na(y))) %>%
    ungroup() %>%
    mutate(low_obs = num_observed < (min_obs *(1/(1-miss_pct))) * (miss_pct) ) %>%
    select(-num_observed)

  train_df %>%
    filter(low_obs == TRUE) %>%
    select(column) %>%
    unique() %>% nrow() -> num_guarded_cols

  message("Guarding ", num_guarded_cols, " columns.")

  test_df <- test_df |>
    mutate(
      row_id            = row_number(),
      orig_y = y
    )
  n_total     <- nrow(train_df)
  n_orig_na   <- sum(is.na(train_df$orig_y))
  n_target_na <- floor(miss_pct * n_total)       # want 95% missing total
  n_to_mask   <- n_target_na - n_orig_na

  if (n_to_mask < 0) {
    stop("Already more than 95% missing in train_df; nothing to do.")
  }
  if (n_to_mask > (n_total - n_orig_na)) {
    stop("Not enough non-missing values to reach 95% missing.")
  }

  mask_ids <-
    train_df |>
    filter(!is.na(orig_y)) |>
    filter(!low_obs) |>
    slice_sample(n = n_to_mask) |>
    pull(row_id)
  # train -> set masks_ids to Na. inverse to test
  train_df <- train_df |>
    mutate(
      y =
        if_else(
          row_id %in% mask_ids,
          NA_real_,
          orig_y
        )
    ) |>
    select(-row_id, -orig_y)

  test_df <- test_df |>
    mutate(
      y =
        if_else(
          row_id %in% mask_ids,
          orig_y,
          NA_real_
        )
    ) |>
    select(-row_id, -orig_y)
  #=================================================================
  # recheck the amount the missing per column

  train_df %>%
    group_by(column) %>%
    summarize(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    filter(na.sum < min_obs) %>%
    select(column) -> cols_to_remove
  message("Removing ", length(unique(cols_to_remove$column)),
          " columns for containing less observations",
          " than required")
  train_df %<>%
    filter(!(column %in% cols_to_remove$column)) %>%
    select(-low_obs) %>%
    arrange(row, column)
  test_df %<>%
    filter(!(column %in% cols_to_remove$column)) %>%
    arrange(row, column)
  message("Number of resulting columns: ", length(unique(train_df$column)))

  #-------------------
  train_df %>%
    group_by(column) %>%
    summarise(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    head(10) %>% print()
  # check that they do not overlap
  stopifnot(sum( (!is.na(train_df$y)) &
                   (!(is.na(test_df$y)))) == 0)
  # Add extra 20% random missing in train_df ---------------

  # print
  message(paste0("Percentage of observations in training set:",
                 round(100*mean(!is.na(train_df$y)),1),
                 "% and in the test set: ",
                 round(100*mean(!is.na(test_df$y)),1),
                 "%."))
  message(paste0("Response data matrix dimension is:",
                 length(unique(train_df$row)),
                 "x",
                 length(unique(train_df$column))))

  # Finalize test_df (drop NA departures) -----------------
  test_df <- test_df %>%
    filter(!is.na(y))

    #train_df  %<>% rename(location=column, time=row)
    #test_df   %<>% rename(location=column, time=row)
  # rename the rows and columns back to time and location for consistency

  # Save splits --------------------------------------------
  file_prefix <- paste0(
    "./article_results/bixi/data/splits/",
    round(miss_pct * 100),
    "percent_",
    timestamp,
    "_"
  )
  saveRDS(train_df, file = paste0(file_prefix, "train.rds"))
  saveRDS(test_df,  file = paste0(file_prefix, "test.rds"))
}
######################


#' Load BIXI model data and build masks/splits
#'
#' @param time_cov Logical; if TRUE, use time-varying covariates
#' @return A list with matrices X, Y, masks, and splits.
prepare_bixi_data <- function(miss_p = 0.8,
                              timestamp = "25Sep") {

  # Read pre-saved train/test splits --------------
  file_prefix <- paste0(
    "./article_results/bixi/data/splits/",
    round(100*miss_p),
    "percent_",
    timestamp,
    "_"
  )
  if(file.exists(paste0(file_prefix, "train.rds"))){
    stop("File already exists: ", paste0(file_prefix, "train.rds"))
  }
  train_df <- readRDS(paste0(file_prefix, "train.rds"))
  test_df  <- readRDS(paste0(file_prefix, "test.rds"))


  # Build covariate matrix X ----------------------
  X <- train_df %>%
    select(row, starts_with("x_")) %>%
    group_by(row) %>%
    slice(1) %>%
    ungroup() %>%
    select(-row)

  Z <- train_df %>%
    select(column, starts_with("z_")) %>%
    group_by(column) %>%
    slice(1) %>%
    ungroup() %>%
    select(-column)

  # Build response matrix Y -----------------------
  Y <- reshape2::dcast(
    train_df,
    row ~ column,
    value.var = "y"
  ) %>%
    select(-row) %>%
    as.matrix()
  # Print the missing value percentages
  require(glue)
  Y %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    unlist(use.names = FALSE) ->
    col_na

  Y %>%
    t() %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    unlist(use.names = FALSE) ->
    row_na

  message(glue(
    "Observed‐value counts (pct):\n",
    " Columns   : min = {min(col_na)}({round(100*min(col_na)/dim(Y)[1],1)}%), ",
    " max = {max(col_na)}({round(100*max(col_na)/dim(Y)[1],1)}%)\n",
    " Rows      : min = {min(row_na)}({round(100*min(row_na)/dim(Y)[2],1)}%), ",
    " max = {max(row_na)}({round(100*max(row_na)/dim(Y)[2],1)}%)\n",
    " Train     : {round(100*mean(!is.na(Y)),1)}%"
  ))
  # Observation mask -----------------------------
  obs_mask <- (!is.na(Y)) * 1
  # print(sum(obs_mask == 1) / length(obs_mask))

  # Identify test entries via merge --------------
  mixed <- train_df %>%
    select(row, column, y) %>%
    merge(
      select(test_df, row, column, y),
      by = c("row", "column"),
      all.x = TRUE
    ) %>%
    as.data.frame() %>%
    arrange(row, column) %>%
    mutate(
      # missing is true for observed training observations + missing
      # missing is false for the test set
      missing = !(is.na(y.x) & !is.na(y.y))
    )
  # print(sum(mixed$missing) / nrow(mixed))

  test_mask <- reshape2::dcast(
    mixed,
    row ~ column,
    value.var = "missing"
  ) %>%
    select(-row) %>%
    { colnames(.) <- NULL; . } %>%
    as.matrix() * 1
  # print(sum(1 - test_mask) / length(test_mask))
  # obs mask is 1 for training  and 0 for missing
  # print(obs_mask[1:5, 1:5])
  # test mask is 0 for test and 1 otherwise
  # print(test_mask[1:5, 1:5])

  # Validation mask via MC split -----------------
  # valid mask is 0 for validation and 1 otherwise
  valid_mask <- IMR:::mask_train_test_split(obs_mask, testp = 0.2)

  # Assemble core model.dat -----------------------
  model_dat <- list(
    X     = as.matrix(X),
    Z     = as.matrix(Z),
    Y     = Y,
    masks = list(
      tr_val = obs_mask,
      test   = test_mask,
      valid  = valid_mask
    )
  )

  # Prepare test matrix ---------------------------
  test_mat <- train_df %>%
    select(row, column, y) %>%
    merge(
      select(test_df, row, column, y),
      by = c("row", "column"),
      all.x = TRUE
    ) %>%
    as.data.frame() %>%
    arrange(row, column) %>%
    select(row, column, y.y) %>%
    reshape2::dcast(
      row ~ column,
      value.var = "y.y"
    ) %>%
    select(-row) %>%
    as.matrix() %>%
    IMR::as.Incomplete()

  message(glue("Test      : {round(100*sum(test_mat!=0)/length(test_mat),1)}%"))
  # message(glue("Test*train {sum((test_mat*Y),na.rm=T)}"))
  # Compile splits --------------------------------
  Xqr <- qr(model_dat$X)
  Zqr <- qr(model_dat$Z)
  model_dat$Xq <- qr.Q(Xqr)
  model_dat$Xr <- qr.R(Xqr)
  model_dat$Zq <- qr.Q(Zqr)
  model_dat$Zr <- qr.R(Zqr)
  model_dat$Y.inc <- IMR::as.Incomplete(model_dat$Y)

  model_dat$splits <- list(
    train = IMR::as.Incomplete(model_dat$Y * valid_mask),
    valid = IMR::as.Incomplete(model_dat$Y * (1 - valid_mask)),
    test  = test_mat
  )

  # print(length(model_dat$splits$train@x) / length(model_dat$Y))
  # print(length(model_dat$splits$test@x))
  # print(length(model_dat$splits$valid@x))

  return(model_dat)
}

#--------------------------------------------------------------------------------------------
