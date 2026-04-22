# required libraries:
# magrittr, tidyverse

mutate_bixi_file <- function(data,
                           suffix,
                           miss_pct,
                           timestamp,
                           prefix = "",
                           file_override = FALSE,
                           create_folder = FALSE,
                           out_dir = "./article_results/bixi/data/splits/",
                           read = FALSE) {

  if (!dir.exists(out_dir)) {
    if(create_folder){
      dir.create(out_dir, recursive = TRUE)
    }else stop("Folder does not exist")
  }

  filename <- paste0(
    round(miss_pct * 100), "percent_",
    timestamp, "_"
  )

  if (prefix != "") {
    filename <- paste0(filename, prefix, "_")
  }

  filename <- paste0(filename, suffix, ".rds")
  full_path <- file.path(out_dir, filename)

  if (file.exists(full_path)) {
    if(read) return(readRDS(full_path))
    if (!file_override) {
      stop("File already exists: ", full_path)
    }
    message("Warning: Overwriting existing file: ", full_path)
  }
  if(read){
    stop("Attempting to read nonexisting file.", full_path)
  }else
    saveRDS(data, file = full_path)
  #invisible(full_path)
}



preprocess_bixi_data <- function(miss_pct = 0.5,
                                 timestamp = format(Sys.Date(), "%d%b"),
                                 seed = 2025,
                                 prefix = "",
                                 file_override = FALSE,
                                 decreasing_train = FALSE,
                                 create_folder = FALSE,
                                 out_dir = "./article_results/bixi/data/splits/",
                                 train_n_steps = 5,
                                 train_stepsize = 0.05) {
  require(BKTR)
  require(corrr)
  # Set seed for reproducibility ---------------------------
  set.seed(seed)

  # Load raw data ------------------------------------------
  bixi_data <- BixiData$new()
  data_df <- bixi_data$data_df

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
  z_vars <- c("mean_temp_c", "total_precip_mm", "holiday", "max_temp_f", "humidity")


  data_df %<>%
    rename(column = location, row = time, y = nb_departure) %>%
    arrange(row, column) %>%
    rename_with(~ paste0("z_", .x), .cols = setdiff(
      names(.),
      c(z_vars, "row", "column", "y")
    )) %>%
    rename_with(~ paste0("x_", .x), .cols = any_of(z_vars))


  # check correlations with the response variable ------
  data_df %>%
    as.data.frame() %>%
    dplyr::select(-row, -column) %>%
    as.matrix() %>%
    replace_na(0) %>%
    corrr::correlate() %>%
    corrr::stretch() %>%
    as.data.frame() %>%
    filter(x == "y") %>%
    arrange(desc(abs(r))) %>%
    mutate(r = round(r, 2)) %>%
    print()
  #-------------------------------------------------------------------------
  # the following 8 stations have less than 5 observations so I'll discard them.
  low.obs.stations <- c(
    "6194 - Métro Atwater (Atwater / Ste-Catherine)",
    "6019 - Métro Sherbrooke (de Rigaud / Berri)",
    "6036 - de la Commune / St-Sulpice",
    "6181 - Clark / Rachel",
    "6157 - de Brébeuf / du Mont-Royal",
    "6227 - de l'Esplanade / Laurier",
    "6136 - Métro Laurier (Rivard / Laurier)",
    "6184 - Métro Mont-Royal (Rivard / du Mont-Royal)"
  )
  data_df %<>%
    filter(!(column %in% low.obs.stations))
  #----------------------------------------------------------------------------------
  # this is the part where the test set was fixed no matter the missing percentage
  # I'll remove for now since we no longer need it
  do.not.run <- FALSE
  if (do.not.run) {
    test50 <- readRDS("article_results/bixi/data/splits/50percent_25Sep_test.rds")
    data_df %<>%
      left_join(
        test50 %>%
          select(row, column) %>%
          mutate(flag = TRUE),
        by = c("row", "column")
      ) %>%
      mutate(u = ifelse(isTRUE(flag), NA, y)) %>%
      select(-flag)
  }

  # Initialize train/test -----------------------------------


  # Determine dimensions & thresholds -----------------------
  num_rows <- length(unique(data_df$row))
  num_columns <- length(unique(data_df$column))
  # but consider that we remove 95% or missing_rate of each column
  # min_obs <- min_obs / ( (1 - miss_pct) * 0.8)
  # the 0.8 is for the CV with 5 folds
  message(paste0(
    "Response data matrix dimension is:", num_rows,
    "x", num_columns
  ))

  min_obs <- 0
  if (min_obs > 0) {
    data_df %>%
      group_by(column) %>%
      summarize(na.sum = sum(!is.na(y))) %>%
      arrange(na.sum) %>%
      filter(na.sum < min_obs) %>%
      select(column) -> cols_to_remove
    message(
      "Removing the following columns for containing less observations",
      " than required ", cols_to_remove
    )
  }
  # Identify columns with enough data ----------------------

  data_df %>%
    group_by(column) %>%
    summarise(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    head(10) %>%
    print()
  #
  ## Sample half of them for masking [seed] ---------------
  # actually no, use all of eligible columns + those
  set.seed(seed)
  # ====================================================================
  # we will leave the columns with very low number of observations
  # we will choose those with min_obs*90
  # train <- train_df; test <- test_df
  train_df <-
    data_df |>
    mutate(
      row_id = row_number(),
      orig_y = y
    )

  test_df <- data_df |>
    mutate(
      row_id = row_number(),
      orig_y = y
    )
  n_total <- nrow(train_df)
  n_orig_na <- sum(is.na(train_df$orig_y))
  n_target_na <- floor(miss_pct * n_total) # want 95% missing total
  n_to_mask <- n_target_na - n_orig_na

  if (n_to_mask < 0) {
    stop("Already more than 95% missing in train_df; nothing to do.")
  }
  if (n_to_mask > (n_total - n_orig_na)) {
    stop("Not enough non-missing values to reach 95% missing.")
  }

  # mask_ids are the ids to set 0 and be a test set.
  mask_ids <-
    train_df |>
    filter(!is.na(orig_y)) |>
    # filter(!low_obs) |>
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
    select(-orig_y)

  test_df <- test_df |>
    mutate(
      y =
        if_else(
          row_id %in% mask_ids,
          orig_y,
          NA_real_
        )
    ) |>
    select(-orig_y)



  if (decreasing_train) {
    train_seq <- seq(1-miss_pct, by = -train_stepsize,
                     length.out = train_n_steps)
    train_seq <- round(train_seq*100)
    all_trains <- list()
    all_trains[[1]] <- train_df
    current_train <- train_df
    new_prefix = paste0(prefix,"_train", train_seq[1])
    mutate_bixi_file(current_train, "train", miss_pct, timestamp, new_prefix, file_override,
                   create_folder, out_dir)
    for (i in 2:length(train_seq))
    {
      train_seq[2]
      n_total <- nrow(current_train)
      n_orig_na <- sum(is.na(current_train$y))
      n_to_mask <- floor(0.05 * n_total)
      #mean(is.na(new_train$y))


      mask_ids <-
        current_train |>
        filter(!is.na(y)) |>
        slice_sample(n = n_to_mask) |>
        pull(row_id)
      # train -> set masks_ids to Na. inverse to test
      all_trains[[i]] <- current_train |>
        mutate(
          y =
            if_else(
              row_id %in% mask_ids,
              NA_real_,
              y
            )
        )
      current_train <- all_trains[[i]]
      new_prefix = paste0(prefix,"_train", train_seq[i])
      mutate_bixi_file(current_train, "train", miss_pct, timestamp, new_prefix, file_override,
                     create_folder, out_dir)
    }
  }
  # -----------------------------------------------------------------


  # =================================================================
  # recheck the amount the missing per column
  if (min_obs > 0) {
    train_df %>%
      group_by(column) %>%
      summarize(na.sum = sum(!is.na(y))) %>%
      arrange(na.sum) %>%
      filter(na.sum < min_obs) %>%
      select(column) -> cols_to_remove
    print(cols_to_remove$column)
    message(
      "Removing ", length(unique(cols_to_remove$column)),
      " columns for containing less observations",
      " than required"
    )
    train_df %<>%
      filter(!(column %in% cols_to_remove$column)) %>%
      select(-low_obs) %>%
      arrange(row, column)
    test_df %<>%
      filter(!(column %in% cols_to_remove$column)) %>%
      arrange(row, column)
    message("Number of resulting columns: ", length(unique(train_df$column)))
  }

  #-------------------
  train_df %>%
    group_by(column) %>%
    summarise(na.sum = sum(!is.na(y))) %>%
    arrange(na.sum) %>%
    head(10) %>%
    print()
  # check that they do not overlap
  stopifnot(sum((!is.na(train_df$y)) &
    (!(is.na(test_df$y)))) == 0)
  # Add extra 20% random missing in train_df ---------------

  # print
  message(paste0(
    "Percentage of observations in training set:",
    round(100 * mean(!is.na(train_df$y)), 1),
    "% and in the test set: ",
    round(100 * mean(!is.na(test_df$y)), 1),
    "%."
  ))
  message(paste0(
    "Response data matrix dimension is:",
    length(unique(train_df$row)),
    "x",
    length(unique(train_df$column))
  ))

  # Save splits --------------------------------------------
  mutate_bixi_file(train_df, "train", miss_pct, timestamp, prefix, file_override,
                 create_folder, out_dir)
  mutate_bixi_file(test_df, "test", miss_pct, timestamp, prefix, file_override,
                 create_folder, out_dir)

}
######################


#' Load BIXI model data and build masks/splits
#'
#' @param time_cov Logical; if TRUE, use time-varying covariates
#' @return A list with matrices X, Y, masks, and splits.
prepare_bixi_data <- function(miss_p = 0.5,
                              timestamp = "25Sep", seed = NULL,
                              val_prop = 0.2,
                              prefix = "",
                              train_prefix = "",
                              x_keep = c(
                                "x_humidity", "x_max_temp_f", "x_mean_temp_c",
                                "x_total_precip_mm", "x_holiday"
                              ),
                              bktr_variables = FALSE,
                              file_dir = "./article_results/bixi/data/splits/",
                              ...
                              # these parameters are sent to the kernel generation
) {
  # Read pre-saved train/test splits --------------
  if(train_prefix != "")
    new_prefix <- paste0(prefix, "_train", train_prefix)
  train_df <- mutate_bixi_file(NULL, "train", miss_p, timestamp, new_prefix,
                               out_dir = file_dir, read=TRUE)
  test_df <- mutate_bixi_file(NULL, "test", miss_p, timestamp, prefix,
                               out_dir = file_dir, read=TRUE)

  # Build covariate matrix X ----------------------
  X <- train_df %>%
    dplyr::select(row, starts_with("x_")) %>%
    group_by(row) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-row) %>%
    dplyr::select(all_of(x_keep))

  Z <- train_df %>%
    dplyr::select(column, starts_with("z_")) %>%
    group_by(column) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-column)

  if (bktr_variables) {
    X <- X %>%
      dplyr::select(x_mean_temp_c, x_total_precip_mm)
    Z <- Z %>%
      dplyr::select(z_area_park)
  }

  # Build response matrix Y -----------------------
  Y <- reshape2::dcast(
    train_df,
    row ~ column,
    value.var = "y"
  ) %>%
    dplyr::select(-row) %>%
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
  #-------
  # test dataset
  test_set <- reshape2::dcast(
    test_df,
    row ~ column,
    value.var = "y"
  ) %>%
    dplyr::select(-row) %>%
    as.matrix() %>%
    IMR::as_incomplete()


  message(glue("Test      : {round(100*sum(test_set!=0)/length(test_set),1)}%"))


  kernels <- generate_similarity_bixi(miss_p, timestamp, prefix = prefix,
                                      train_prefix = train_prefix,
                                      file_dir =file_dir,...)
  output <- list()
  output$modd <- IMR::prepare_data(Y,
    X = as.matrix(X), Z = as.matrix(Z),
    similarity_rows = kernels$temporal,
    similarity_cols = kernels$spatial,
    seed = seed, val_prop = val_prop
  )

  output$test_mask <- IMR::as_incomplete((test_set != 0) * 1)
  output$test <- test_set
  output$X <- X
  output$Z <- Z
  output$col_names <- colnames(Y)
  output$row_names <- train_df$row %>% unique()



  return(output)
}

#--------------------------------------------------------------------------------------------
