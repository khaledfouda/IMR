library(dplyr)
library(tidyr)
library(IMR)
library(tidyverse)
library(magrittr)
library(data.table)
library(lubridate)
library(BKTR)

generate_bixi_sample <- function(split_id = 1, train_size = 65) {

  miss_pct <- 0.25
  train_stepsize <- 0.05
  seed <- 2025 + split_id

  set.seed(seed)
  bixi_data <- BKTR::BixiData$new()
  data_df <- bixi_data$data_df

  data_df <- data_df |>
    dplyr::group_by(time) |>
    dplyr::filter(!all(is.na(nb_departure))) |>
    dplyr::ungroup() |>
    dplyr::group_by(location) |>
    dplyr::filter(!all(is.na(nb_departure))) |>
    dplyr::ungroup()

  z_vars <- c("mean_temp_c", "total_precip_mm", "holiday", "max_temp_f", "humidity")

  data_df <- data_df |>
    dplyr::rename(column = location, row = time, y = nb_departure) |>
    dplyr::arrange(row, column)

  current_cols <- names(data_df)
  cols_to_z <- setdiff(current_cols, c(z_vars, "row", "column", "y"))

  data_df <- data_df |>
    dplyr::rename_with(\(x) paste0("z_", x), dplyr::all_of(cols_to_z)) |>
    dplyr::rename_with(\(x) paste0("x_", x), dplyr::any_of(z_vars))

  low_obs_stations <- c(
    "6194 - Métro Atwater (Atwater / Ste-Catherine)",
    "6019 - Métro Sherbrooke (de Rigaud / Berri)",
    "6036 - de la Commune / St-Sulpice",
    "6181 - Clark / Rachel",
    "6157 - de Brébeuf / du Mont-Royal",
    "6227 - de l'Esplanade / Laurier",
    "6136 - Métro Laurier (Rivard / Laurier)",
    "6184 - Métro Mont-Royal (Rivard / du Mont-Royal)"
  )

  data_df <- data_df |>
    dplyr::filter(!(column %in% low_obs_stations))

  train_df <- data_df |> dplyr::mutate(row_id = dplyr::row_number(), orig_y = y)
  test_df <- data_df |> dplyr::mutate(row_id = dplyr::row_number(), orig_y = y)

  n_total <- nrow(train_df)
  n_orig_na <- sum(is.na(train_df$orig_y))
  n_target_na <- floor(miss_pct * n_total)
  n_to_mask <- n_target_na - n_orig_na

  mask_ids <- train_df |>
    dplyr::filter(!is.na(orig_y)) |>
    dplyr::slice_sample(n = n_to_mask) |>
    dplyr::pull(row_id)

  train_df <- train_df |>
    dplyr::mutate(y = dplyr::if_else(row_id %in% mask_ids, NA_real_, orig_y)) |>
    dplyr::select(-orig_y)

  test_df <- test_df |>
    dplyr::mutate(y = dplyr::if_else(row_id %in% mask_ids, orig_y, NA_real_)) |>
    dplyr::select(-orig_y)

  train_seq <- round(seq(1 - miss_pct, by = -train_stepsize, length.out = 5) * 100)
  target_idx <- which(train_seq == train_size)

  if(target_idx > 1) {
    for (i in 2:target_idx) {
      step_mask_count <- floor(train_stepsize * nrow(train_df))

      step_mask_ids <- train_df |>
        dplyr::filter(!is.na(y)) |>
        dplyr::slice_sample(n = step_mask_count) |>
        dplyr::pull(row_id)

      train_df <- train_df |>
        dplyr::mutate(y = dplyr::if_else(row_id %in% step_mask_ids, NA_real_, y))
    }
  }

  train_df <- train_df |> dplyr::mutate(row = as.Date(row))
  test_df <- test_df |> dplyr::mutate(row = as.Date(row))

  X <- train_df |>
    dplyr::select(row, dplyr::starts_with("x_")) |>
    dplyr::distinct(row, .keep_all = TRUE) |>
    dplyr::select(-row) |>
    dplyr::select(x_mean_temp_c, x_total_precip_mm)

  Z <- train_df |>
    dplyr::select(column, dplyr::starts_with("z_")) |>
    dplyr::distinct(column, .keep_all = TRUE) |>
    dplyr::select(-column) |>
    dplyr::select(z_area_park)

  Y <- train_df |>
    dplyr::select(row, column, y) |>
    tidyr::pivot_wider(names_from = column, values_from = y) |>
    dplyr::select(-row) |>
    as.matrix()

  test_set <- test_df |>
    dplyr::select(row, column, y) |>
    tidyr::pivot_wider(names_from = column, values_from = y) |>
    dplyr::select(-row) |>
    as.matrix() |>
    IMR::as_incomplete()

  train_df_pos <- train_df |> dplyr::rename(location = column, time = row)

  bkdat <- BKTR::BixiData$new()
  bkdat$temporal_positions_df %<>%
    filter(time %in% train_df_pos$time)

  p_lgth <- BKTR::KernelParameter$new(value = 7, is_fixed = TRUE)
  se_lgth <- BKTR::KernelParameter$new(value = 6.427, is_fixed = TRUE)
  per_lgth <- BKTR::KernelParameter$new(value = 1.039, is_fixed = TRUE)
  temporal_kernel <- BKTR::KernelSE$new(lengthscale = se_lgth) *
    BKTR::KernelPeriodic$new(lengthscale = per_lgth, period_length = p_lgth)
  temporal_kernel$set_positions(bkdat$temporal_positions_df)

  bkdat$spatial_positions_df %<>%
    filter(location %in% train_df_pos$location)
  sp_lgth <- BKTR::KernelParameter$new(value = 0.018, is_fixed = TRUE)
  spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 5,lengthscale = sp_lgth)
  spatial_kernel$set_positions(bkdat$spatial_positions_df)

  list(
    Y = Y,
    X = X,
    Z = Z,
    test = test_set,
    spatial_distance = as.matrix(spatial_kernel$distance_matrix),
    temporal_positions = as.matrix(temporal_kernel$distance_matrix)
  )
}

n <- 100; m <- 150

bixi_example <- generate_bixi_sample()

rind <- sample.int(nrow(bixi_example$Y), n)
cind <- sample.int(ncol(bixi_example$Y), m)

Bixi_sample <- list( Y = bixi_example$Y[rind,cind],
                     test = bixi_example$test[rind,cind],
                     X = as.matrix(bixi_example$X[rind,]),
                     Z = as.matrix(bixi_example$Z[cind,]),
                     spatial_distance = bixi_example$spatial_distance[cind,cind],
                     temporal_distance = bixi_example$temporal_positions[rind,rind])

usethis::use_data(Bixi_sample, overwrite = TRUE, compress = "xz")
