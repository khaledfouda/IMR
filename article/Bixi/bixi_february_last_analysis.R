library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(magrittr)
source("./article/Bixi/data_bixi.R")
source("./other_models/SoftImpute_cv.R")
source("./article/bixi/fit_bktr.R")
source("./article/bixi/helpers_bixi.R")
#------------------------------------------------
# we begin by generating the data
# we generate 10 test sets of 30%
# or just keep the original 20% missing as it's convenient
original_missing_pct <- 0.1163287
total_miss <- 0.25
test_pct <- total_miss - original_missing_pct
test_pct

# ================================================================
# we now train >>

model_combn <- expand.grid(
  similarity = c(TRUE, FALSE),
  covariates = F, #c(T, F),
  Intercepts = T, #c(T, F),
  stringsAsFactors = FALSE
)

total_results <- data.frame()
train_seq <- round(seq(1 - total_miss, by = -.05, length.out = 5) * 100)
vars_to_keep <- c(
  "model", "lambda_m", "error.test", "corr.test", "error.train",
  "time0", "time1", "time2"
)
# =================================================================
# we read BKTR
seed <- 4000
results <- data.frame()
# prefix=1
# train_size = train_seq[1]
for (prefix in 1:10) {
  for (train_size in train_seq) {
    dat <- prepare_bixi_data(total_miss, "Feb_last", seed,
      prefix = prefix,
      train_prefix = train_size,
      val_prop = 0.2,
      bktr_variables = TRUE,
      file_dir = "./article/bixi/data/splits2/",
      temporal = model_combn$kernels[1],
      spatial = model_combn$kernels[1],
      temporal_jitter = TRUE,
      spatial_jitter = TRUE,
      jitter_kappa_max = 1e3,
      jitter_tau_max = 1e-2
    )

    bktr_out <- BKTR_Bixi_Wrapper(
      dat = dat,
      miss = total_miss,
      timestamp = "Feb_last",
      prefix = prefix,
      train_prefix = train_size,
      file_dir = "./article/bixi/data/splits2/",
      seed = seed,
      return_fit = TRUE,
      burn_in_iter = 1000,
      sampling_iter = 500,
      test_error = IMR:::error_metric$rmse
    )


    results %<>% rbind(as.data.frame(bktr_out$results[vars_to_keep]) %>%
      mutate(
        metric = "RMSE",
        train_size = train_size,
        rank_estim = bktr_out$fit$fit$rank_decomp
      ))

    bktr_out <- BKTR_Bixi_Wrapper(
      dat = dat,
      miss = total_miss,
      timestamp = "Feb_last",
      prefix = prefix,
      train_prefix = train_size,
      file_dir = "./article/bixi/data/splits2/",
      seed = seed,
      return_fit = TRUE,
      burn_in_iter = 1000,
      sampling_iter = 500,
      test_error = IMR:::error_metric$rel.rmse
    )


    results %<>% rbind(as.data.frame(bktr_out$results[vars_to_keep]) %>%
      mutate(
        metric = "RRMSE",
        train_size = train_size,
        rank_estim = bktr_out$fit$fit$rank_decomp
      ))
  }
  print(paste0("prefix=", prefix))
}
saveRDS(results, "./article/Bixi/data/final_results/BKTR_results_final_25pct_2.rds")
# ====================================================================
results <- readRDS("./article/Bixi/data/final_results/BKTR_results_final_25pct_2.rds")


#=====================================================================
results2 <- readRDS("./article/Bixi/data/final_results/IMR_results_final_25pct_2_2.rds") %>%
  rbind(readRDS("./article/Bixi/data/final_results/IMR_results_final_25pct_2_3.rds")) %>%
  rename(
    time0 = time,
    time1 = time2.1,
    time2 = time2.2
  ) %>%
  group_by(model, train_size, prefix, metric) %>%
  # summarise_all(mean) %>%
  slice_min(test, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(model, train_size, metric, test, time0, time1, time2, rank_estim)

results2 %<>% rbind(
  results %>%
    dplyr::select(model, train_size, metric, error.test, time0, time1, time2, rank_estim) %>%
    rename(test = error.test)
)


results2 %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(train_size, test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time0 = round(time0, 2))

# ==============================================================

readRDS("./article/Bixi/data/final_results/IMR_results_final_25pct_2.rds") %>%
  filter(metric == "RMSE") %>%
  select(test, model, prefix, lambda_laplace, rank_M, train_size) %>%
  group_by(model, train_size, prefix) %>%
  slice_min(test, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  view()

# =========================================================
results3 <- readRDS("./article/Bixi/data/final_results/IMR_results_onefit_25pct_2_3.rds") %>%
  # rbind(readRDS("./article/Bixi/data/final_results/IMR_results_onefit_25pct_2_3.rds")) %>%
  rename(
    time0 = time,
    time1 = time2.1,
    time2 = time2.2
  ) %>%
  dplyr::select(model, train_size, metric, test, time0, time1, time2, rank_estim)

results3 %<>% rbind(
  results %>%
    dplyr::select(model, train_size, metric, error.test, time0, time1, time2, rank_estim) %>%
    rename(test = error.test)
)


results3 %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(train_size, test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time0 = round(time0, 2)) %>%
  select(model, train_size, test, time0, rank_estim) %>% kable()
# =============================================
# article table >
require(gt)
results3 |>
  filter(metric == "RRMSE") %>%
  #mutate(time0 = if_else(model=="BKTR", time0/60, time0))%>%
  group_by(model, train_size, metric) %>%
  summarise_all(c(mean=mean, sd=sd)) %>%
  as.data.frame() %>%
  arrange(train_size, test_mean) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  select(model, train_size, contains("test"), contains("time0")) %>%
  mutate(
    time_mean = as.numeric(time0_mean),
    model = case_match(
      model,
      "similarity" ~ "IMR-Sim",
      "original" ~ "IMR-Identity",
      .default = model
    )
  ) |>
  select(-time0_mean) %>%
  mutate(time = paste0(round(time_mean,2)," (",round(time0_sd,2),")")) %>%
  #mutate(time = if_else(model=="BKTR", paste(time,"min"), paste(time, "s"))) |>
  mutate(test = paste0(round(test_mean,4),"  (",round(test_sd,4),")")) %>%
  select(-time0_sd, -test_sd) %>%
  pivot_wider(
    names_from = model,
    values_from = c(test, time, time_mean, test_mean)
  ) |>
  mutate(
    row_id = row_number(),
    speedup_sim = (time_mean_BKTR*60) / `time_mean_IMR-Sim`,
    speedup_ident = (time_mean_BKTR*60) / `time_mean_IMR-Identity`,
    diff_sim = (`test_mean_IMR-Sim` - test_mean_BKTR) / test_mean_BKTR,
    diff_ident = (`test_mean_IMR-Identity` - test_mean_BKTR) / test_mean_BKTR
  ) |>
  select(
    row_id,
    train_size,
    test_BKTR, time_BKTR,
    test_mean_BKTR, time_mean_BKTR,
    `test_mean_IMR-Sim`, `time_mean_IMR-Sim`,
    `test_mean_IMR-Identity`, `time_mean_IMR-Identity`,
    `test_IMR-Sim`, `time_IMR-Sim`, speedup_sim, #diff_sim,
    `test_IMR-Identity`, `time_IMR-Identity`, speedup_ident#, diff_ident
  ) %>%
  arrange(train_size) %>%
  gt() |>
  tab_spanner(
    "Test RRMSE",
    c(test_BKTR, `test_IMR-Sim`, `test_IMR-Identity`)
  ) %>%
  tab_spanner(
    "Computation Time",
    c(time_BKTR, `time_IMR-Sim`, speedup_sim, `time_IMR-Identity`, speedup_ident)
  ) %>%
  cols_label(
    train_size = "Train Size",
    test_BKTR = "BKTR",
    time_BKTR = "BKTR",
    `test_IMR-Sim` = "IMR-S",
    `time_IMR-Sim` = "IMR-S",
    speedup_sim = "Speedup",
    #diff_sim = "Diff %",
    `test_IMR-Identity` = "IMR-N",
    `time_IMR-Identity` = "IMR-N",
    speedup_ident = "Speedup",
    #diff_ident = "Diff %"
  ) |>
  fmt_number(
    columns = starts_with("speedup"),
    decimals = 0,
    pattern = "×{x}"
  ) |>
  fmt_number(
    columns = starts_with("test"),
    decimals = 4,
  ) |>
  fmt_number(
    columns = contains("train_size"),
    decimals = 0,
    pattern = '{x}\\%'
  ) |>
  tab_style(
    style = cell_fill(color = "gray94"), # Light gray
    locations = cells_body(
      rows = row_id %% 2 == 0
    )
  ) |>
  tab_style(
   style = cell_text(weight = "bold"),
   locations = list(cells_column_labels(),cells_column_spanners(),
                    cells_body(c(train_size, `time_IMR-Sim`,`time_IMR-Identity`) ))
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners(spanners = everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = test_BKTR,
      rows = test_mean_BKTR <= `test_mean_IMR-Sim` & test_mean_BKTR <= `test_mean_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `test_IMR-Sim`,
      rows = `test_mean_IMR-Sim` <= test_mean_BKTR & `test_mean_IMR-Sim` <= `test_mean_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `test_IMR-Identity`,
      rows = `test_mean_IMR-Identity` <= test_mean_BKTR & `test_mean_IMR-Identity` <= `test_mean_IMR-Sim`
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `speedup_sim`,
      rows = `time_mean_IMR-Sim` <= time_mean_BKTR
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `speedup_ident`,
      rows = `time_mean_IMR-Identity` <= time_mean_BKTR
    )) |>
  cols_hide(c(row_id,contains("Speed"))) %>%
      cols_hide(columns = c(contains("_mean"), contains("_sd"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")

table %>%
  tab_header(
    title ="Comparison of Test RRSME"
  ) %>%
  cols_hide(columns = c(contains("time"), contains("speed"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


results3 |>
  filter(metric == "RRMSE") %>%
  select(-metric) |>
  rename(time=time0) |>
  mutate(time = as.numeric(time)) |>
  #mutate(time0 = if_else(model=="BKTR", time0/60, time0))%>%
  group_by(model, train_size) %>%
  summarise_all(c(mean=mean, sd=sd)) %>%
  as.data.frame() %>%
  arrange(train_size, test_mean) %>%
  group_by(train_size) |>
  mutate(
    ref_time = time_mean[model=="BKTR"],
    speedup = if_else(model=="BKTR", NA, ref_time/time_mean),
    is_best_rmse = test_mean == min(test_mean),
    is_fast = !is.na(speedup) & speedup > 1
  ) |>
  ungroup() |>
  mutate(
    model = case_match(
      model,
      "simulated" ~ "IMR-Similarity",
      "none" ~ "IMR-Identity",
      .default = model
    )
  ) |>
  mutate(
    rmse_str = sprintf("%.4f (%.4f)", test_mean, test_sd),
    time_str = if_else(
      model == "BKTR",
      sprintf("%.2f (%.2f) min", time_mean/60, time_sd/60),
      sprintf("%.2f (%.2f) s", time_mean, time_sd)
    )
  ) |>
  arrange(train_size, factor(model, levels = c("BKTR", "IMR-Similarity", "IMR-Identity"))) |>
  select(train_size, model, rmse_str, time_str, speedup, is_best_rmse, is_fast) |>
  gt(groupname_col = "train_size") |>
  cols_label(
    model = "Model",
    rmse_str = "Test RMSE",
    time_str = "Time",
    speedup = "Speedup"
  ) |>
  fmt_number(
    columns = speedup,
    decimals = 0,
    pattern = "×{x}"
    #na_str = "—"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = rmse_str,
      rows = is_best_rmse
    )
  ) |>

  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = speedup,
      rows = is_fast
    )
  ) |>
  cols_hide(columns = c(is_best_rmse, is_fast, speedup)) |>
  cols_align(align = "center", columns = everything()) |>
  cols_align(align = "left", columns = model)





  gt() |>
  # tab_spanner(
  #   label = "BKTR",
  #   columns = c(test_BKTR, time_BKTR)
  # ) |>
  # tab_spanner(
  #   label = "IMR-Similarity",
  #   columns = c(`test_IMR-Sim`, `time_IMR-Sim`, speedup_sim, diff_sim)
  # ) |>
  # tab_spanner(
  #   label = "IMR-Identity",
  #   columns = c(`test_IMR-Identity`, `time_IMR-Identity`, speedup_ident, diff_ident)
  # ) |>
  tab_spanner(
    "Test RRMSE",
    c(test_BKTR, `test_IMR-Sim`, `test_IMR-Identity`)
  ) %>%
  tab_spanner(
    "Computation Time",
    c(time_BKTR, `time_IMR-Sim`, speedup_sim, `time_IMR-Identity`, speedup_ident)
  ) %>%
  cols_label(
    train_size = "Train Size",
    test_BKTR = "BKTR",
    time_BKTR = "BKTR (min)",
    `test_IMR-Sim` = "IMR-Similarity",
    `time_IMR-Sim` = "IMR-Similarity (s)",
    speedup_sim = "Speedup",
    #diff_sim = "Diff %",
    `test_IMR-Identity` = "IMR-Identity",
    `time_IMR-Identity` = "IMR-Identity (s)",
    speedup_ident = "Speedup",
    #diff_ident = "Diff %"
  ) |>
  fmt_number(
    columns = starts_with("speedup"),
    decimals = 0,
    pattern = "×{x}"
  ) |>
  fmt_number(
    columns = contains("train_size"),
    decimals = 0,
    pattern = '{x}\\%'
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(cells_column_labels(),cells_column_spanners(),
                     cells_body(train_size))
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners(spanners = everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = test_BKTR,
      rows = test_mean_BKTR <= `test_mean_IMR-Sim` & test_mean_BKTR <= `test_mean_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `test_IMR-Sim`,
      rows = `test_mean_IMR-Sim` <= test_mean_BKTR & `test_mean_IMR-Sim` <= `test_mean_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `test_IMR-Identity`,
      rows = `test_mean_IMR-Identity` <= test_mean_BKTR & `test_mean_IMR-Identity` <= `test_mean_IMR-Sim`
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `speedup_sim`,
      rows = `time_mean_IMR-Sim` <= time_mean_BKTR
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `speedup_ident`,
      rows = `time_mean_IMR-Identity` <= time_mean_BKTR
    )) |>
  cols_hide(contains("Speed")) %>%
  cols_hide(columns = c(contains("_mean"), contains("_sd"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")

table %>%
  tab_header(
    title ="Comparison of Test RRSME"
  ) %>%
  cols_hide(columns = c(contains("time"), contains("speed"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")
