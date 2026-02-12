library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(magrittr)
source("./article/bixi/data_bixi.R")
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
  kernels = c("simulated", "none"),
  covariates = c(T, F),
  Intercepts = c(T, F),
  stringsAsFactors = FALSE
)
model_combn <- model_combn[c(7, 8), ]

total_results <- data.frame()
train_seq <- round(seq(1 - total_miss, by = -.05, length.out = 5) * 100)
vars_to_keep <- c(
  "model", "lambda_M", "error.test", "corr.test", "error.train",
  "time0", "time1", "time2"
)
# =================================================================
# we read BKTR
seed <- 4000
results <- data.frame()
# prefix=1
# train_size = train_seq[1]
for (prefix in 1:6) {
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
saveRDS(results, "./article/Bixi/data/final_results/BKTR_results_final_25pct.rds")
# ====================================================================

results2 <- readRDS("./article/Bixi/data/final_results/IMR_results_final_25pct_2.rds") %>%
  rename(
    time0 = time,
    time1 = time2.1,
    time2 = time2.2
  ) %>%
  group_by(model, train_size, prefix, metric) %>%
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
results3 <- readRDS("./article/Bixi/data/final_results/IMR_results_onefit_25pct.rds") %>%
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
  select(model, train_size, test, time0)
# =============================================
# article table >

results3 |>
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(train_size, test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time0 = time0/60)%>%#round(time0, 2)) %>%
  select(model, train_size, test, time0) %>%
  mutate(
    time = as.numeric(time0),
    model = case_match(
      model,
      "simulated" ~ "IMR-Sim",
      "none" ~ "IMR-Identity",
      .default = model
    )
  ) |>
  select(-time0) %>%

  pivot_wider(
    names_from = model,
    values_from = c(test, time)
  ) |>
  mutate(
    speedup_sim = time_BKTR / `time_IMR-Sim`,
    speedup_ident = time_BKTR / `time_IMR-Identity`,
    diff_sim = (`test_IMR-Sim` - test_BKTR) / test_BKTR,
    diff_ident = (`test_IMR-Identity` - test_BKTR) / test_BKTR
  ) |>
  select(
    train_size,
    test_BKTR, time_BKTR,
    `test_IMR-Sim`, `time_IMR-Sim`, speedup_sim, diff_sim,
    `test_IMR-Identity`, `time_IMR-Identity`, speedup_ident, diff_ident
  ) %>%
  mutate(`time_IMR-Sim` = round(`time_IMR-Sim`*60,2),
         `time_IMR-Identity` = round(`time_IMR-Identity`*60,2)
         ) %>%
  arrange(train_size) %>%
  gt() |>
  tab_spanner(
    label = "BKTR",
    columns = c(test_BKTR, time_BKTR)
  ) |>
  tab_spanner(
    label = "IMR-Similarity",
    columns = c(`test_IMR-Sim`, `time_IMR-Sim`, speedup_sim, diff_sim)
  ) |>
  tab_spanner(
    label = "IMR-Identity",
    columns = c(`test_IMR-Identity`, `time_IMR-Identity`, speedup_ident, diff_ident)
  ) |>
  fmt_number(
    columns = contains("test"),
    decimals = 4
  ) |>
  fmt_number(
    columns = contains("time"),
    decimals = 2,
    pattern = "{x}"
  ) |>
  fmt_number(
    columns = contains("train_size"),
    decimals = 0,
    pattern = '{x}%'
  ) |>
  fmt_number(
    columns = starts_with("speedup"),
    decimals = 0,
    pattern = "×{x}"
  ) |>
  fmt_percent(
    columns = starts_with("diff"),
    decimals = 2,
    force_sign = TRUE
  ) |>
  cols_label(
    train_size = "Train Size",
    test_BKTR = "RRSME",
    time_BKTR = "Time (min)",
    `test_IMR-Sim` = "RRSME",
    `time_IMR-Sim` = "Time (s)",
    speedup_sim = "Speedup",
    diff_sim = "Diff %",
    `test_IMR-Identity` = "RRSME",
    `time_IMR-Identity` = "Time (s)",
    speedup_ident = "Speedup",
    diff_ident = "Diff %"
  ) |>
  #tab_style(
  #  style = cell_text(weight = "bold"),
  #  locations = cells_column_labels()
  #) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners(spanners = everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = test_BKTR,
      rows = test_BKTR <= `test_IMR-Sim` & test_BKTR <= `test_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = `test_IMR-Sim`,
      rows = `test_IMR-Sim` <= test_BKTR & `test_IMR-Sim` <= `test_IMR-Identity`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = `test_IMR-Identity`,
      rows = `test_IMR-Identity` <= test_BKTR & `test_IMR-Identity` <= `test_IMR-Sim`
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = `speedup_sim`,
      rows = `time_IMR-Sim` <= time_BKTR
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", color = "#2f4f4f"),
    locations = cells_body(
      columns = `speedup_ident`,
      rows = `time_IMR-Identity` <= time_BKTR
    )
  )-> table
table %>%
  tab_header(
    title = "Comparsion of Training Time"
    # subtitle = "Comparison of Test RRSME and Training Time vs BKTR Baseline"
  ) %>%
  cols_hide(columns = c(contains("time"), contains("speed"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")

table %>%
  tab_header(
    title = "Comparison of Test RRSME"
    # subtitle = "Comparison of Test RRSME and Training Time vs BKTR Baseline"
  ) |>
  cols_hide(columns = c(contains("test"), contains("Diff"))) %T>%
  gtsave("article/Bixi/data/final_results/table1.tex")
