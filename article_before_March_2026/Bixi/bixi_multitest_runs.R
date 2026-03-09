library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(dplyr)
library(magrittr)
source("./article_results/bixi/data_bixi.R")
source("./other_models/SoftImpute_cv.R")
source("./article_results/bixi/fit_bktr.R")
source("./article_results/bixi/helpers_bixi.R")
#-------------------------------------------------------------------------------
#' we will create multiple test set with the same size but with different seeds
#' First: we will run BKTR on all of them and save them to a file.
#' Start with 50 different test sets.
#------------------------------------------------------------
# step 1: generate the data
num_testsets = 20
timestamp = "Jan"
original_missing_pct = 0.1163287
test_pct <- 0.2
total_miss = test_pct + original_missing_pct
seed = 2025
for(prefix in 1:num_testsets){
  seed = seed + 1
  preprocess_bixi_data(total_miss, timestamp, seed, prefix,T)
}
#----------------------------------------------------------------
# step 2: train BKTR on all of them

for(prefix in 1:8){
  o = fit_BKTR_to_Bixi(total_miss, timestamp, prefix, seed = 0)

  dat <- prepare_bixi_data(total_miss, timestamp, seed = seed, prefix = prefix,
                           val_prop = vals, temporal = "original", spatial="none",
                           temporal_jitter = T, spatial_jitter = T,
                           bktr_variables = TRUE,
                           kappa_max = 1e3, tau_max = 1e-2)
  bktr_res <- BKTR_Bixi_Wrapper(dat, timestamp, seed, total_miss, T, prefix)
  print(bktr_res$results$error.test)
  print(bktr_res$results$time)
}
#-------------------------------------------------------------
# step 3: train our model + soft-impute >>
bktr_res <- data.frame()
for(prefix in c(1:8,14,15)){

  dat <- prepare_bixi_data(total_miss, timestamp, seed = seed, prefix = prefix,
                           val_prop = vals, temporal = "original", spatial="none",
                           temporal_jitter = T, spatial_jitter = T,
                           bktr_variables = TRUE,
                           kappa_max = 1e3, tau_max = 1e-2)
  res <- BKTR_Bixi_Wrapper(dat, timestamp, seed, total_miss, T, prefix)
  bktr_res %<>% rbind(data.frame(temporal = NA,
                                 spatial=NA,
                                 model    = "BKTR",
                                 val_size = vals,
                                 intercept = NA,
                                 covariate = NA,
                                 shared    = NA,
                                 lambda_m = NA,
                                 rank     = res$fit$fit$rank_decomp,
                                 prefix  = prefix,
                                 seed     = NA, #corrtarg = cc,
                                 error= res$results$error.test,
                                 counter = NA,
                                 time1 = res$results$time[2],
                                 time2 = res$results$time[3],
                                 mape = NA
                                 )
                      );
}

si_res <- readRDS("./article_results/bixi/data/results_jan26/SIfit_20testsplits.rds")
IMR_res <- readRDS("./article_results/bixi/data/results_jan26/imrfit_20testsplits.rds") %>%
  mutate(model = "IMR (spatial=Identity)")
IMRsim_res <- readRDS("./article_results/bixi/data/results_jan26/imrfit_20testsplits_simulated2.rds") %>%
  mutate(model = "IMR (spatial=Matern)")
#--------
bktr_res %>%
  rbind(si_res %>% mutate(model = "SoftImpute")) %>%
  rbind(IMR_res) %>% #-> results_df
  rbind(IMRsim_res) ->results_df
  dplyr::select(-temporal, -spatial, -intercept, -covariate, -shared, -seed, -counter, -mape) %>%
  group_by(model) %>%
  summarize_all(mean) %>%
  arrange(error)

require(gt)


summary_table <- results_df |>
  summarise(
    # Error Statistics
    Avg_RMSE = mean(error),
    SD_RMSE  = sd(error),

    # Complexity Statistics
    Avg_Rank = mean(rank),
    SD_Rank = sd(rank),

    # Efficiency Statistics (using time2 = elapsed time)
    Avg_Time_Sec = mean(time2),

    .by = model
  ) |>
  arrange(Avg_RMSE) # Sort by best accuracy

# Create a beautiful table for the report
summary_table |>
  gt() |>
  tab_header(
    title = "Model Performance",
    subtitle = "Averaged over 10 Train-Test Splits and 10 random replications for each"
  ) |>
  fmt_number(
    columns = c(Avg_RMSE, SD_RMSE),
    decimals = 5
  ) |>
  fmt_number(
    columns = c(Avg_Rank, SD_Rank, Avg_Time_Sec),
    decimals = 1
  ) |>
  cols_label(
    Avg_RMSE = "RMSE (AVG)",
    SD_RMSE = "RMSE (SD)",
    Avg_Rank = "Rank (AVG)",
    SD_Rank = "Rank (SD)",
    Avg_Time_Sec = "Time (s)"
  ) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(columns = model)
  )



efficiency_table <- summary_table |>
  mutate(
    # How much better is the RMSE compared to SoftImpute?
    # (Assuming SoftImpute is the 'worst' error, used as baseline)
    RMSE_Improvement_Pct = (1 - (Avg_RMSE / max(Avg_RMSE))) * 100,

    # How many times slower is the model compared to SoftImpute?
    # (Assuming SoftImpute is the fastest)
    Time_Factor = Avg_Time_Sec / min(Avg_Time_Sec)
  ) |>
  arrange(Avg_RMSE) %>%
  select(model, -Avg_RMSE, RMSE_Improvement_Pct, -Avg_Time_Sec, Time_Factor)

# Formatting the Efficiency Table
efficiency_table |>
  gt() |>
  tab_header(
    title = "Performance Relative to SoftImpute"#,
    # subtitle = "Relative to SoftImpute"
  ) |>
  # fmt_number(columns = c(Avg_RMSE), decimals = 4) |>
  fmt_number(columns = c(RMSE_Improvement_Pct), decimals = 2, pattern = "+{x}%") |>
  fmt_number(columns = c(Time_Factor), decimals = 1) |>
  cols_label(
    RMSE_Improvement_Pct = "RMSE reduction by",
    Time_Factor = "Slower (x) times"
  ) |>
  tab_style(
    style = cell_text(color = "red", weight = "bold"),
    locations = cells_body(
      columns = Time_Factor,
      rows = Time_Factor > 100 # Highlight extremely slow models
    )
  )


#================================================

dat <- prepare_bixi_data(total_miss, timestamp, seed, val_prop = 0.2, prefix="test",
                  spatial = "original", temporal="original",
                  spatial_jitter = TRUE,
                  temporal_jitter = TRUE,
                  kappa_max = 1e3,
                  tau_max =1e-1)




o = prepare_bixi_data(test_pct, "25Sep")

as.matrix(o$modd$obs_mask) * 1 -> m1
as.matrix(o$test_mask) * 1 -> m2

# m1 = 0 when test or missing
# m2 = 0 when train or missing

sum(m1 == 0) / length(m1)
sum(m2 == 0) / length(m2)

(-sum(m2 == 1) + sum(m1 == 0)) / length(m1)




mean( m2 * m1) # m2 = 1 and m1 = 1 # both test and observed? this is wrong!!
mean(())


sum(o$test_mask@x)/length(o$test_mask) #-
sum(o$modd$obs_mask@x)/length(o$test_mask)

is.na(train_df$y) -> mask1
length(test_df$y)

s0
