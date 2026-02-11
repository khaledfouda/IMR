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
test_pct <- 0.2
total_miss <- test_pct + original_missing_pct -> miss_p
total_miss <- 0.35
test_pct <- total_miss - original_missing_pct


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
vars_to_keep = c('model', 'lambda_M', 'error.test', 'corr.test', 'error.train',
                 'time0', 'time1', 'time2')
# =================================================================
# we read BKTR
seed <- 4000
results <- data.frame()
prefix=1
train_size = train_seq[1]
for (prefix in 1:7) {
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

    bktr_out <- BKTR_Bixi_Wrapper(dat=dat,
                                  miss = total_miss,
                                 timestamp = "Feb_last",
                                 prefix = prefix,
                                 train_prefix = train_size,
                                 file_dir = "./article/bixi/data/splits2/",
                                 seed = seed,
                                 return_fit = TRUE,
                                 burn_in_iter = 1000,
                                 sampling_iter = 500,
                                 test_error = IMR:::error_metric$rmse)


  results %<>% rbind(as.data.frame(bktr_out$results[vars_to_keep]) %>%
                       mutate(metric="RMSE",
                                                       train_size=train_size,
                                                       rank_estim = bktr_out$fit$fit$rank_decomp))

  bktr_out <- BKTR_Bixi_Wrapper(dat=dat,
                                miss = total_miss,
                                timestamp = "Feb_last",
                                prefix = prefix,
                                train_prefix = train_size,
                                file_dir = "./article/bixi/data/splits2/",
                                seed = seed,
                                return_fit = TRUE,
                                burn_in_iter = 1000,
                                sampling_iter = 500,
                                test_error = IMR:::error_metric$rel.rmse)


  results %<>% rbind(as.data.frame(bktr_out$results[vars_to_keep]) %>%
                       mutate(metric="RRMSE",
                                                               train_size=train_size,
                                                               rank_estim = bktr_out$fit$fit$rank_decomp))

  }
  print(paste0("prefix=",prefix))
}

bktr_out %>% names()

#====================================================================

results2 <- readRDS("./article/Bixi/data/final_results/IMR_results_final.rds") %>%
  rename(time0 = time,
         time1 = time2.1,
         time2 = time2.2) %>%
  dplyr::select(model, train_size, metric, test, time0, time1, time2, rank_estim)

results2 %<>% rbind(
results %>%
  dplyr::select(model, train_size, metric, error.test, time0, time1, time2, rank_estim) %>%
  rename(test = error.test) )


results2 %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(train_size, test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  mutate(time0 = round(time0, 2))

#==============================================================
