library(devtools)

# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(magrittr)
source("./article_results/bixi/data_bixi.R")
source("./other_models/SoftImpute_cv.R")
source("./article_results/bixi/fit_bktr.R")
source("./article_results/bixi/helpers_bixi.R")
#------------------------------------------------

miss = 0.5
seed = 2026

dat <- prepare_bixi_data(miss, "25Sep", seed = seed,
                         val_prop = 0.1,
                         x_keep = c("x_max_temp_f",
                                    "x_total_precip_mm","x_holiday" ))
#--------------------------------------------------------
#-- change similarity matrices


generate_similarity_bixi(0.55, "25Sep", temporal = "none", spatial = "none",
                         temporal_jitter=T, spatial_jitter = FALSE,
                         ell=1.3,
                         #matern_range = function(x) .001,
                         # matern_range =  function(D, q = 0.5) {
                         #   v <- D[upper.tri(D)]
                         #   v <- v[is.finite(v) & v > 0]
                         #   unname(quantile(v, probs = q, type = 7))
                         # },
                         kappa_max = 1e3) -> ks
ks$spatial[1:10,1:10]
dat$modd$similarity_cols <- ks$spatial
dat$modd$similarity_rows <- ks$temporal


ks$temporal[1:10,1:10]
#--------------


hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                       dat$modd$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(1, 0.1, 0.01)
hparam$laplace$min <- 0
hparam$laplace$max <- 5
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
hparam$rank$max <- 15
hparam$beta$n_streaks <- hparam$gamma$n_streaks <- 2
hparam$beta$max <- 0.31
hparam$gamma$max <- 0.2
hparam$beta$step_sizes <- 1e-7
hparam$gamma$step_sizes <- 0.2 / 300


initialize_parallel_workers(9)

bench::bench_time(res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                             trace=2, hpar=hparam, intercept_row = T,
                                             intercept_col = T, ls_initial = T,
                                             seed = seed, warm_start = NULL, maxit=600,
                                             shared_information = T, thresh=1e-5,
                                             num_cores = 0)) -> timee
round(lubridate::time_length(timee[2], "seconds"), 2)
s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
s0$res$error.test

vestim <- IMR::reconstruct_partial(res0$best_fit, dat$modd, dat$modd$y_valid, T)
error_metric$rmse(dat$modd$y_valid@x, vestim@x)



rec <- IMR::reconstruct(res$best_fit, dat$modd,shared_information = T)

rec$beta
rec$gamma
res$best_fit$beta0

#=====================================================================
# Q: choosing the kernel (no covariates needed here; use val_prop = 0.2)
{
B = c(5, 10, 20, 25)
miss = 0.5
temporal_kernels = c("original")
spatial_kernels = c("original")
combined <- expand.grid(
  temporal = temporal_kernels,
  spatial  = spatial_kernels,
  stringsAsFactors = FALSE
)
all_res <- data.frame()
for(i in B){
initialize_parallel_workers(9)
  seed = 2025 #+ i
  for(j in 1:nrow(combined)){
    temporal = combined$temporal[j]
    spatial = combined$spatial[j]
    dat <- prepare_bixi_data(miss, "25Sep", seed = seed,
                             val_prop = 0.05)
    dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL
    kernels <- generate_similarity_bixi(miss, "25Sep", temporal = temporal, spatial=spatial,
                                        temporal_jitter = F, spatial_jitter = F, ell_t = 1.3,
                                        kappa_max = 1e3, cor.target = 0.5, matern_scale  = i)
    dat$modd$similarity_cols <- kernels$spatial
    dat$modd$similarity_rows <- kernels$temporal
    if(spatial == "none") dat$modd$similarity_cols = NULL
    if(temporal == "none") dat$modd$similarity_rows = NULL
    hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                           dat$modd$similarity_cols, 0, 0)
    hparam$laplace$step_sizes <- c(0.1, 0.01)
    hparam$laplace$min <- 0
    hparam$laplace$max <- 1
    hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
    hparam$rank$max <- 15
    res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                 trace=1, hpar=hparam, intercept_row = T,
                                 intercept_col = T, ls_initial = T,
                                 seed = seed, warm_start = NULL, maxit=600,
                                 shared_information = T, thresh=1e-6,
                                 num_cores = 0)
    s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
    all_res %<>% rbind (data.frame(temporal = temporal, spatial=spatial, error=s0$res$error.test,i=i))
  print(all_res)
  }
  print(i)
}

all_res %>% arrange(error) %>% mutate(error = round(error,6))
}

saveRDS(all_res,"./article_results/bixi/data/results_jan26/kernel_results.rds")

all_res %>%

  group_by(temporal, spatial) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  # as.data.frame() %>%
  # mutate(error = round(error,5)) %>%
  arrange((error))

library(dplyr)
library(tidyr)
library(knitr)

all_res %>%
  group_by(temporal, spatial) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  filter(spatial != "original", temporal != "original") %>%
  mutate(error = round(error, 4)) %>%
  as.data.frame() %>%
  pivot_wider(names_from = spatial, values_from = error) %>%
  # arrange(error) %>%
  kable( caption = "Mean error by temporal and spatial kernel (original excluded).")

library(dplyr)
library(DT)

tab <- all_res %>%
  group_by(temporal, spatial) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  filter(spatial != "original", temporal != "original") %>%
  mutate(error = round(error, 4)) %>%
  arrange(error)

min_err <- min(tab$error, na.rm = TRUE)

datatable(tab, rownames = FALSE, options = list(pageLength = nrow(tab), dom = "t")) %>%
  formatStyle(
    "error",
    backgroundColor = styleEqual(min_err, "khaki"),
    fontWeight = styleEqual(min_err, "bold")
  )
#=============================================================
#' date: January 13
#' Purpose: verifying the effect of the random seed on the train-test validation
#' No covariaties, with intercept, original similarity, no jitter.
#' 1) variable  seed x val_size.  show distribution of error.
#' 2) Fix training size; change val_size and see test error. replicate with different seeds
#' 3) Test 1) with soft-impute.
#' -----------------------------------------------------------
# we begin by generating the data. important here:
# you do not generate a new test with each seed. you generate a new validation with each seed
timestamp = "Jan"
original_missing_pct = 0.1163287
test_pct <- 0.3
total_miss = test_pct + original_missing_pct -> miss_p
seed = 2025
#preprocess_bixi_data(total_miss, timestamp, seed, "",F)

{
  #B = c(5, 10, 20, 25)
  with_intercept = c(F)
  miss = total_miss
  temporal = "none"
  spatial = "none"
  n_seeds = 100
  val_size = c(.05, .1, .2, .3)


  all_res <- data.frame()
  for(i in 1:n_seeds){
    IMR::initialize_parallel_workers(9)
    seed = 2025 + i
    for(vals in val_size){
      dat <- prepare_bixi_data(miss, timestamp, seed = seed,
                               val_prop = vals, temporal = temporal, spatial=spatial,
                               temporal_jitter = T, spatial_jitter = T,
                               kappa_max = 1e4, tau_max = 1e-2)

      dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL

      #kernels <- generate_similarity_bixi(miss, "25Sep", temporal = temporal, spatial=spatial,
       #                                   temporal_jitter = T, spatial_jitter = T)
      #dat$modd$similarity_cols <- kernels$spatial
      #dat$modd$similarity_rows <- kernels$temporal

      if(spatial == "none") dat$modd$similarity_cols = NULL
      if(temporal == "none") dat$modd$similarity_rows = NULL
      hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                             dat$modd$similarity_cols, 0, 0)
      hparam$laplace$step_sizes <- c(0.1, 0.01)
      hparam$laplace$min <- 0
      hparam$laplace$max <- 0.5
      hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
      hparam$rank$max <- 15

      for(intercept in with_intercept){

      res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                   trace=1, hpar=hparam, intercept_row = intercept,
                                   intercept_col = intercept, ls_initial = T,
                                   seed = seed, warm_start = NULL, maxit=600,
                                   shared_information = T, thresh=1e-6,
                                   num_cores = 0)
      s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
      s0.1 <- output_wrapper_bixi(res0$best_fit, dat,T, IMR:::error_metric$mape)
      all_res %<>% rbind (data.frame(temporal = temporal,
                                     spatial=spatial,
                                     model    = "IMR",
                                     val_size = vals,
                                     intercept = intercept,
                                     lambda_m = res0$best_fit$lambda_laplace,
                                     rank     = s0$res$rank_M,
                                     seed     = seed,
                                     error=s0$res$error.test,
                                     mape = s0.1$res$error.test))

      }
      #-- fit soft-impute
      # fit_si <- simpute.cv(
      #   y_full = dat$modd$Y,
      #   y_train = dat$modd$y_train,
      #   y_valid = dat$modd$y_valid,
      #   n.lambda = 80,
      #   trace=FALSE,
      #    maxit = 800,
      #   test_error = IMR:::error_metric$rmse
      # )
      # s1 <- output_wrapper_bixi(fit_si$fit, dat)
      # s1.0 <- output_wrapper_bixi(fit_si$fit, dat, T, IMR:::error_metric$mape)
      # all_res %<>% rbind (data.frame(temporal = temporal,
      #                                spatial = spatial,
      #                                model    = "Simpute",
      #                                val_size = vals,
      #                                intercept = NA,
      #                                lambda_m = fit_si$lambda,
      #                                rank     = s1$res$rank_M,
      #                                seed     = seed,
      #                                error=s1$res$error.test,
      #                                mape = s1.0$res$error.test))
      #---
      #print(all_res)
    }
    # all_res %>% arrange(error) %>% mutate(error = round(error,6),
    #                                       mape  = round(mape, 6)) %>% print()

    all_res %>%
      dplyr::select(-temporal, -spatial) %>%
      group_by(model, intercept, val_size) %>%
      dplyr::summarise_all(summ_res ) %>% print()

      print(i)
  }

  all_res %>% arrange(error) %>% mutate(error = round(error,6))
}
saveRDS(all_res,"./article_results/bixi/data/results_jan26/seed_vs_valsize_2.rds")
 summ_res <- function(x)
   paste0(round(median(x),4), " (", round(sd(x),4), ")[",round(min(x),4), ",",
          round(max(x),4), "]")

#--- read and analyze
readRDS("./article_results/bixi/data/results_jan26/seed_vs_valsize_1.rds") %>%
  mutate(model = ifelse(model == "IMR", "IMR Matern", model)) %>%
  rbind(
    readRDS("./article_results/bixi/data/results_jan26/seed_vs_valsize_2.rds") %>%
      mutate(model = ifelse(model == "IMR", "IMR Identity", model))
  ) -> all_res


all_res %>%
  dplyr::select(-temporal, -spatial) %>%
  group_by(model, intercept, val_size) %>%
  dplyr::summarise_all(mean) %>%
  ungroup() %>%
  arrange(error) %>%
  transmute(model, intercept, val_size, order = 1:n()) %>%
  left_join(all_res %>%
              dplyr::select(-temporal, -spatial) %>%
              group_by(model, intercept, val_size) %>%
              dplyr::summarise_all(summ_res) %>%
              ungroup(),
            c("model", "intercept", "val_size")) %>%
  arrange(order) %>%
  dplyr::select(-order, -seed, -lambda_m, -mape)
 #%>%


#----------------------------------------------------------------------------------------------
# 2)


{
  #B = c(5, 10, 20, 25)
  miss = total_miss
  temporal = "none"
  spatial = "none"
  n_seeds = 100
  val_size = c(.05, .1, .2, .3)
  max_val = max(val_size)

  all_res_s2 <- data.frame()
  for(i in 1:n_seeds){

    initialize_parallel_workers(9)
    seed = 2025 + i

    dat <- prepare_bixi_data(miss, timestamp, seed = seed,
                             val_prop = 0.3, temporal = temporal, spatial=spatial,
                             temporal_jitter = T, spatial_jitter = T,
                             kappa_max = 1e4, tau_max = 1e-2)
    dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL
    if(spatial == "none") dat$modd$similarity_cols = NULL
    if(temporal == "none") dat$modd$similarity_rows = NULL

    hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                           dat$modd$similarity_cols, 0, 0)
    hparam$laplace$step_sizes <- c(0.1, 0.01)
    hparam$laplace$min <- 0
    hparam$laplace$max <- 0.5
    hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
    hparam$rank$max <- 15

    validation_mask <- as.matrix(dat$modd$valid_mask)
    validation_set <- dat$modd$y_valid

    for(vals in val_size){
      if(vals != max_val){
        new_ratio = vals / max_val
        new_mask <- IMR::mask_train_test_split(validation_mask, new_ratio, seed)
        dat$modd$valid_mask <- new_mask
        dat$modd$y_valid <- validation_set * new_mask
      }

      for(intercept in with_intercept){

        res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                     trace=1, hpar=hparam, intercept_row = intercept,
                                     intercept_col = intercept, ls_initial = T,
                                     seed = seed, warm_start = NULL, maxit=600,
                                     shared_information = T, thresh=1e-6,
                                     num_cores = 0)
        s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
        s0.1 <- output_wrapper_bixi(res0$best_fit, dat,T, IMR:::error_metric$mape)
        all_res_s2 %<>% rbind (data.frame(temporal = temporal,
                                       spatial=spatial,
                                       val_size = vals,
                                       intercept = intercept,
                                       lambda_m = res0$best_fit$lambda_laplace,
                                       rank     = s0$res$rank_M,
                                       seed     = seed,
                                       error=s0$res$error.test,
                                       mape = s0.1$res$error.test))

      }

      #print(all_res_s2)
    }
    # all_res_s2 %>% arrange(error) %>% mutate(error = round(error,6),
    #                                       mape  = round(mape, 6)) %>% print()

    all_res_s2 %>%
      dplyr::select(-temporal, -spatial) %>%
      group_by(intercept,  val_size) %>%
      summarise_all(summ_res ) %>% print()

    print(i)
  }

  all_res_s2 %>% arrange(error) %>% mutate(error = round(error,6))
}
saveRDS(all_res_s2,"./article_results/bixi/data/results_jan26/fixedtrain_vs_valsize_3.rds")

summ_res <- function(x)
  paste0(round(mean(x),4), " (", round(sd(x),4), ")[",round(min(x),4), ",",
         round(median(x), 4), ",",
         round(max(x),4), "]")
readRDS("./article_results/bixi/data/results_jan26/fixedtrain_vs_valsize.rds_2") %>%
  # mutate(model = ifelse(model == "IMR", "IMR Matern", model)) %>%
  rbind(
    readRDS("./article_results/bixi/data/results_jan26/fixedtrain_vs_valsize_3.rds") #%>%
      # mutate(model = ifelse(model == "IMR", "IMR Identity", model))
  ) -> all_res2


all_res2 %>%
  #dplyr::select(-temporal, -spatial) %>%
  group_by(temporal, spatial, intercept, val_size) %>%
  dplyr::summarise_all(median) %>%
  ungroup() %>%
  arrange(error) %>%
  transmute(temporal, spatial, intercept, val_size, order = 1:n()) %>%
  left_join(all_res2 %>%
              #dplyr::select(-temporal, -spatial) %>%
              group_by(temporal, spatial, intercept, val_size) %>%
              dplyr::summarise_all(summ_res) %>%
              ungroup(),
            c("temporal", "spatial", "intercept", "val_size")) %>%
  arrange(order) %>%
  dplyr::select(-order, -seed, -lambda_m, -mape)
#=====================================================================


clean_res <- all_res2 |>
  # Create a distinct identifier for the 3 model variations
  mutate(
    model_id = case_when(
      temporal == "original" & spatial == "original" & intercept == TRUE  ~ "Orig/Orig (Int)",
      temporal == "original" & spatial == "original" & intercept == FALSE ~ "Orig/Orig (No Int)",
      temporal == "none" & spatial == "none" & intercept == FALSE     ~ "None/None (No Int)",
      TRUE ~ "Other" # Catch-all for safety
    ),
    # Ensure val_size is treated as a factor for boxplots, but numeric for lines
    val_size_fct = factor(val_size)
  ) |>
  # Filter to keep only the variations of interest
  filter(model_id != "Other")


summary_table <- clean_res |>
  summarise(
    mean_error = mean(error, na.rm = TRUE),
    sd_error   = sd(error, na.rm = TRUE),
    min_error  = min(error, na.rm = TRUE),
    max_error  = max(error, na.rm = TRUE),
    n_reps     = n(),
    .by = c(model_id, val_size) # Group by Model and Validation Size
  ) |>
  arrange(model_id, val_size)

# Print table
print(summary_table)

# Plot 1: Trend of Mean Error (with Confidence Intervals)
# This answers: Does changing validation size systematically shift the error?
p1 <- clean_res |>
  ggplot(aes(x = val_size, y = error, color = model_id)) +
  # Add jittered points to see raw data distribution
  geom_jitter(alpha = 0.2, width = 0.01) +
  # Add trend lines (mean)
  #geom_smooth(method = "loess", se = TRUE, alpha = 0.1) +
  labs(
    title = "Impact of Validation Size on Test RMSE",
    subtitle = "Comparing Model Variations",
    x = "Validation Proportion",
    y = "Test RMSE",
    color = "Model Variation"
  ) +
  theme_minimal();p1

# Plot 2: Stability Analysis (Boxplots)
# This answers: Does variance increase/decrease with validation size?
p2 <- clean_res |>
  ggplot(aes(x = val_size_fct, y = error, fill = model_id)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  facet_wrap(~model_id, scales = "fixed") +
  labs(
    title = "Error Distribution by Validation Size",
    x = "Validation Proportion",
    y = "Test RMSE",
    fill = "Model Variation"
  ) +
  theme_minimal() +
  theme(legend.position = "none"); p2

# Display plots
print(p1)
print(p2)

# Fit a linear model: Error ~ Model * Validation_Size
fit_mod <- lm(error ~ model_id * val_size + 1, data = clean_res)

# Check the ANOVA table to see significance of main effects and interactions
anova(fit_mod)
summary(fit_mod)
#=====================================================================
# Q: choosing val_prop
B = 30
miss=0.5
val_prop = c(0.05, 0.1, 0.2, 0.3, 0.4)
temporal_kernel = "RBF"
spatial_kernel = "Matern"
all_res2 <- data.frame()
initialize_parallel_workers(9)
for(i in 1:B){
  seed = 2025 + i
  for(j in 1:length(val_prop)){
    dat <- prepare_bixi_data(miss, "25Sep", seed = seed,
                             val_prop = val_prop[j])
    #dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL
    kernels <- generate_similarity_bixi(miss, "25Sep", temporal = temporal_kernel,
                                        spatial=spatial_kernel,
                                        temporal_jitter = T, spatial_jitter = T, ell = 1.3,
                                        kappa_max = 1e3)
    dat$modd$similarity_cols <- kernels$spatial
    dat$modd$similarity_rows <- kernels$temporal
    hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                           dat$modd$similarity_cols, 0, 0)
    hparam$laplace$step_sizes <- c(1, 0.1, 0.01)
    hparam$laplace$min <- 0
    hparam$laplace$max <- 5
    hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
    hparam$rank$max <- 15
    res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                 trace=1, hpar=hparam, intercept_row = T,
                                 intercept_col = T, ls_initial = T,
                                 seed = seed, warm_start = NULL, maxit=600,
                                 shared_information = T, thresh=1e-6,
                                 num_cores = 0)
    s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
    all_res2 %<>% rbind (data.frame(val_prop = val_prop[j], error=s0$res$error.test))
  print(j)
  }
  print(i)
}
saveRDS(all_res2, "./article_results/bixi/data/results_jan26/val_prop_results.rds")
all_res2 %>%
  group_by(val_prop) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  as.data.frame() %>%
  # mutate(error = round(error,5)) %>%
  arrange((error)) %>%
  mutate(
    val_prop = paste0(round(val_prop*100, 0),"%"),
    error = round(error, 4)
  ) %>%
  rename(`validation size` = val_prop) %>%
  arrange(error) -> tab2;tab2

min_err2 <- min(tab2$error, na.rm = TRUE)

datatable(tab2, rownames = FALSE, options = list(pageLength = nrow(tab2), dom = "t")) %>%
  formatStyle(
    "error",
    backgroundColor = styleEqual(min_err2, "khaki"),
    fontWeight = styleEqual(min_err2, "bold")
  )

#=================================================================================================
# choosing the correlation target
temporal = "RBF"
spatial = "Matern"
seed = 2025
initialize_parallel_workers(9)
dat <- prepare_bixi_data(miss, "25Sep", seed = seed,
                         val_prop = 0.05)
dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL
corrt = 0.08

kernels <- generate_similarity_bixi(miss, "25Sep", temporal = temporal, spatial=spatial,
                                    temporal_jitter = T, spatial_jitter = T, ell = 1.3,
                                    kappa_max = 1e3, cor.target = corrt)
dat$modd$similarity_cols <- kernels$spatial
dat$modd$similarity_rows <- kernels$temporal
hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                       dat$modd$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(1, 0.1, 0.01)
hparam$laplace$min <- 0
hparam$laplace$max <- 5
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
hparam$rank$max <- 15
res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                             trace=2, hpar=hparam, intercept_row = T,
                             intercept_col = T, ls_initial = T,
                             seed = seed, warm_start = NULL, maxit=600,
                             shared_information = T, thresh=1e-6,
                             num_cores = 0)
s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
print(corrt)
print(s0$res$error.test)
#=====================================================================
res.bktr <- BKTR_Bixi_Wrapper(dat, "25Sep", 2025, miss,return_fit = T)
res.bktr$results$error.test
#======================================================================
#--- main tests
# Intercepts, no covariates 2 scenarios
for(B in 1:100){


intercept = c(TRUE, FALSE)
covariates = c(NA, TRUE, FALSE)
temporal = "RBF"
spatial = "none"
seed = 2025 + B
corrt = 0.08
val_prop = 0.05


combined <- expand.grid(
  intercept = intercept,
  covariates  = covariates,
  stringsAsFactors = FALSE
)
all_res3 <- list()
all_res4 <- data.frame()
initialize_parallel_workers(9)
#for(i in 1:B){
for(j in 1:nrow(combined)){
  shared <- combined$covariates[j]
  dat <- prepare_bixi_data(miss, "25Sep", seed = seed, val_prop = val_prop)
  if(is.na(shared)){
    dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL
    shared = TRUE
  }
  kernels <- generate_similarity_bixi(miss, "25Sep", temporal = temporal, spatial=spatial,
                                      temporal_jitter = T, spatial_jitter = T, ell = 1.3,
                                      kappa_max = 1e3, cor.target = corrt)
  dat$modd$similarity_cols <- NULL#kernels$spatial
  dat$modd$similarity_rows <- kernels$temporal
  hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                         dat$modd$similarity_cols, 0, 0)
  hparam$laplace$step_sizes <- c(0.1, 0.01)
  hparam$laplace$min <- 0
  hparam$laplace$max <- 1
  hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
  hparam$rank$max <- 15
  hparam$beta$step_sizes <- hparam$gamma$step_sizes <-
    step_size <- function(min_val, max_val, n=60) {
      step = (max_val - min_val) / (n - 1L)
      print(step)
      return(step)
    }


  res0 <- IMR:::imr.cv(dat$modd,
                               trace=2, hpar=hparam, intercept_row = combined$intercept[j],
                               intercept_col = combined$intercept[j], ls_initial = T,
                               seed = seed, warm_start = NULL, maxit=600,
                               shared_information = shared, thresh=1e-6,
                               num_cores = 0)
  s0 <- output_wrapper_bixi(res0$best_fit, dat,shared)
  all_res3[j] <- list(intercept = combined$intercept[j], covariate=combined$covariates[j],res=s0)

  all_res4 %<>% rbind (data.frame(intercept = combined$intercept[j], covariate=combined$covariates[j],
                                  error=s0$res$error.test))
}
print(B)
}
s0$res
saveRDS(s0,"./article_results/bixi/data/results_jan26/last_results_list1.rds")
saveRDS(all_res4,"./article_results/bixi/data/results_jan26/last_results_table.rds")

all_res4 <- readRDS("./article_results/bixi/data/results_jan26/last_results_table.rds")
all_res4 %>%
  mutate(model = paste0(if_else(intercept, "Intercept+",""),
                        if_else(is.na(covariate), "",
                                if_else(covariate, "+covariate(shared)",
                                "+covariate(separate)")))) %>%
  select(model, error) %>%
  rbind( data.frame(model="BKTR", error=res.bktr$results$error.test)) %>%
  arrange(error)


library(dplyr)
library(DT)

tab3 <- all_res4 %>%
  mutate(model = paste0(
    if_else(intercept, "Intercept + ", ""),
    if_else(is.na(covariate), "",
            if_else(covariate, "covariate(shared)", "covariate(separate)"))
  )) %>%
  select(model, error) %>%
  bind_rows(tibble(model = "BKTR", error = res.bktr$results$error.test)) %>%
  mutate(
    model = if_else(model == "", "(none)", model),
    error = round(error, 4)
  ) %>%
  arrange(error)

min_err3 <- min(tab3$error, na.rm = TRUE)

datatable(tab3, rownames = FALSE, options = list(pageLength = nrow(tab3), dom = "t")) %>%
  formatStyle(
    "error",
    backgroundColor = styleEqual(min_err3, "khaki"),
    fontWeight = styleEqual(min_err3, "bold")
  )


#==================================================================
hparam$laplacian_row <- hparam$laplacian_col <- NULL
dat$modd$similarity_rows <- dat$modd$similarity_cols <- NULL

bench::bench_time(res <- IMR:::imr.cv(dat$modd,
                                              trace=2, hpar=hparam, intercept_row = T,
                                              intercept_col = T, ls_initial = T,
                                              seed = seed, warm_start = NULL, maxit=900,
                                              thresh = 1e-6,
                                              shared_information = T,
                                              num_cores = 0)) -> timee
round(lubridate::time_length(timee[2], "minutes"), 2)
s <- output_wrapper_bixi(res$best_fit, dat,T)
s$res$error.test

res$best_fit$gamma
res$best_fit$beta
s$rec$beta
s$rec$gamma


dat$modd$similarity_cols[1:5,1:5]


distance.mat = temporal_kernel$distance_matrix

distance.mat = spatial_kernel$distance_matrix %>% as.matrix()

distance.mat = fields::rdist(as.matrix(1:579))
kernel =   fields::Matern(distance.mat, smoothness = 3/2)

svd(K)$d

min(distance.mat[upper.tri(distance.mat)])
d <- distance.mat
dvec <- d[upper.tri(d)]
dvec <- dvec[dvec > 0]
range0 <- median(dvec)

K <- fields::Matern(d, smoothness = 3/2, range = range0)
tau2 <- 1e-3
K <- K + diag(tau2, nrow(K))
solve(K)[1:5,1:5]
Ki <- chol2inv(chol(K))
Ki[1:5,1:5]

e <- eigen(K, symmetric=TRUE, only.values=TRUE)$values
range(e)
kappa(K)
max(abs(Ki - t(Ki)))


distance.mat[1:10,1:10]
solve(kernel)[1:10,1:10]


e <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
range(e)
kappa(K)

spatial_kernel[upper.tri(spatial_kernel)] = t(spatial_kernel)[upper.tri(spatial_kernel)]

L = chol(K)

solve(L)[1:10,1:10]


mask_train_valid_split <- function(
    obs_mask,
    valid_p = 0.2,
    seed = NULL,
    min_train_row = 1,
    min_train_col = 1
) {
  stopifnot(is.matrix(obs_mask))
  stopifnot(is.numeric(valid_p), valid_p >= 0, valid_p <= 1)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  obs_mask <- obs_mask != 0
  n_rows <- nrow(obs_mask)
  n_cols <- ncol(obs_mask)

  obs_ij <- which(obs_mask, arr.ind = TRUE)
  n_obs <- nrow(obs_ij)
  if (n_obs == 0) {
    return(matrix(0L, nrow = n_rows, ncol = n_cols))
  }

  target_valid <- floor(n_obs * valid_p)

  row_left <- rowSums(obs_mask)
  col_left <- colSums(obs_mask)

  # Random order over observed cells
  ord <- sample.int(n_obs)
  obs_ij <- obs_ij[ord, , drop = FALSE]

  valid_linear <- integer(0)

  for (k in seq_len(n_obs)) {
    if (length(valid_linear) >= target_valid) break

    r <- obs_ij[k, 1]
    c <- obs_ij[k, 2]

    # Only hold out if the remaining training counts stay >= minimums
    if ((row_left[r] - 1L) >= min_train_row && (col_left[c] - 1L) >= min_train_col) {
      valid_linear <- c(valid_linear, r + (c - 1L) * n_rows)
      row_left[r] <- row_left[r] - 1L
      col_left[c] <- col_left[c] - 1L
    }
  }

  if (length(valid_linear) < target_valid) {
    warning(
      sprintf(
        "Could only select %d/%d validation cells (valid_p=%.3f) under min_train_row=%d, min_train_col=%d.",
        length(valid_linear), target_valid, valid_p, min_train_row, min_train_col
      ),
      call. = FALSE
    )
  }

  valid_mask <- matrix(0L, nrow = n_rows, ncol = n_cols)
  valid_mask[valid_linear] <- 1L
  valid_mask
}

obs <- as.matrix(dat$modd$obs_mask)
valid1 = mask_train_valid_split(as.matrix(dat$modd$obs_mask), 0.05, 2025, 5, 5)
apply( obs * (1-valid2), 1, sum) %>% summary()
apply( valid2, 2, sum) %>% summary()



valid2 = IMR::mask_train_test_split(as.matrix(dat$modd$obs_mask), 0.05, 2025)

