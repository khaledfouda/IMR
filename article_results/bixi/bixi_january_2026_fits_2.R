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
      fit_si <- simpute.cv(
        y_full = dat$modd$Y,
        y_train = dat$modd$y_train,
        y_valid = dat$modd$y_valid,
        n.lambda = 80,
        trace=FALSE,
         maxit = 800,
        test_error = IMR:::error_metric$rmse
      )
      s1 <- output_wrapper_bixi(fit_si$fit, dat)
      s1.0 <- output_wrapper_bixi(fit_si$fit, dat, T, IMR:::error_metric$mape)
      all_res %<>% rbind (data.frame(temporal = temporal,
                                     spatial = spatial,
                                     model    = "Simpute",
                                     val_size = vals,
                                     intercept = NA,
                                     lambda_m = fit_si$lambda,
                                     rank     = s1$res$rank_M,
                                     seed     = seed,
                                     error=s1$res$error.test,
                                     mape = s1.0$res$error.test))
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
  #saveRDS(all_res_s2,"./article_results/bixi/data/results_jan26/fixedtrain_vs_valsize_3.rds")
}

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
  dplyr::summarise_all(mean) %>%
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
# analyze the results above
{
clean_res <- all_res |>
  mutate(
    # Create a clear label for each model variation
    model_id = case_when(
      model == "IMR Matern" & intercept == TRUE  ~ "IMR Matern (Int)",
      model == "IMR Matern" & intercept == FALSE ~ "IMR Matern (No Int)",
      model == "IMR Identity" ~ "IMR Identity (No Int)", # Handling if it exists in full data
      model == "Simpute" ~ "Simpute",
      TRUE ~ paste(model, intercept) # Catch-all
    ),
    val_size_fct = factor(val_size)
  )

# Preview grouping
clean_res |> count(model_id)

# Plot: Error Distribution by Validation Size
clean_res |>
  ggplot(aes(x = val_size_fct, y = error, fill = model_id)) +
  # Boxplot shows the IQR (spread due to seed)
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.5) +
  labs(
    title = "Impact of Random Seed on Model Error",
    subtitle = "Spread (Box Height) represents variance due to different Train/Valid splits",
    x = "Validation Proportion",
    y = "Test RMSE",
    fill = "Model Variation"
  ) +
  theme_minimal()

# Calculate stability metrics across the 100 seeds
seed_analysis <- clean_res |>
  summarise(
    mean_error = mean(error, na.rm = TRUE),
    sd_error   = sd(error, na.rm = TRUE), # This is the "Seed Effect" magnitude
    cv_percent = (sd(error) / mean(error)) * 100, # Coefficient of Variation
    .by = c(model_id, val_size)
  ) |>
  arrange(val_size, mean_error)

print(seed_analysis)


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
}
#=====================================================================
#' 1)  why is the intercept worse? not sure
#' 2)  do we have the best kernels? > original/none
#' 3)  why are not much better than soft-impute
#' 4) what about the covariates?
#' --------------------------------------------------------------------
#' data setup:
all_runs <- data.frame(); counter = 1;

timestamp = "Jan"
original_missing_pct = 0.1163287
test_pct <- 0.2
total_miss = test_pct + original_missing_pct -> miss_p
seed = 2026
#for (cc in c(.05, .1, .2, .3, .4, .5)){
{
intercept = FALSE
covariate = TRUE
shared = FALSE
temporal = "original"
spatial = "none"
vals = 0.2
prefix = 15
for(seed in  1:10){
for(prefix in c(1:8, 14, 15)){
  initialize_parallel_workers(9)

#for(shared in c(T,F)){
#  for(intercept in c(T,F)){
if(seed == 4){
  if(prefix %in% 1:1)
    next
}

dat <- prepare_bixi_data(total_miss, timestamp, seed = seed, prefix = prefix,
                        val_prop = vals, temporal = temporal, spatial=spatial,
                        temporal_jitter = T, spatial_jitter = T,
                        bktr_variables = TRUE,
                        kappa_max = 1e3, tau_max = 1e-2)
if(spatial == "none") dat$modd$similarity_cols = NULL
if(temporal == "none") dat$modd$similarity_rows = NULL

if(!covariate)
  dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <- dat$modd$Zq <- dat$modd$Zr <- NULL

hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                       dat$modd$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(0.1, 0.01)
hparam$laplace$min <- 0
hparam$laplace$max <- 2
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
hparam$rank$max <- 15

hparam$beta$max = .5
hparam$gamma$max = .3

hparam$beta$step_sizes = step_size(0, 0.5, 40)
hparam$gamma$step_sizes = step_size(0, 0.3, 40)
# hparam$beta$step_sizes <- hparam$gamma$step_sizes <-
#   step_size <- function(min_val, max_val, n=40) {
#     step = (max_val - min_val) / (n - 1L)
#     print(step)
#     return(step)
#   }


ttime <- bench::bench_time( res0 <- IMR:::imr.cv(dat$modd,
                             trace=2, hpar=hparam, intercept_row = intercept,
                             intercept_col = intercept, ls_initial = T,
                             seed = seed, warm_start = NULL, maxit=800,
                             shared_information = shared, thresh=1e-6,
                             num_cores = 0))

s0 <- output_wrapper_bixi(res0$best_fit, dat,shared)
s0.1 <- output_wrapper_bixi(res0$best_fit, dat,shared, IMR:::error_metric$mape)
all_runs %<>% rbind(data.frame(temporal = temporal,
               spatial=spatial,
               model    = "IMR",
               val_size = vals,
               intercept = intercept,
               covariate = covariate,
               shared    = shared,
               lambda_m = res0$best_fit$lambda_laplace,
               rank     = s0$res$rank_M,
               prefix  = prefix,
               seed     = seed, #corrtarg = cc,
               error=s0$res$error.test,
               counter = counter,
               time1 = lubridate::time_length(ttime[1], "seconds"),
               time2 = lubridate::time_length(ttime[2], "seconds"),
               mape = s0.1$res$error.test)); counter=counter+1
all_runs %>% print()
}
saveRDS(all_runs, "./article_results/bixi/data/results_jan26/imrfit_20testsplits.rds")
}
}
#=======================
#' we will test the results from BKTR to see what we get
all_runs <- data.frame(); counter = 1;

timestamp = "Jan"
original_missing_pct = 0.1163287
test_pct <- 0.2
total_miss = test_pct + original_missing_pct -> miss_p -> miss
seed = 2026
prefix = 1
bktr_res <- BKTR_Bixi_Wrapper(dat, timestamp, seed, total_miss, T, prefix)
print(bktr_res$results$error.test)
bktr_res$results$time / 60
#for (cc in c(.05, .1, .2, .3, .4, .5)){

  intercept = TRUE
  covariate = TRUE
  shared = TRUE
  temporal = "original"
  spatial = "none"
  vals = 0.2

  dat <- prepare_bixi_data(total_miss, timestamp, seed = seed,prefix = prefix,
                           val_prop = vals, temporal = temporal, spatial=spatial,
                           temporal_jitter = T, spatial_jitter = T,
                           bktr_variables = TRUE,
                           kappa_max = 1e3, tau_max = 1e-2, cor.target=cc)







