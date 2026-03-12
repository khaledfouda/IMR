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

generate_data <- FALSE
if (generate_data) {
  for (prefix in 11:50) {
    preprocess_bixi_data(total_miss, "Feb_last", 2025+prefix, prefix,
      file_override = T,
      decreasing_train = TRUE, create_folder = TRUE,
      train_n_steps = 5, train_stepsize = .05,
      out_dir = "./article/bixi/data/splits2/"
    )
  }
}


# ================================================================
# we now train >>

model_combn <- expand.grid(
  similarity = c(TRUE, FALSE),
  covariates = F, #c(T, F),
  Intercepts = F, #c(T, F),
  stringsAsFactors = FALSE
)
#model_combn <- model_combn[c(7, 8), ]
train_seq <- round(seq(1 - total_miss, by = -.05, length.out = 5) * 100)
# ====================================================================
#rep = 1; prefix = 1; train_size = train_seq[1]; i=1
#---
convergence <- IMR::imr_convergence(maxit = 600, thresh=1e-5, trace=FALSE, ls_initial = TRUE)
grid <- IMR::imr_tune_grid(laplace = c(0,10,80, 5), rank = c(5,30,1, 5)); grid

#rep = 1; prefix = 1; i = 1; train_size=train_seq[1]

total_results <- all_res <-  data.frame()
for (rep in 1:10) {
  seed <- 4000 + rep
  for (prefix in 1:50) {
    for (train_size in train_seq) {
      for (i in 1:nrow(model_combn)) {
        dat <- prepare_bixi_data(total_miss, "Feb_last", seed,
          prefix = prefix,
          train_prefix = train_size,
          val_prop = 0.2,
          bktr_variables = TRUE,
          file_dir = "./article/bixi/data/splits2/",
          temporal = ifelse(model_combn$similarity[i],"simulated", "none"),
          spatial = ifelse(model_combn$similarity[i],"simulated", "none"),
          temporal_jitter = TRUE,
          spatial_jitter = TRUE,
          jitter_kappa_max = 1e3,
          jitter_tau_max = 1e-2
        )
        model_data <- dat$modd
        #print(model_data)

        model_data <- update(model_data,
                             row_covariates = model_combn$covariates[i],
                             col_covariates = model_combn$covariates[i],
                             shared_beta = TRUE, shared_gamma = TRUE,
                             intercept_row = model_combn$Intercepts[i],
                             intercept_col = model_combn$Intercepts[i],
                             row_similarity = model_combn$similarity[i],
                             col_similarity = model_combn$similarity[i]); print(model_data)

        grid <- IMR::imr_set_grid_limits(model_data, grid, default_rank=2, convergence=convergence,
                                         verbose=1,bisection_iter = 5); grid
        start <- Sys.time()
        bench::bench_time(fitimr <- IMR::imr_tune(model_data, grid, final_fit = TRUE,
                                                  fast_laplace = FALSE,
                                                  laplace_log_scale = FALSE,
                                                  convergence=convergence, n_cores=7,
                                                  seed = seed, verbose=1)) -> time.imr
        time <- Sys.time() - start

        s0 <- output_wrapper_bixi(fitimr$fit, dat, shared_information = T)

        res <- data.frame(
          test = s0$res$error.test,
          time = time,
          model = paste0(
            ifelse(model_combn$similarity[i], "similarity", "original"),
            ifelse(model_combn$covariates[i], "+covariates", ""),
            ifelse(model_combn$Intercepts[i], "+Intercept", "")
          ),
          total_miss = total_miss,
          prefix = prefix,
          test_pct = test_pct,
          #lambda_beta = fitimr$fit$meta$lambdas["beta"],
          #lambda_gamma = fitimr$fit$meta$lambdas["gamma"],
          lambda_laplace = fitimr$fit$meta$lambdas["M"],
          rank_estim = fitimr$fit$meta$rank,
          rank_M = s0$res$rank_M,
          # rank_beta = s0$res$rank_beta,
          # rank_gamma = s0$res$rank_gamma,
          # sparsity_beta = s0$res$sparsity_beta,
          # sparsity_gamma = s0$res$sparsity_gamma,
          time2.1 = as.numeric(time.imr[1]),
          time2.2 = as.numeric(time.imr[2]),
          metric = "RMSE",
          train_size = train_size
        );res

        s0 <- output_wrapper_bixi(fitimr$fit, dat,
          shared_information = T,
          test_error = IMR:::error_metrics$rrmse
        )
        res <- rbind(
          res, data.frame(
            test = s0$res$error.test,
            time = time,
            model = paste0(
              ifelse(model_combn$similarity[i], "similarity", "original"),
              ifelse(model_combn$covariates[i], "+covariates", ""),
              ifelse(model_combn$Intercepts[i], "+Intercept", "")
            ),
            total_miss = total_miss,
            prefix = prefix,
            test_pct = test_pct,
            # lambda_beta = fitimr$fit$params$lambda_beta,
            # lambda_gamma = fitimr$fit$params$lambda_gamma,
            lambda_laplace = fitimr$fit$meta$lambdas["M"],
            rank_estim = fitimr$fit$meta$rank,
            rank_M = s0$res$rank_M,
            # rank_beta = s0$res$rank_beta,
            # rank_gamma = s0$res$rank_gamma,
            # sparsity_beta = s0$res$sparsity_beta,
            # sparsity_gamma = s0$res$sparsity_gamma,
            time2.1 = as.numeric(time.imr[1]),
            time2.2 = as.numeric(time.imr[2]),
            metric = "RRMSE",
            train_size = train_size
          )
        )
        print(res)
        total_results <- rbind(total_results, res)
      }
      saveRDS(total_results, "./article/Bixi/data/final_results/IMR_results_final_25pct_2_4.rds")
    }
    print(paste0("Rep:", rep, ", prefix: ", prefix))
  }
}

total_results %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time = round(time, 2))

# =================================================================
# we now fit the same on BKTR
seed <- 4000
fit_bktr <- FALSE
if(fit_bktr){
  for (prefix in 1:50) {
    for (train_size in train_seq) {
      bktr_out <- fit_BKTR_to_Bixi(total_miss,
        "Feb_last",
        prefix = prefix,
        train_prefix = train_size,
        file_dir = "./article/bixi/data/splits2/",
        seed = seed,
        burn_in_iter = 1000,
        sampling_iter = 500
      )
      print(train_size)
    }
    print(paste0("prefix=", prefix))
  }
}
# ==============================================================
# we now fit IMR with fixed hyperparameters >>
seed <- 4000
convergence <- IMR::imr_convergence(maxit = 2000, thresh=1e-7, trace=FALSE, ls_initial = TRUE)

total_results %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, prefix, metric) %>%
  slice_min(test, n=1) %>%
  ungroup() %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(test) %>%
  #mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  #mutate(time = round(time, 2)) %>%
  group_by(model) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  select(model, lambda_laplace, rank_estim) %>%
  transmute(similarity = model == "similarity",
            lambda=lambda_laplace,
            rank = round(rank_estim)) ->
  hps;hps
  # group_by()

model_combn %<>%
  left_join(hps, "similarity")
#
# model_combn %<>%
#   mutate(
#     lambda = if_else(similarity, 1.2140, 0.7456),
#     rank = if_else(similarity, 12, 10)
#   )
total_results <- data.frame()
for (prefix in 1:50) {
  for (train_size in train_seq) {
    for (i in 1:nrow(model_combn)) {
      dat <- prepare_bixi_data(total_miss, "Feb_last", seed,
        prefix = prefix,
        train_prefix = train_size,
        val_prop = 0.2,
        bktr_variables = TRUE,
        file_dir = "./article/bixi/data/splits2/",
        temporal = ifelse(model_combn$similarity[i],"simulated", "none"),
        spatial = ifelse(model_combn$similarity[i],"simulated", "none"),
        temporal_jitter = TRUE,
        spatial_jitter = TRUE,
        jitter_kappa_max = 1e3,
        jitter_tau_max = 1e-2
      )
      model_data <- dat$modd
      model_data <- update(model_data,
                           row_covariates = model_combn$covariates[i],
                           col_covariates = model_combn$covariates[i],
                           intercept_row = model_combn$Intercepts[i],
                           intercept_col = model_combn$Intercepts[i],
                           row_similarity = model_combn$similarity[i],
                           col_similarity = model_combn$similarity[i]); print(model_data)


      start <- Sys.time()

      bench::bench_time(fitimr <- IMR::imr_fit(
        data = model_data,
        rank = model_combn$rank[i],
        lambda_m = model_combn$lambda[i],
        lambda_beta = 0,
        lambda_gamma = 0,
        convergence = convergence
      )) -> time.imr

      time <- Sys.time() - start

      s0 <- output_wrapper_bixi(fitimr, dat, shared_information = T)
      # s0$res

      res <- data.frame(
        test = s0$res$error.test,
        time = time,
        model = paste0(
          ifelse(model_combn$similarity[i], "similarity", "original"),
          ifelse(model_combn$covariates[i], "+covariates", ""),
          ifelse(model_combn$Intercepts[i], "+Intercept", "")
        ),
        total_miss = total_miss,
        prefix = prefix,
        test_pct = test_pct,
        lambda_laplace = model_combn$lambda[i],
        rank_estim = length(fitimr$coefficients$d > 0),
        rank_M = length(fitimr$coefficients$d > 0),
        time2.1 = as.numeric(time.imr[1]),
        time2.2 = as.numeric(time.imr[2]),
        metric = "RMSE",
        train_size = train_size
      )

      s0 <- output_wrapper_bixi(fitimr, dat,
        shared_information = T,
        test_error = IMR:::error_metrics$rrmse
      )
      res <- rbind(
        res, data.frame(
          test = s0$res$error.test,
          time = time,
          model =paste0(
            ifelse(model_combn$similarity[i], "similarity", "original"),
            ifelse(model_combn$covariates[i], "+covariates", ""),
            ifelse(model_combn$Intercepts[i], "+Intercept", "")
          ),
          total_miss = total_miss,
          prefix = prefix,
          test_pct = test_pct,
          lambda_laplace = model_combn$lambda[i],
          rank_estim = length(fitimr$coefficients$d > 0),
          rank_M = length(fitimr$coefficients$d > 0),
          time2.1 = as.numeric(time.imr[1]),
          time2.2 = as.numeric(time.imr[2]),
          metric = "RRMSE",
          train_size = train_size
        )
      )
      print(res)
      total_results <- rbind(total_results, res)
    }
    saveRDS(total_results, "./article/Bixi/data/final_results/IMR_results_onefit_25pct_2_4.rds")
  }
}
