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

for (prefix in 1:10) {
  preprocess_bixi_data(total_miss, "Feb_last", 2025, prefix,
    file_override = FALSE,
    decreasing_train = TRUE, create_folder = TRUE,
    train_n_steps = 5, train_stepsize = .05,
    out_dir = "./article/bixi/data/splits2/"
  )
}




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
#====================================================================

for (rep in 1:10) {
  seed <- 4000 + rep
  for (prefix in 1:10) {
    for (train_size in train_seq) {
      for (i in 1:nrow(model_combn)) {
        dat <- prepare_bixi_data(total_miss, "Feb_last", seed,
          prefix = prefix,
          train_prefix = train_size,
          val_prop = 0.2,
          bktr_variables = TRUE,
          file_dir = "./article/bixi/data/splits2/",
          temporal = model_combn$kernels[i],
          spatial = model_combn$kernels[i],
          temporal_jitter = TRUE,
          spatial_jitter = TRUE,
          jitter_kappa_max = 1e3,
          jitter_tau_max = 1e-2
        )

        if (!model_combn$covariates[i]) {
          dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <-
            dat$modd$Zq <- dat$modd$Zr <- NULL
        }

        hparam <- IMR::get_imr_default_hparams()
        hparam$beta
        hparam$rank
        hparam$laplace

        hparam$rank$min <- 1

        hparam$laplace$min <- 0
        hparam$laplace$max <- 2
        hparam$laplace$step_sizes <- c(0.1)


        start <- Sys.time()
        bench::bench_time(fitimr <- IMR:::imr.cv_21(dat$modd,
          intercept_row = model_combn$Intercepts[i],
          intercept_col = model_combn$Intercepts[i],
          hpar = hparam,
          thresh = 1e-4, maxit = 800,
          shared_information = TRUE,
          final_thresh = 1e-6, final_maxit = 1000,
          # init_thresh = 1e-4, init_maxit = 500,
          trace = 1, num_cores = 9, seed = seed
        )) -> time.imr
        time <- Sys.time() - start

        s0 <- output_wrapper_bixi(fitimr$fit, dat, shared_information = T)
        # s0$res

        res <- data.frame(
          test = s0$res$error.test,
          time = time,
          model = paste0(
            model_combn$kernels[i],
            ifelse(model_combn$covariates[i], "+covariates", ""),
            ifelse(model_combn$Intercepts[i], "+Intercept", "")
          ),
          total_miss = total_miss,
          prefix = prefix,
          test_pct = test_pct,
          #lambda_beta = fitimr$fit$params$lambda_beta,
          #lambda_gamma = fitimr$fit$params$lambda_gamma,
          lambda_laplace = fitimr$fit$params$lambda_laplace,
          rank_estim = fitimr$fit$params$rank,
          rank_M = s0$res$rank_M,
          #rank_beta = s0$res$rank_beta,
          #rank_gamma = s0$res$rank_gamma,
          #sparsity_beta = s0$res$sparsity_beta,
          #sparsity_gamma = s0$res$sparsity_gamma,
          time2.1 = as.numeric(time.imr[1]),
          time2.2 = as.numeric(time.imr[2]),
          metric = "RMSE",
          train_size = train_size
        )

        s0 <- output_wrapper_bixi(fitimr$fit, dat,
          shared_information = T,
          test_error = IMR:::error_metric$rel.rmse
        )
        res <- rbind(
          res, data.frame(
            test = s0$res$error.test,
            time = time,
            model = paste0(
              model_combn$kernels[i],
              ifelse(model_combn$covariates[i], "+covariates", ""),
              ifelse(model_combn$Intercepts[i], "+Intercept", "")
            ),
            total_miss = total_miss,
            prefix = prefix,
            test_pct = test_pct,
            #lambda_beta = fitimr$fit$params$lambda_beta,
            #lambda_gamma = fitimr$fit$params$lambda_gamma,
            lambda_laplace = fitimr$fit$params$lambda_laplace,
            rank_estim = fitimr$fit$params$rank,
            rank_M = s0$res$rank_M,
            #rank_beta = s0$res$rank_beta,
            #rank_gamma = s0$res$rank_gamma,
            #sparsity_beta = s0$res$sparsity_beta,
            #sparsity_gamma = s0$res$sparsity_gamma,
            time2.1 = as.numeric(time.imr[1]),
            time2.2 = as.numeric(time.imr[2]),
            metric = "RRMSE",
            train_size = train_size
          )
        )
        print(res)
        total_results <- rbind(total_results, res)
      }
      saveRDS(total_results, "./article/Bixi/data/IMR_results_final2.rds")
    }
  }
}

total_results %>%
  filter(metric == "RRMSE") %>%
  group_by(model, train_size, metric) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(train_size, test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time = round(time, 2))

# =================================================================
# we now fit the same on BKTR
  seed <- 4000
  for (prefix in 5:10) {
    for (train_size in train_seq) {
      bktr_out <- fit_BKTR_to_Bixi(total_miss,
                                   "Feb_last",
                                   prefix = prefix,
                                   train_prefix = train_size,
                                   file_dir = "./article/bixi/data/splits2/",
                                   seed = seed,
                                   burn_in_iter = 1000,
                                   sampling_iter = 500)
      print(train_size)
    }
    print(paste0("prefix=",prefix))
  }
#==============================================================
