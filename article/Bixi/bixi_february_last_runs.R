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



generate_similarity_bixi(0.55, "25Sep", temporal = "none", spatial = "simulated",
                         temporal_jitter=T, spatial_jitter = T,
                         jitter_kappa_max = 1e3
                         ) -> ks


model_combn <- expand.grid(
  kernels = c("simulated", "none"),
  covariates = c(T,F),
  Intercepts = c(T, F),
  stringsAsFactors = FALSE
)
model_combn <- model_combn[c(7,8),]

#i = 6
# prefix=6

total_results <- data.frame()
total_results <- readRDS("./article/Bixi/data/IMR_results_final.rds")


for(rep in 1: 10){
seed = 4000 + rep
for(prefix in 1:50){
for(i in 1:nrow(model_combn)){


original_missing_pct = 0.1163287
test_pct <- 0.2
total_miss = test_pct + original_missing_pct -> miss_p


dat <- prepare_bixi_data(total_miss, "Jan", seed,
                         prefix = prefix,
                         val_prop = 0.2,
                        bktr_variables = TRUE,
                        temporal = model_combn$kernels[i],
                        spatial = model_combn$kernels[i],
                        temporal_jitter = TRUE,
                        spatial_jitter = TRUE,
                        jitter_kappa_max=1e3,
                        jitter_tau_max = 1e-2)

if(!model_combn$covariates[i])
  dat$X <- dat$Z <- dat$modd$Xq <- dat$modd$Xr <-
  dat$modd$Zq <- dat$modd$Zr <- NULL

hparams <- IMR::get_imr_default_hparams()
hparam$beta
hparam$rank
hparam$laplace

hparam$rank$min = 1

hparam$laplace$min = 0
hparam$laplace$max = 2
hparam$laplace$step_sizes = c(0.1)

hparam$beta$max = .5#0.5
hparam$gamma$max = .1#0.2
hparam$gamma$value = 0
hparam$beta$value = .17
hparam$beta$length = hparam$gamma$length = 10

start = Sys.time()
bench::bench_time(fitimr <- IMR:::imr.cv_21(dat$modd, intercept_row = model_combn$Intercepts[i],
                          intercept_col = model_combn$Intercepts[i],
                          hpar = hparam,
                          thresh = 1e-4, maxit = 800,
                          shared_information = TRUE,
                          final_thresh = 1e-6, final_maxit = 1000,
                          #init_thresh = 1e-4, init_maxit = 500,
                          trace = 2, num_cores = 9, seed = seed)) -> time.imr
time  = Sys.time()-start

s0 <- output_wrapper_bixi(fitimr$fit, dat,shared_information = T)
#s0$res

res <- data.frame(
  test = s0$res$error.test,
  time = time,
  model = paste0(model_combn$kernels[i],
                 ifelse(model_combn$covariates[i], "+covariates",""),
                 ifelse(model_combn$Intercepts[i], "+Intercept","")),
  total_miss = total_miss,
  prefix = prefix,
  test_pct = test_pct,
  lambda_beta = fitimr$fit$params$lambda_beta,
  lambda_gamma = fitimr$fit$params$lambda_gamma,
  lambda_laplace = fitimr$fit$params$lambda_laplace,
  rank_estim   = fitimr$fit$params$rank,
  rank_M  = s0$res$rank_M,
  rank_beta = s0$res$rank_beta,
  rank_gamma = s0$res$rank_gamma,
  sparsity_beta = s0$res$sparsity_beta,
  sparsity_gamma = s0$res$sparsity_gamma,
  time2.1 = as.numeric(time.imr[1]),
  time2.2 = as.numeric(time.imr[2]),
  metric = "RMSE",
  run_id = 3
)

s0 <- output_wrapper_bixi(fitimr$fit, dat,shared_information = T,
                          test_error = IMR:::error_metric$rel.rmse)
res <- rbind(
  res, data.frame(
    test = s0$res$error.test,
    time = time,
    model = paste0(model_combn$kernels[i],
                   ifelse(model_combn$covariates[i], "+covariates",""),
                   ifelse(model_combn$Intercepts[i], "+Intercept","")),
    total_miss = total_miss,
    prefix = prefix,
    test_pct = test_pct,
    lambda_beta = fitimr$fit$params$lambda_beta,
    lambda_gamma = fitimr$fit$params$lambda_gamma,
    lambda_laplace = fitimr$fit$params$lambda_laplace,
    rank_estim   = fitimr$fit$params$rank,
    rank_M  = s0$res$rank_M,
      rank_beta = s0$res$rank_beta,
    rank_gamma = s0$res$rank_gamma,
    sparsity_beta = s0$res$sparsity_beta,
    sparsity_gamma = s0$res$sparsity_gamma,
    time2.1 = as.numeric(time.imr[1]),
    time2.2 = as.numeric(time.imr[2]),
    metric = "RRMSE",
    run_id = 3
  )
)
print(res)
total_results <- rbind(total_results, res)
}
  saveRDS(total_results, "./article/Bixi/data/IMR_results_final.rds")
}
}
#================================================
total_results %>%
  filter(metric == "RRMSE") %>%
  group_by(model) %>%
  summarise_all(mean) %>%
  as.data.frame() %>%
  arrange(test) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  mutate(time = round(time, 2)) %>%
  filter(model %in% c("simulated+covariates",
                      "simulated",
                      "none"))
