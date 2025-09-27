library(tidyverse)
library(devtools)
load_all()
library(magrittr)
source("./article_results/bixi/preprocess_data.R")
source("./other_models/SoftImpute_cv.R")

preprocess_bixi_data(.99, "25Sep", 2025)

dat <- prepare_bixi_data(0.99, "25Sep")

dat$splits$test@x %>% length()


bktr.fit <- fit_BKTR_to_Bixi(0.6, "25Sep", 2025)
BKTR_Bixi_Wrapper(dat, "25Sep", 2025, 0.6)

hpar <- IMR::get_imr_default_hparams()
hpar$beta$n.lambda = 30
hpar$beta$lambda_max = 0.1
hpar$gamma$lambda_max = 0.1
hpar$gamma$n.lambda = 30
hpar$M$n.lambda = 80

for(miss in c(.6,.7,.8,.9,.95, .55,.65, .75, .85, .99)){
  dat <- prepare_bixi_data(miss, "25Sep")
  bktr.fit <- fit_BKTR_to_Bixi(miss, "25Sep", 2025)
future::plan(future::sequential)
future::plan(future::multisession, workers=9)

fit.imr <- IMR::imr.cv(Y = dat$Y.inc,
                       X = dat$Xq,
                       Z = dat$Zq,
                       intercept_row = T,
                       intercept_col = T ,
                       hpar = hpar,
                       verbose = 1,
                       seed = 2025)
saveRDS(fit.imr, paste0("./article_results/bixi/data/imr_",round(100*miss),"_fit.rds"))
}

total_res <- data.frame()
for(miss in c(.6,.7,.8,.9,.95, .55,.65, .75, .85, .99)){
  dat <- prepare_bixi_data(miss, "25Sep")

  res.bktr <- BKTR_Bixi_Wrapper(dat, "25Sep", 2025, miss)
  fit.imr <- readRDS(paste0("./article_results/bixi/data/imr_",round(100*miss),"_fit.rds"))
  out <- IMR:::reconstruct(fit.imr$fit, dat)
  prepare_output_bixi(
    NULL,
    X           = dat$X,
    estim.test  = out$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test    = dat$splits$test@x,
    obs.train   = dat$Y.inc@x,
    beta.estim  = out$beta,
    M.estim     = out$M,
    time_per_fit = fit.imr$time_per_fit,
    total_num_fits = fit.imr$total_num_fits
  ) -> res.imr
  total_res <- rbind(total_res,
  res.imr[-10] %>%
    rbind(res.bktr[-c(13,1:3)]) %>%
    cbind(miss=miss,model=c("IMR", "BKTR"))
                     )


}

total_res %>%
  as.data.frame() %>%
  mutate(error.test = as.numeric(error.test),
         corr.test = as.numeric(corr.test),
         error.train = as.numeric(error.train),
         miss = as.numeric(miss),
         sparsity = as.numeric(sparsity)) %>%
  group_by(miss) %>%
  mutate(diff = error.test[1] - error.test[2]) %>%
  ungroup() %>%
  as.data.frame() %>%
  dplyr::arrange(miss, model) %>%
  dplyr::select(-time, -time_per_fit, -total_num_fits) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

out <- IMR:::reconstruct(fit.imr$fit, dat)

prepare_output_bixi(
  time        = NULL,
  X           = dat$Z,
  estim.test  = out$estimates[as.matrix(dat$splits$test != 0)],
  estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
  obs.test    = dat$splits$test@x,
  obs.train   = dat$Y.inc@x,
  beta.estim  = out$gamma,
  M.estim     = out$M,
  time_per_fit = fit.imr$time_per_fit,
  total_num_fits = fit.imr$total_num_fits
)





fit <- simpute.cv(
  y_train   = dat$splits$train,
  y_valid   = dat$splits$valid,
 # W_valid   = dat$masks$valid,
  y_full =  dat$Y.inc,
  n.lambda  = hpar$M$n.lambda,
  trace     = FALSE,
  print.best= FALSE,
  tol       = 5,
  thresh    = 1e-6,
  rank.init = hpar$M$rank.init,
  rank.limit= hpar$M$rank.limit,
  rank.step = hpar$M$rank.step,
  maxit     = 600,
  seed      = NULL
)
fit <- IMR:::reconstruct(fit$fit, dat)
  prepare_output_bixi(
    time        = fit$time,
    X           = NULL,
    estim.test  = fit$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = fit$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test    = dat$splits$test@x,
    obs.train   = dat$Y.inc@x,
    M.estim     = fit$estimates,
    total_num_fits = fit$total_num_fits,
    time_per_fit = fit$time_per_fit
  )


# BKTR .0870
# rows .0877
# both .0876
# interc .0878
# none  .088
#Simpute .101


