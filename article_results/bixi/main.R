library(tidyverse)
library(devtools)
load_all()
library(magrittr)
source("./article_results/bixi/preprocess_data.R")
source("./other_models/SoftImpute_cv.R")

 preprocess_bixi_data(.6, "25Sep", 2025)

dat <- prepare_bixi_data(0.6, "25Sep")
bktr.fit <- fit_BKTR_to_Bixi(0.6, "25Sep", 2025)


hpar <- IMR::get_imr_default_hparams()
hpar$beta$n.lambda = 60
hpar$gamma$n.lambda = 60
hpar$M$n.lambda = 80

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


