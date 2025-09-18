library(devtools)
load_all()

require(tidyverse)
source("./notes/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")
source("./other_models/MCCI.R")


sim1_res <- function(dat, fit, ortho=TRUE, error_metric=IMR::error_metric$rmse){
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  has.beta <- "beta" %in% names(fit)
  has.M    <- "M" %in% names(fit) | all(c("u","d","v") %in% names(fit))
  has.intercept <- "beta0" %in% names(fit) & !is.null(fit$beta0)
  estimates <- 0
  out <- data.frame(beta=NA, M=NA, theta=NA, test=NA, rank=NA)
  # check M
  if(all(c("u", "d", "v") %in% names(fit)))
    fit$M <- fit$u %*% (fit$d * t(fit$v))
  if(has.M) {
    estimates <- fit$M
    out$M <- error_metric(fit$M, dat$M)
  }
  if(has.beta){
    fit$beta <- solve(dat$fit_data$X$R) %*% fit$beta
    estimates <- estimates + dat$X %*% fit$beta
    out$beta <- error_metric(fit$beta, dat$beta)
  }
  if(has.intercept){
    estimates <- estimates + fit$beta0 %*% matrix(1,1,ncol(dat$Y))
  }


  stopifnot(all(estimates!=0))
  out$theta <- error_metric(estimates, dat$theta)
  test.obs <- dat$Y == 0
  out$test <- IMR::error_metric$rel.rmse(estimates[test.obs], dat$theta[test.obs])
  out$rank <- qr(estimates)$rank
  return(out)
}




dat <-
  generate_simulated_data(1000, 1000, 10, 20, 0, 0.8,
                          sparsity_beta = .5, sparsity_gamma = 0,
                          prepare_for_fitting = T,mv_coeffs = T,seed = 2025)

future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()

# hpar$M$lambda_factor <- 1
# hpar$M$n.lambda <- 20
# hpar$M$rank.init <- 2
# hpar$M$rank.step <- 2
# hpar$M$rank.max  <- 15
# hpar$M$lambda_max <- 80
# hpar$beta$n.lambda <- 10
# hpar$gamma$n.lambda <- 10
fit.imr <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                     Z = dat$fit_data$Z$Q,intercept_row = F,
                     hpar = hpar, seed = 2025,
                     intercept_col = F, verbose=2)
sim1_res(dat, fit.imr$fit)

fit.imri <-  IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                         Z = NULL,intercept_row = T,
                         hpar = hpar, seed = 2025,
                         intercept_col = F, verbose=2)
sim1_res(dat, fit.imri$fit)

fit.nlrr <- IMR::imr.fit_no_low_rank()


quick_camc_simu_res(dat, fit32$fit)



fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                    trace=T#tol = 5, n.lambda=80, rank.init = 2,
                    #rank.step = 1,
                    #lambda0_fun = IMR::get_lambda_M_max
)
quick_camc_simu_res(dat, fitsi$fit)



sim1_res(dat, fit32$fit)
sim1_res(dat, o)



fitmmci <- MCCI.cv(dat$Y, dat$X, dat$mask, numCores = 9)
quick_camc_simu_res(dat, fitmmci$fit, T)
