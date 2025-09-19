library(devtools)
load_all()

require(tidyverse)
library(magrittr)
source("./notes/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")
source("./other_models/MCCI.R")
source("./other_models/Naive.R")

sim1_res <- function(dat, fit, name="", ortho=TRUE, error_metric=IMR::error_metric$rmse){
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  has.beta <- "beta" %in% names(fit) & ! is.null(fit$beta)
  has.gamma <- "gamma" %in% names(fit) & ! is.null(fit$gamma)
  has.M    <- "M" %in% names(fit) | all(c("u","d","v") %in% names(fit))
  has.intercept <- "beta0" %in% names(fit) & !is.null(fit$beta0)
  estimates <- 0
  out <- data.frame(model=name, beta=NA, gamma=NA, M=NA, theta=NA, test=NA, rank=NA)
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
  if(has.gamma){
    fit$gamma <- fit$gamma %*% solve(t(dat$fit_data$Z$R))
    estimates <- estimates + fit$gamma %*% t(dat$Z)
    out$gamma <- error_metric(fit$gamma, dat$gamma)
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


one_loop_fit <- function(dat, hpar = IMR::get_imr_default_hparams(), seed=2025){


 # fit.imr <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
 #                      Z = dat$fit_data$Z$Q,intercept_row = F,
 #                      hpar = hpar, seed = seed, ls_initial = FALSE,
 #                      intercept_col = F, verbose=1)

  fit.nlrr <-  IMR:::nlrr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                           Z = dat$fit_data$Z$Q,intercept_row = F,
                           hpar = hpar, seed = seed,
                           intercept_col = F, verbose=1)

  fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                      trace=T,
                      tol = hpar$M$early.stopping, seed = seed,
                      n.lambda = hpar$M$n.lambda
                      #, n.lambda=80, rank.init = 2,
                      #rank.step = 1,
                      #lambda0_fun = IMR::get_lambda_M_max
  )
  fitmmci <- MCCI.cv(dat$Y, dat$X, dat$mask, numCores = 9, seed = seed)
  fitnai <- naive.fit(dat$Y, dat$X)

  rbind(
  #sim1_res(dat, fit.imr$fit, "IMR"),
  sim1_res(dat, fit.nlrr$fit, "IMR"),
  sim1_res(dat, fitsi$fit, "SoftImpute"),
  sim1_res(dat, fitmmci$fit, "MCCI"),
  sim1_res(dat, fitnai, "Naive")
  )
}





future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()

# hpar$M$lambda_factor <- 1
hpar$M$n.lambda <- 30
# hpar$M$rank.init <- 2
# hpar$M$rank.step <- 2
# hpar$M$rank.max  <- 15
# hpar$M$lambda_max <- 80
# hpar$beta$n.lambda <- 10
# hpar$gamma$n.lambda <- 10
hpar$M$early.stopping <- 3

results4 <- data.frame()

for(i in 1:30){
  seed = 2025 + i
  dat <-
    generate_simulated_data(400, 400, 10, 20, 0, 0.8,
                            sparsity_beta = .5, sparsity_gamma = 0.0,
                            prepare_for_fitting = T,mv_coeffs = T,seed = seed)
  start <- Sys.time()
  results4 %<>% rbind(
  one_loop_fit(dat, hpar, seed)
  )
  results4 %>%
    dplyr::group_by(model) %>%
    dplyr::summarise_all(mean) %>%
    print()
  print(Sys.time() - start)
  print(i)

}

saveRDS(results4, "./notes/saved/results4400.rds")

results4 %>%
  dplyr::group_by(model) %>%
  summarize_all(function(x) paste0(round(mean(x),3),"(",round(sd(x),3),")") )

#----- combine results
results3 %>%
  dplyr::select(-gamma) %>%
  dplyr::group_by(model) %>%
  summarize_all(function(x) paste0(round(mean(x),3),"(",round(sd(x),3),")") ) %>%
  mutate(n = 800) %>%
  pivot_longer(c("beta", "M", "theta", "test", "rank"),
               names_to = c("metric")) %>% rbind(
                 results2 %>%
                   dplyr::select(-gamma) %>%
                   dplyr::group_by(model) %>%
                   summarize_all(function(x) paste0(round(mean(x),3),"(",round(sd(x),3),")") ) %>%
                   mutate(n = 600) %>%
                   pivot_longer(c("beta", "M", "theta", "test", "rank"),
                                names_to = c("metric"))
               ) -> tab

require(kableExtra)

metric_labels <- c(
  beta = "RMSE($\\beta$)",
  M    = "RMSE($M$)",
  theta    = "RMSE($\\Theta$)",
  test = "Test error",
  rank       = "Rank($\\Theta$)"
)

# (Optional) method display order (put yours here to match the paper)
method_order <- c(
  "IMR", "SoftImpute", "MCCI",  "Naive"
)

wide <- tab %>%
  mutate(metric = recode(metric, !!!metric_labels)) %>%
  mutate(model = factor(model, levels = method_order)) %>%
  arrange(n, model) %>%
  select(n, model, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

# How many rows per panel (for pack_rows)
panel_counts <- wide %>% count(n, name = "rows")
panel_index  <- stats::setNames(panel_counts$rows,
                                paste0("$n = m = ", panel_counts$n, "$"))

# Build table
kbl(
  wide %>% select(-n),
  booktabs  = TRUE,
  escape    = FALSE,  # keep LaTeX math in headers
 # col.names = c("Method", setdiff(names(wide), c("n", "method"))),
  align     = c("l", rep("r", ncol(wide) - 2)),
  caption   = paste("<span style='color:#000'>",
                    "RMSE, test error, estimated ranks, and their standard errors (in parenthes",
  "under model $\\Theta=X\\beta+M$ and with rank(M)=10, number of covariates = 20,",
  "true rank of $\\Theta=30$, 20% observation rate of $\\Theta$, (n,m)=(600,600),(800,800),",
  "and 50% of $\\beta$ is 0, at random.</span>")
) |>
  pack_rows(index = panel_index) |>
  kable_styling(latex_options = c("hold_position", "striped"), font_size = 12) |>
  row_spec(0, color = "#000") |>
  row_spec(1:nrow(wide), color = "#000")



#--------------------------------

fit <- fit.nlrr$fit
sim1_res(dat, fit2)

fit2 <- IMR::imr.fit(dat$fit_data$Y_full, dat$fit_data$X$Q,
                     lambda_M = fit.nlrr$lambda_M, r = fit.nlrr$rank_M,
                     lambda_beta = fit.nlrr$lambda_beta, warm_start = fit,
                     trace = T)



fit.nlrr$v <- fitsi$fit$v


fit0 <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q,
                                 lambda_beta=0.4, intercept_row = T,
                                 intercept_col = T)
fit0$gamma0


partial_crossprod(matrix(1,nrow(dat$Y),1),
                  matrix(fit0$gamma0, nrow=1),
                  dat$fit_data$valid@i, dat$fit_data$valid@p)
quick_camc_simu_res(dat, fit32$fit)

ttrain <- tvalid <- fit.nlrr$fit$resid
ttrain[dat$fit_data$train==0] = 0
tvalid[dat$fit_data$valid==0] = 0
fitc <- simpute.cv(ttrain, tvalid, fit.nlrr$fit$resid, trace=T, tol=15)
fitc$fit$beta <- fit.nlrr$fit$beta
fitc$fit$beta0 <- fit.nlrr$fit$beta0
sim1_res(dat, fitc$fit)


quick_camc_simu_res(dat, fitsi$fit)

fitmm <- IMR::imr.cv_M(ttrain, tvalid, Y_full = fit.nlrr$fit$resid, trace = T)
fitmm$fit$beta <- fit.nlrr$fit$beta
sim1_res(dat, fitmm$fit)

sim1_res(dat, fit32$fit)
sim1_res(dat, o)



quick_camc_simu_res(dat, fitmmci$fit, T)
