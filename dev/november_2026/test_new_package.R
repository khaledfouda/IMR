source("./article/imr-reproducible-code/code/simulation/helper.R")


# devtools::build()
devtools::document()
devtools::load_all()

tp <- function(p, t, d = 5) mean(round(p, 5) == 0 & t == 0)
fp <- function(p, t, d = 5) mean(round(p, 5) == 0 & t != 0)


convergence <- imr_convergence()

grid <- imr_tune_grid(rank = 4,
  metric = "rmse")

metrics <- function(out, simd, metric = "rmse", digits=4) c(
  theta = evaluate(out$estimates, simd$theta, metric = metric),
  beta = evaluate(out$beta, simd$beta, metric = metric),
  gamma = evaluate(out$gamma, simd$gamma, metric = metric),
  beta0 = evaluate(out$beta0, simd$beta0, metric = metric),
  gamma0 = evaluate(out$gamma0, simd$gamma0, metric = metric)
  #beta_nonzero = evaluate(out$beta, as_incomplete(simd$beta), metric = metric),
  #gamma_nonzero = evaluate(out$gamma, as_incomplete(simd$gamma), metric = metric)
  #beta_tp = tp(out$beta, simd$beta),
  #beta_fp = fp(out$beta, simd$beta),
  #gamma_tp = tp(out$gamma, simd$gamma),
  #gamma_fp = fp(out$gamma, simd$gamma
  )  %>% round(digits)


# one scenario: simulate, fit, identify, and return the two rows of metrics
run_scenario <- function(intercepts, center, scale) {
  simd <- generate_simulated_data(600, 800, r = 4, p = 3, q = 3, sparsity = 0.8,
    sparsity_beta = 0, sparsity_gamma = 0,
    snr = 3, mechanism = "mcar",
    identified = TRUE,
    row_intercept = intercepts,
    col_intercept = intercepts,
    seed = 2025)

  fit1 <- imr(simd$Y, simd$X, simd$Z, row_intercept = intercepts,
    col_intercept = intercepts, val_prop = 0.2,  center = center, scale = scale,
    seed = 2025, grid = grid, fast_nuclear = TRUE,
    maxit = 2000, thresh = 1e-7, n_cores = 8, verbose = 1)

  out1 <- reconstruct(fit1$fit, fit1$data)

  fitid1 <- imr_identify(fit1$fit, fit1$data)
  out2 <- reconstruct(fitid1, fit1$data)

  data.frame(
    scenario = sprintf("intercepts=%s, center=%s, scale=%s",
      intercepts, center, scale),
    model = c("imr", "imr_identified"),
    rbind(metrics(out1, simd), metrics(out2, simd))
  )
}


scenarios <- data.frame(
  intercepts = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  center     = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
  scale      = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
)

results <- do.call(rbind, Map(run_scenario,
  scenarios$intercepts, scenarios$center, scenarios$scale))
rownames(results) <- NULL
results
