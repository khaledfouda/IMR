# generate data
nr <- 300
nc <- 500
dat <- generate_simulated_data(nr, nc, 2, 1, 1, snr = 0)
dat_shared   <- generate_simulated_data(nr, nc, 2, 1, 1, shared = TRUE, snr = 0)

test_that("imr_fit performs basic fitting and reconstruction", {
  # Create data object
  test_data <- imr_data(dat$Y, dat$X, dat$Z, val_prop = 0.3)

  expect_output(print(test_data))

  # Fit model and reconstruct
  fit <- imr_fit(test_data, 2)
  rec <- reconstruct(fit, test_data, trace = FALSE)

  # Check spearman correlation (should be 1 due to the trivial example)
  expect_equal(
    round(evaluate(rec$estimates, dat$theta, "spearman")),
    1,
    info = "Simple fit function test didn't result in 100% correlation"
  )

  # Test with warmstart
  fit_warm <- imr_fit(test_data, 2, warm_start = fit)
  rec_warm <- reconstruct(fit_warm, test_data, trace = FALSE)
  expect_equal(round(evaluate(rec_warm$estimates, dat$theta, "spearman")), 1)

  # Test that evaluate() returns numeric data
  out <- evaluate(rec_warm$estimates, dat$theta)
  expect_true(all(vapply(out, is.numeric, logical(1))))

  # Test partial reconstruction
  rec_partial <- reconstruct_partial(fit_warm, test_data, test_data$Y@i, test_data$Y@p)
  expect_length(rec_partial, length(test_data$Y@x))

  # Invalid config update should error
  expect_error(update(test_data, row_similarity = TRUE))

  # Test warmstart with a different number of eigenvalues
  # Using expect_no_error ensures it executes completely without crashing
  expect_no_error(imr_fit(test_data, 3, warm_start = fit_warm))
  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("the project option is removed from the fitting and tuning functions", {
  funs <- list(imr_fit, imr_solver, imr_tune, imr_tune_lasso, imr_tune_nuclear_fast,
               imr_tune_nuclear_slow, imr_set_grid_limits, imr_lambda_max)
  for (f in funs) expect_false("project" %in% names(formals(f)))

  small <- generate_simulated_data(20, 30, 2, 1, 1, sparsity = 0.5, snr = 0, seed = 1)
  small_data <- imr_data(small$Y, small$X, small$Z)
  expect_error(imr_fit(small_data, 2, project = TRUE), "unused argument")
  fit <- suppressMessages(imr_fit(small_data, 2, convergence = imr_convergence(maxit = 5)))
  expect_false("orthogonal_components" %in% names(fit$meta))
})

test_that("imr_tune and imr_set_grid_limits reject the removed project argument", {
  # Argument matching fails before any argument is evaluated, so no data are needed
  # (and no random numbers are drawn that would shift the state for later tests).
  expect_error(imr_tune(NULL, NULL, project = TRUE), "unused argument")
  expect_error(imr_set_grid_limits(NULL, NULL, project = TRUE), "unused argument")
})

test_that("imr_fit handles shared covariates and intercepts", {
  test_data <- imr_data(dat_shared$Y, dat_shared$X, dat_shared$Z, val_prop = 0.3)

  # Scenario 1: Shared covariates without intercepts
  data_shared_cov <- update(
    test_data,
    row_covariates = TRUE, col_covariates = TRUE, low_rank_component = TRUE,
    shared_beta = TRUE, shared_gamma = TRUE, row_similarity = FALSE,
    col_similarity = FALSE, row_intercept = FALSE, col_intercept = FALSE
  )

  fit_shared <- imr_fit(data_shared_cov, 2)
  rec_shared <- reconstruct(fit_shared, data_shared_cov, trace = FALSE)

  expect_equal(rec_shared$beta, dat_shared$beta, tolerance = 0.1)
  expect_equal(rec_shared$gamma, dat_shared$gamma, tolerance = 0.1)

  # Scenario 2: Intercepts only [intercept should approximate xbeta]
  data_intercepts <- update(
    test_data,
    row_covariates = FALSE, col_covariates = FALSE, low_rank_component = TRUE,
    shared_beta = FALSE, shared_gamma = FALSE, row_similarity = FALSE,
    col_similarity = FALSE, row_intercept = TRUE, col_intercept = TRUE
  )

  fit_intercepts <- imr_fit(data_intercepts, 2)

  expect_output(print(fit_intercepts))
  expect_output(print(summary(fit_intercepts)))
  expect_no_error(imr_convergence())
  expect_length(coef(fit_intercepts)$beta0, nr)

  rec_intercepts <- reconstruct(fit_intercepts, data_intercepts, trace = FALSE)
  corrcoef_beta <- evaluate(rec_intercepts$beta0, dat_shared$X %*% dat_shared$beta, "spearman")
  corrcoef_gamma <- evaluate(rec_intercepts$gamma0, as.vector(dat_shared$gamma %*% t(dat_shared$Z)), "spearman")
  expect_gte(min(corrcoef_beta, corrcoef_gamma), 0.5)

  # Scenario 3: Intercepts with warm_start - should return similar output
  fit_intercepts_warm <- imr_fit(data_intercepts, 2, warm_start = fit_intercepts, training = TRUE)
  rec_intercepts_warm <- reconstruct(fit_intercepts_warm, data_intercepts, trace = FALSE)

  corrcoef_beta_warm <- evaluate(rec_intercepts_warm$beta0, dat_shared$X %*% dat_shared$beta, "spearman")
  corrcoef_gamma_warm <- evaluate(rec_intercepts_warm$gamma0, as.vector(dat_shared$gamma %*% t(dat_shared$Z)), "spearman")
  expect_gte(min(corrcoef_beta_warm, corrcoef_gamma_warm), 0.5)

  expect_message(
    imr_fit(data_intercepts, 2, convergence = imr_convergence(5, 0.001, FALSE, FALSE))
  )

  # giving less eigenvalues than needed. expect no output
  expect_output(imr_fit(data_intercepts, 5, warm_start = fit_intercepts,
                convergence = imr_convergence(5, 0.001, TRUE, FALSE)))
})

test_that("mc_with_means correctly imputes missing values", {
  y_input <- matrix(c(
    1, NA,  3,
    4,  5, NA,
    0,  8,  9
  ), nrow = 3, byrow = TRUE)

  y_expected <- matrix(c(
    1.0, 4.25, 3.00,
    4.0, 5.00, 5.25,
    5.5, 8.00, 9.00
  ), nrow = 3, byrow = TRUE)

  y_out <- as.matrix(mc_with_means(y_input))
  expect_equal(y_out, y_expected)
})

# The tests below seed the global RNG through generate_simulated_data(seed = ).
# They stay at the end of the file so that they do not change the random state
# seen by the earlier tests, which use the unseeded file-level data.
# They check the solver on a given design, so they turn scaling off: without
# intercepts, data$X and data$Z are then sim$X and sim$Z.

test_that("imr_fit objective never increases with lasso penalties on beta and Gamma", {
  sim <- generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, sparsity_beta = 0.5,
                                 sparsity_gamma = 0.5, snr = 2, seed = 11)
  data <- imr_data(sim$Y, sim$X, sim$Z, scale = FALSE)
  obs <- sim$Y != 0
  lam <- c(m = 1, beta = 0.3, gamma = 0.3)
  # dense brute-force objective of the 20 x 30 problem
  objective <- function(fit) {
    cf <- fit$coefficients
    theta <- sim$X %*% cf$beta + cf$gamma %*% t(sim$Z) + cf$u %*% (cf$d * t(cf$v))
    0.5 * sum((obs * (sim$Y - theta))^2) + lam[["m"]] * sum(cf$d) +
      lam[["beta"]] * sum(abs(cf$beta)) + lam[["gamma"]] * sum(abs(cf$gamma))
  }
  # the iterations are deterministic, so the k-th iterate is the fit stopped at maxit = k
  obj <- vapply(1:20, function(k) {
    fit <- suppressMessages(imr_fit(
      data, rank = 4, lambda_m = lam[["m"]], lambda_beta = lam[["beta"]],
      lambda_gamma = lam[["gamma"]],
      convergence = imr_convergence(maxit = k, thresh = 1e-300, ls_initial = FALSE)
    ))
    objective(fit)
  }, numeric(1))
  # tolerance 1e-10 * obj[1] (obj[1] is about 67): room for rounding in the in-place
  # residual updates only; the observed decrease per iteration is at least 0.08
  expect_true(all(diff(obj) <= 1e-10 * obj[1]))
})

test_that("imr_fit residuals equal the brute-force masked residual after lasso Gamma updates", {
  sim <- generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, sparsity_beta = 0.5,
                                 sparsity_gamma = 0.5, snr = 2, seed = 11)
  data <- imr_data(sim$Y, sim$X, sim$Z, scale = FALSE)
  obs <- sim$Y != 0
  fit <- suppressMessages(imr_fit(
    data, rank = 4, lambda_m = 1, lambda_beta = 0.3, lambda_gamma = 0.3,
    convergence = imr_convergence(maxit = 10, ls_initial = FALSE)
  ))
  cf <- fit$coefficients
  theta <- sim$X %*% cf$beta + cf$gamma %*% t(sim$Z) + cf$u %*% (cf$d * t(cf$v))
  expect_length(fit$residuals@x, sum(obs))
  # tolerance 1e-12 (relative): ten iterations of in-place updates accumulate rounding of
  # about 1e-14; any missed or doubled update is of the order of the residual itself
  expect_equal(as.matrix(fit$residuals), obs * (sim$Y - theta), tolerance = 1e-12)
})

test_that("imr_fit Gamma estimate satisfies the lasso KKT conditions at convergence", {
  sim <- generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, sparsity_beta = 0.5,
                                 sparsity_gamma = 0.5, snr = 2, seed = 11)
  data <- imr_data(sim$Y, sim$X, sim$Z, scale = FALSE)
  lam_gamma <- 0.3
  fit <- suppressMessages(imr_fit(
    data, rank = 4, lambda_m = 1, lambda_beta = 0.3, lambda_gamma = lam_gamma,
    convergence = imr_convergence(maxit = 20000, thresh = 1e-10, ls_initial = FALSE)
  ))
  expect_true(fit$meta$converged)
  cf <- fit$coefficients
  theta <- sim$X %*% cf$beta + cf$gamma %*% t(sim$Z) + cf$u %*% (cf$d * t(cf$v))
  # negative gradient of the loss in Gamma, from the dense brute-force residual
  grad <- ((sim$Y != 0) * (sim$Y - theta)) %*% sim$Z
  active <- cf$gamma != 0
  # both sets must be non-empty, otherwise one of the conditions below is vacuous
  expect_true(any(active) && any(!active))
  # tolerance 1e-3: stopping at a relative change of 1e-10 leaves a stationarity residual
  # of about 6e-5; a penalty on the wrong scale (e.g. lambda * LZ) would leave about 0.3
  expect_lt(max(abs(grad[active] - lam_gamma * sign(cf$gamma[active]))), 1e-3)
  # no tolerance: the observed margin, lambda - max |grad| = 0.02, is far above the
  # convergence error
  expect_true(all(abs(grad[!active]) <= lam_gamma))
})

test_that("with lambda_gamma = 0 the Gamma estimate solves the least-squares normal equations", {
  sim <- generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, sparsity_beta = 0.5,
                                 sparsity_gamma = 0.5, snr = 2, seed = 11)
  data <- imr_data(sim$Y, sim$X, sim$Z, scale = FALSE)
  fit <- suppressMessages(imr_fit(
    data, rank = 4, lambda_m = 1, lambda_beta = 0.3, lambda_gamma = 0,
    convergence = imr_convergence(maxit = 20000, thresh = 1e-10, ls_initial = FALSE)
  ))
  expect_true(fit$meta$converged)
  cf <- fit$coefficients
  # no soft-thresholding: the step reduces to a plain gradient step, so no entry is zeroed
  expect_true(all(cf$gamma != 0))
  theta <- sim$X %*% cf$beta + cf$gamma %*% t(sim$Z) + cf$u %*% (cf$d * t(cf$v))
  grad <- ((sim$Y != 0) * (sim$Y - theta)) %*% sim$Z
  # tolerance 1e-3: the observed stationarity residual at a relative change of 1e-10 is about 2e-5
  expect_lt(max(abs(grad)), 1e-3)
})
