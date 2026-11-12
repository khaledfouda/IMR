# Fits without a low-rank component that start at zero used to stop after one
# iteration: the stopping rule skipped every block whose previous norm was zero,
# and it ignored the intercepts.

# least-squares fitted values at the observed entries, for the design built from
# row and column intercepts and the per-column (X) and per-row (Z) covariate effects
ls_fitted <- function(Y, X = NULL, Z = NULL) {
  Y <- as_incomplete(Y)
  i <- Y@i + 1
  j <- rep(seq_len(ncol(Y)), diff(Y@p))
  rows <- outer(i, seq_len(nrow(Y)), "==") * 1
  cols <- outer(j, seq_len(ncol(Y)), "==") * 1
  design <- cbind(rows, cols)
  if (!is.null(X)) for (k in seq_len(ncol(X))) design <- cbind(design, X[i, k] * cols)
  if (!is.null(Z)) for (k in seq_len(ncol(Z))) design <- cbind(design, Z[j, k] * rows)
  qr.fitted(qr(design), Y@x)
}

test_that("a fit without M that starts at zero iterates to the least-squares solution", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  data <- update(imr_data(sim$Y, sim$X, sim$Z), low_rank_component = FALSE,
                 row_intercept = TRUE, col_intercept = TRUE)

  # default convergence settings: more than the single iteration of the old rule
  fit <- suppressMessages(imr_fit(data))
  expect_gt(fit$meta$n_iter, 1)

  fit <- suppressMessages(imr_fit(data, convergence = imr_convergence(maxit = 1e5, thresh = 1e-14)))
  expect_true(fit$meta$converged)
  fitted <- reconstruct_partial(fit, data, data$Y@i, data$Y@p)
  # tolerance 1e-6 (relative): the fit stops at a relative change of 1e-14;
  # the one-iteration fit of the old rule is off by 21% here
  expect_equal(fitted, ls_fitted(sim$Y, sim$X, sim$Z), tolerance = 1e-6)
})

test_that("an intercepts-only fit iterates to the two-way least-squares solution", {
  sim <- generate_simulated_data(15, 20, 2, 0, 0, sparsity = 0.4, snr = 1, seed = 22)
  data <- update(imr_data(sim$Y), low_rank_component = FALSE,
                 row_intercept = TRUE, col_intercept = TRUE)

  fit <- suppressMessages(imr_fit(data, convergence = imr_convergence(maxit = 1e5, thresh = 1e-14)))
  expect_gt(fit$meta$n_iter, 1)
  expect_true(fit$meta$converged)
  fitted <- reconstruct_partial(fit, data, data$Y@i, data$Y@p)
  expect_equal(fitted, ls_fitted(sim$Y), tolerance = 1e-6)
})

test_that("the lambda_m grid maximum is the operator norm of the converged null residual", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 23)
  data <- update(imr_data(sim$Y, sim$X, sim$Z, val_prop = 0.2, seed = 1),
                 row_intercept = TRUE, col_intercept = TRUE)
  tight <- imr_convergence(maxit = 1e5, thresh = 1e-14, ls_initial = FALSE)

  grid <- imr_set_grid_limits(data, imr_tune_grid(beta = 0, gamma = 0, nuclear = c(0, NA, 5, 2)),
                              convergence = tight, verify_iter = 0)
  # independent reference: residual of the least-squares null model on the training entries
  null_resid <- data$y_train
  null_resid@x <- null_resid@x - ls_fitted(data$y_train, sim$X, sim$Z)
  expected <- svd(as.matrix(null_resid))$d[1]
  expect_equal(grid$nuclear$max, expected, tolerance = 1e-6)
})

# For lambda_m at or above sigma_1 of the residual of the model without M, M = 0
# solves the problem (KKT condition of the nuclear norm). The iteration only
# approaches it: a cold fit there used to run about thresh^(-1/2) iterations, or
# to reach maxit, and to return a small nonzero d.

# sigma_1 of the residual of the least-squares model with both intercepts and no M
critical_lambda <- function(Y, X = NULL, Z = NULL) {
  resid <- as_incomplete(Y)
  resid@x <- resid@x - ls_fitted(Y, X, Z)
  svd(as.matrix(resid))$d[1]
}

test_that("a cold fit above the critical lambda_m returns the model without M at once", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  critical <- critical_lambda(sim$Y, sim$X, sim$Z)

  expect_no_message(
    fit <- imr_fit(data, rank = 2, lambda_m = 1.01 * critical,
                   convergence = imr_convergence(maxit = 1e5, thresh = 1e-14))
  )
  expect_equal(fit$meta$n_iter, 0)
  expect_true(all(fit$coefficients$d == 0))
  fitted <- reconstruct_partial(fit, data, data$Y@i, data$Y@p)
  # tolerance 1e-6 (relative): the fit of the model without M stops at a relative
  # change of 1e-14, as in the tests above
  expect_equal(fitted, ls_fitted(sim$Y, sim$X, sim$Z), tolerance = 1e-6)
})

test_that("a cold fit below the critical lambda_m still iterates to a nonzero M", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  critical <- critical_lambda(sim$Y, sim$X, sim$Z)

  fit <- suppressMessages(imr_fit(data, rank = 2, lambda_m = 0.9 * critical,
                                  convergence = imr_convergence(maxit = 1e4, thresh = 1e-8)))
  expect_gt(fit$meta$n_iter, 0)
  expect_gt(fit$coefficients$d[1], 0)
})

test_that("the fit at the lambda_m grid maximum is the model without M", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 23)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-6)
  grid <- imr_tune_grid(beta = 0, gamma = 0, nuclear = c(0, NA, 5, 2))
  # with covariates and intercepts, then with Y alone. The maximum comes from
  # another SVD call than the one of the fit: the two agree only up to rounding
  for (full in c(TRUE, FALSE)) {
    data <- imr_data(sim$Y, if (full) sim$X, if (full) sim$Z, val_prop = 0.2, seed = 1)
    data <- update(data, row_intercept = full, col_intercept = full)
    lambda_max <- imr_set_grid_limits(data, grid, convergence = conv,
                                      verify_iter = 0)$nuclear$max
    expect_no_message(
      fit <- imr_fit(data, rank = 2, lambda_m = lambda_max, convergence = conv, training = TRUE)
    )
    expect_equal(fit$meta$n_iter, 0)
    expect_true(all(fit$coefficients$d == 0))
  }
})

test_that("imr_tune() takes the first point of the lambda_m path from the model without M", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 23)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-6)
  data <- update(imr_data(sim$Y, sim$X, sim$Z, val_prop = 0.2, seed = 1),
                 row_intercept = TRUE, col_intercept = TRUE)
  grid <- imr_tune_grid(beta = 0, gamma = 0, nuclear = c(0, NA, 4, 2), rank = c(2, 4, 2, 2))
  grid <- imr_set_grid_limits(data, grid, convergence = conv, verify_iter = 0)

  # reference: the validation error of the model fitted without M on the same split
  null_data <- update(imr_data(sim$Y, sim$X, sim$Z, val_prop = 0.2, seed = 1),
                      low_rank_component = FALSE, row_intercept = TRUE, col_intercept = TRUE)
  null_fit <- imr_fit(null_data, convergence = conv, training = TRUE)
  valid <- data$y_valid
  expected <- get_metric("rmse")(reconstruct_partial(null_fit, null_data, valid@i, valid@p),
                                 valid@x)

  for (fast in c(TRUE, FALSE)) {
    # the fit that follows a fit with M = 0 starts cold, without a warning
    expect_no_warning(
      out <- suppressMessages(imr_tune(data, grid, final_fit = FALSE, fast_nuclear = fast,
                                       convergence = conv, verbose = 0, n_cores = 1, seed = 1))
    )
    expect_equal(out$history$rank_out[1], 0)
    # the same fit of the model without M on both sides: default tolerance
    expect_equal(out$history$verror[1], expected)
  }
})

test_that("a model with similarity matrices is not screened", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  sr <- imr_similarity("rbf", d = as.matrix(dist(seq_len(nrow(sim$Y)))) / 10)
  sc <- imr_similarity("rbf", d = as.matrix(dist(seq_len(ncol(sim$Y)))) / 10)
  data <- imr_data(sim$Y, sim$X, sim$Z, similarity_rows = sr, similarity_cols = sc)
  # far above sigma_1 of the residual of the model without M: with similarity
  # matrices the penalty is not the nuclear norm, and that bound does not apply
  fit <- suppressMessages(imr_fit(data, rank = 2, lambda_m = 1e3 * critical_lambda(sim$Y),
                                  convergence = imr_convergence(maxit = 5)))
  expect_gt(fit$meta$n_iter, 0)
})

test_that("a fit started from a fit with M = 0 equals the cold fit, without a warning", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  critical <- critical_lambda(sim$Y, sim$X, sim$Z)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-8)

  null <- imr_fit(data, rank = 2, lambda_m = 1.01 * critical, convergence = conv)
  expect_no_warning(
    warm <- suppressMessages(imr_fit(data, rank = 2, lambda_m = 0.5 * critical,
                                     convergence = conv, warm_start = null))
  )
  cold <- suppressMessages(imr_fit(data, rank = 2, lambda_m = 0.5 * critical, convergence = conv))
  # the same computation from the same start: default tolerance
  expect_equal(warm$coefficients, cold$coefficients)
})

test_that("print() and summary() accept a fit with M = 0", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 21)
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  fit <- imr_fit(data, rank = 2, lambda_m = 1.01 * critical_lambda(sim$Y, sim$X, sim$Z),
                 convergence = imr_convergence(maxit = 1e4, thresh = 1e-8))
  expect_true(all(fit$coefficients$d == 0))
  expect_output(print(fit), "Converged")
  expect_output(summary(fit))
})
