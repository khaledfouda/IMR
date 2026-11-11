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
  # the one-iteration fit of the old rule is off by 28% here
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
