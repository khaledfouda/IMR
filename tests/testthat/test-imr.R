# Data and grids for the test file
sim <- generate_simulated_data(25, 30, 2, 1, 1, sparsity = 0.5, snr = 1, seed = 7)
grid_m <- imr_tune_grid(nuclear = c(0, NA, 5, 2), rank = c(2, 6, 2, 2))
grid_all <- imr_tune_grid(beta = c(0, NA, 3), gamma = c(0, NA, 3),
                          nuclear = c(0, NA, 4, 2), rank = c(2, 6, 2, 2))

# The reference for imr(): the calls it replaces, with the same seed and the
# same settings. `...` goes to imr_tune().
manual_imr <- function(Y, X = NULL, Z = NULL, grid, seed, row_intercept = FALSE,
                       col_intercept = FALSE, thresh = 1e-5, ...) {
  convergence <- imr_convergence(maxit = 600, thresh = thresh, ls_initial = TRUE)
  data <- imr_data(Y, X, Z, val_prop = 0.2, seed = seed)
  data <- update(data, row_intercept = row_intercept, col_intercept = col_intercept)
  grid <- imr_set_grid_limits(data, grid, convergence = convergence, training = TRUE)
  out <- imr_tune(data, grid, convergence = convergence, verbose = 0, n_cores = 1,
                  seed = seed, ...)
  out$data <- data
  out$grid <- grid
  out
}

# imr() and its reference run the same computations from the same seed. The
# default tolerance of expect_equal() (1.5e-8) only allows for a BLAS whose
# rounding varies from one call to the next. The run times are not compared.
expect_same_tuning <- function(out, ref) {
  expect_equal(out$params, ref$params)
  expect_equal(out$history, ref$history)
  expect_equal(out$fit$coefficients, ref$fit$coefficients)
  expect_equal(out$grid, ref$grid)
  expect_equal(out$data, ref$data)
}

test_that("imr() equals the manual calls when only the low-rank part is tuned", {
  # the slow path, so that fast_nuclear is seen to reach imr_tune()
  out <- suppressMessages(imr(sim$Y, grid = grid_m, fast_nuclear = FALSE, n_cores = 1,
                              seed = 1, verbose = 0))
  ref <- suppressMessages(manual_imr(sim$Y, grid = grid_m, seed = 1, fast_nuclear = FALSE))
  expect_same_tuning(out, ref)
})

test_that("imr() equals the manual calls when lambda_beta is tuned as well", {
  # a linear lambda_m grid, so that nuclear_log_scale is seen to reach imr_tune()
  out <- suppressMessages(imr(sim$Y, sim$X, grid = grid_all, nuclear_log_scale = FALSE,
                              n_cores = 1, seed = 1, verbose = 0))
  ref <- suppressMessages(manual_imr(sim$Y, sim$X, grid = grid_all, seed = 1,
                                     nuclear_log_scale = FALSE))
  expect_same_tuning(out, ref)
})

test_that("imr() equals the manual calls when both Lasso penalties are tuned", {
  skip_on_cran()
  # both intercepts, and a threshold other than the default
  out <- suppressMessages(imr(sim$Y, sim$X, sim$Z, row_intercept = TRUE, col_intercept = TRUE,
                              grid = grid_all, thresh = 1e-4, n_cores = 1, seed = 1,
                              verbose = 0))
  ref <- suppressMessages(manual_imr(sim$Y, sim$X, sim$Z, grid = grid_all, seed = 1,
                                     row_intercept = TRUE, col_intercept = TRUE,
                                     thresh = 1e-4))
  expect_same_tuning(out, ref)
})

test_that("imr() returns the data object its fit was computed on", {
  out <- suppressMessages(imr(sim$Y, sim$X, col_intercept = TRUE, grid = grid_all,
                              n_cores = 1, seed = 1, verbose = 0))
  y <- out$data$Y
  fitted <- reconstruct_partial(out$fit, out$data, y@i, y@p)
  # fitted values = Y - residuals at the observed entries. The residual is
  # updated in place at every iteration, so it carries the accumulated rounding
  # (about 1e-15 here); 1e-8 leaves room for longer runs.
  expect_equal(fitted, y@x - out$fit$residuals@x, tolerance = 1e-8)
})

test_that("imr() asks for a validation split when it tunes", {
  expect_error(imr(sim$Y, grid = grid_m, val_prop = 0, verbose = 0), "val_prop")
})

test_that("imr() fits once, without a split, when the grid fixes every parameter", {
  grid <- imr_tune_grid(beta = 0.5, gamma = 0.2, nuclear = 1, rank = 3)
  out <- suppressMessages(imr(sim$Y, sim$X, sim$Z, col_intercept = TRUE, grid = grid,
                              seed = 3, verbose = 0))

  data <- update(imr_data(sim$Y, sim$X, sim$Z, val_prop = 0, seed = 3), col_intercept = TRUE)
  fit <- suppressMessages(imr_fit(data, rank = 3, lambda_m = 1, lambda_beta = 0.5,
                                  lambda_gamma = 0.2,
                                  convergence = imr_convergence(maxit = 600, thresh = 1e-5)))
  expect_named(out, c("fit", "data", "grid"))
  expect_false(out$data$meta$split_data)
  # the same computation on both sides: see expect_same_tuning() for the tolerance
  expect_equal(out$fit$coefficients, fit$coefficients)
})

test_that("imr() ignores the grid entries of covariates that are not given", {
  # beta and gamma keep their default ranges, but there is no X and no Z
  out <- suppressMessages(imr(sim$Y, grid = imr_tune_grid(nuclear = 1, rank = 3), seed = 3,
                              verbose = 0))
  fit <- suppressMessages(imr_fit(imr_data(sim$Y, seed = 3), rank = 3, lambda_m = 1))
  expect_named(out, c("fit", "data", "grid"))
  # the same computation on both sides: see expect_same_tuning() for the tolerance
  expect_equal(out$fit$coefficients, fit$coefficients)
})

test_that("imr() selects on the median absolute error when asked to", {
  # final_fit = FALSE returns the selected fit on the training split
  out <- suppressMessages(imr(sim$Y, grid = grid_m, metric = "medae", final_fit = FALSE,
                              n_cores = 1, seed = 2, verbose = 0))
  valid <- out$data$y_valid
  fitted <- reconstruct_partial(out$fit, out$data, valid@i, valid@p)
  expect_identical(out$metric, "MEDAE")
  # the closed form of the metric at the validation entries, computed with the
  # operations of imr_tune(): only the default tolerance of expect_equal()
  expect_equal(out$params$verror, median(abs(valid@x - fitted)))
})

test_that("imr() builds the convergence settings from maxit and thresh", {
  out <- suppressMessages(imr(sim$Y, grid = imr_tune_grid(nuclear = 1, rank = 3),
                              maxit = 50, thresh = 1e-3, verbose = 0))
  expect_identical(out$fit$convergence, imr_convergence(maxit = 50, thresh = 1e-3))
})

test_that("imr() reports its own steps, and only them, at verbose = 1", {
  # a loose threshold, so that no fit reaches maxit and prints the solver's message
  msgs <- capture_messages(imr(sim$Y, grid = grid_m, thresh = 1e-3, n_cores = 1, seed = 1,
                               verbose = 1))
  # the data, the grid limits and the tuning: one line before and one after each
  expect_length(msgs, 6)
  expect_true(all(startsWith(msgs, "imr: ")))
})

test_that("imr() is silent at verbose = 0", {
  # a loose threshold, so that no fit reaches maxit and prints the solver's message
  expect_silent(imr(sim$Y, grid = grid_m, thresh = 1e-3, n_cores = 1, seed = 1, verbose = 0))
})
