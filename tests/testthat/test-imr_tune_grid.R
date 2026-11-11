test_that("imr_tune_grid() creates expected grid objects", {
  # Test default grid
  grid_default <- imr_tune_grid()
  expect_output(print(grid_default))

  # Test custom grid initialized to zero
  grid_zeros <- imr_tune_grid(beta = 0, gamma = 0, nuclear = 0, rank = 0)
  expect_output(print(grid_zeros))

  #  check if all min/max are 0
  grid_vals <- c(
    grid_zeros$beta$min, grid_zeros$gamma$min, grid_zeros$nuclear$min, grid_zeros$rank$min,
    grid_zeros$beta$max, grid_zeros$gamma$max, grid_zeros$nuclear$max, grid_zeros$rank$max
  )
  expect_equal(grid_vals, rep(0, 8))
})

test_that("imr_set_grid_limits() computes numeric max values without NAs", {
  dat <- generate_simulated_data(23, 30, 2, 1, 1, 0.6, snr = 0, seed = 2)
  test_data <- imr_data(dat$Y, dat$X, dat$Z, val_prop = 0.3)
  grid <- imr_tune_grid()
  convergence <- imr_convergence(thresh = 1)
  updated_grid <- imr_set_grid_limits(test_data, grid, verify_iter = 1, refine_iter = 1,
                                      convergence = convergence)

  # Extract the relevant components
  grid_params <- unlist(updated_grid[c("beta", "gamma", "rank", "nuclear")])
  # Check that they are numeric and contain no NAs
  expect_type(grid_params, "double")
  expect_false(anyNA(grid_params))
})

test_that("the default verification keeps the KKT value of lambda_m instead of doubling it", {
  sim <- generate_simulated_data(15, 20, 2, 1, 1, sparsity = 0.4, snr = 1, seed = 23)
  data <- update(imr_data(sim$Y, sim$X, sim$Z, val_prop = 0.2, seed = 1),
                 row_intercept = TRUE, col_intercept = TRUE)
  grid <- imr_tune_grid(beta = 0, gamma = 0, nuclear = c(0, NA, 5, 2))

  kkt <- suppressMessages(imr_set_grid_limits(data, grid, verify_iter = 0))$nuclear$max
  verified <- suppressMessages(imr_set_grid_limits(data, grid))$nuclear$max
  # a fit at exactly the KKT value stops (325 iterations) with its largest singular
  # value still at 0.0086 (zero tolerance 0.0034), so a check made there doubled the
  # bound; 5% above it, the value is 9e-7
  expect_equal(verified, kkt)
})
