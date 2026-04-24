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
  dat <- generate_simulated_data(23, 30, 2, 1, 1, 0.95, snr = 0)
  test_data <- imr_data(dat$Y, dat$X, dat$Z, val_prop = 0.3)
  grid <- imr_tune_grid()
  convergence <- imr_convergence(thresh = 1)
  updated_grid <- imr_set_grid_limits(test_data, grid, bisection_iter = 1, convergence = convergence)

  # Extract the relevant components
  grid_params <- unlist(updated_grid[c("beta", "gamma", "rank", "nuclear")])
  # Check that they are numeric and contain no NAs
  expect_type(grid_params, "double")
  expect_false(anyNA(grid_params))
})
