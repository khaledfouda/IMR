test_that("imr_tune forwards nuclear_log_scale to the lasso scenarios", {
  dat <- generate_simulated_data(30, 40, 2, 1, 1, snr = 1, seed = 5)
  # streaks larger than the grid length, so every lambda_m value is fitted
  grid <- imr_tune_grid(beta = c(0, 1, 2), gamma = c(0, 1, 2),
                        nuclear = c(0.01, 1, 4, 10), rank = c(2, 4, 2, 10))
  lin <- seq(1, 0.01, length.out = 4)
  logs <- exp(seq(log(1), log(0.01), length.out = 4))
  lambda_values <- function(out) sort(unique(out$history$lambda_m), decreasing = TRUE)

  # Scenario 2: only lambda_beta is tuned
  data_x <- imr_data(dat$Y, dat$X, val_prop = 0.2, seed = 1)
  out_lin <- suppressMessages(imr_tune(data_x, grid, final_fit = FALSE, n_cores = 1,
                                       verbose = 0, nuclear_log_scale = FALSE))
  out_log <- suppressMessages(imr_tune(data_x, grid, final_fit = FALSE, n_cores = 1,
                                       verbose = 0, nuclear_log_scale = TRUE))
  expect_equal(lambda_values(out_lin), lin)
  expect_equal(lambda_values(out_log), logs)

  # Scenario 3: lambda_beta and lambda_gamma are tuned in alternation
  data_xz <- imr_data(dat$Y, dat$X, dat$Z, val_prop = 0.2, seed = 1)
  out_lin <- suppressMessages(imr_tune(data_xz, grid, final_fit = FALSE, n_cores = 1,
                                       verbose = 0, nuclear_log_scale = FALSE, tune_maxit = 1))
  expect_equal(lambda_values(out_lin), lin)
})
