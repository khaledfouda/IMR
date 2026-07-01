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
