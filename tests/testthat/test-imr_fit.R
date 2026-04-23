testthat::test_that("imr_fit function", {
  dat <- generate_simulated_data(500, 400, 2, 1, 1, snr = 0)
  # create data object
  data <- imr_data(dat$Y, dat$X, dat$Z,val_prop =  0.3);data
  # fit model
  fit <- imr_fit(data, 2)
  # reconstruct response matrix
  rec <- reconstruct(fit, data, trace = FALSE)
  # get spearman correlation between estimates and original matrix. should be 1.
    testthat::expect_equal(round(evaluate(rec$estimates, dat$theta, "spearman")),
                           1,
                           info = "simple fit function test didn't result in 100% correlation")

  # test with warmstart
  fit <- imr_fit(data, 2, warm_start = fit)
  # reconstruct response matrix
  rec <- reconstruct(fit, data, trace = FALSE)
  # get spearman correlation between estimates and original matrix. should be 1.
    testthat::expect_equal(round(evaluate(rec$estimates, dat$theta, "spearman")),
                           1)

})

testthat::test_that("imr_fit function with shared covariates",{

  dat <- generate_simulated_data(500, 400, 2, 1, 1, shared = TRUE, snr=0)
  data <- imr_data(dat$Y, dat$X, dat$Z,val_prop =  0.3);
  data <- update(data, row_covariates = TRUE, col_covariates = TRUE,low_rank_component = TRUE,
                 shared_beta = TRUE, shared_gamma = TRUE,row_similarity = FALSE,
                 col_similarity = FALSE, row_intercept = FALSE, col_intercept = FALSE)
  # fit model
  fit <- imr_fit(data, 2)
  # reconstruct response matrix
  rec <- reconstruct(fit, data, trace = FALSE)
  expect_equal(rec$beta, dat$beta, tolerance = 0.01)
  expect_equal(rec$gamma, dat$gamma, tolerance = 0.01)

  # we now see that if we can estimate xbeta and gammaz using the intercepts only
  data <- update(data, row_covariates = FALSE, col_covariates = FALSE,low_rank_component = TRUE,
                 shared_beta = FALSE, shared_gamma = FALSE,row_similarity = FALSE,
                 col_similarity = FALSE, row_intercept = TRUE, col_intercept = TRUE)
  # fit model
  fit <- imr_fit(data, 2)
  silent(print(fit))
  rec <- reconstruct(fit, data, trace = FALSE)
  corrcoef1 <- evaluate(rec$beta0, dat$X %*% dat$beta, "spearman")
  corrcoef2 <- evaluate(rec$gamma0, as.vector(dat$gamma %*% t(dat$Z)), "spearman")
  expect_gte(min(corrcoef1, corrcoef2), 0.95)
  # see if you get the same results with warm_start on
  # fit model
  fit <- imr_fit(data, 2, warm_start = fit, training = TRUE)
  rec <- reconstruct(fit, data, trace = FALSE)
  corrcoef1 <- evaluate(rec$beta0, dat$X %*% dat$beta, "spearman")
  corrcoef2 <- evaluate(rec$gamma0, as.vector(dat$gamma %*% t(dat$Z)), "spearman")
  expect_gte(min(corrcoef1, corrcoef2), 0.95)

  #  we check trace=TRUE, and random initializations are working
  # this should give a warning of no convergence due to small number of iteratins
  fit <- silent(suppressWarnings(imr_fit(data, 2, convergence = imr_convergence(5, 0.001, TRUE, FALSE))))


})

testthat::test_that("mc_with_means", {
  Y_input <- matrix(c(
    1, NA,  3,
    4,  5, NA,
    0,  8,  9
  ), nrow = 3, byrow = TRUE)

  Y_expected <- matrix(c(
    1.0, 4.25, 3.00,
    4.0, 5.00, 5.25,
    5.5, 8.00, 9.00
  ), nrow = 3, byrow = TRUE)

  Y_out <- as.matrix(mc_with_means(Y_input))

  expect_equal(Y_out, Y_expected)

})

