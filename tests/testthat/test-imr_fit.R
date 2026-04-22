testthat::test_that("imr_fit function", {
  dat <- generate_simulated_data(500, 400, 2, 1, 1)
  # remove noise for exact recovery
  dat$Y <- dat$theta * dat$mask
  # create data object
  data <- imr_data(dat$Y, dat$X, dat$Z,val_prop =  0.3);data
  # fit model
  fit <- imr_fit(data, 2)
  # reconstruct response matrix
  rec <- reconstruct(fit, data, trace = FALSE)
  # get spearman correlation between estimates and original matrix. should be 1.
  estimated =
    testthat::expect_equal(round(evaluate(rec$estimates, dat$theta, "spearman")),
                           1,
                           info = "simple fit function test didn't result in 100% correlation")


})

testthat::test_that("mc_with_means", {
  Y_input <- matrix(c(
    1, NA,  3,
    4,  5, NA,
    NA,  8,  9
  ), nrow = 3, byrow = TRUE)

  Y_expected <- matrix(c(
    1.0, 4.25, 3.00,
    4.0, 5.00, 5.25,
    5.5, 8.00, 9.00
  ), nrow = 3, byrow = TRUE)

  Y_out <- mc_with_means(Y_input)

  expect_equal(Y_out, Y_expected)

  # 5. Verify no NAs remain (good sanity check)
  expect_false(any(is.na(Y_out)))
})

