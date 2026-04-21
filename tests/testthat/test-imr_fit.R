testthat::test_that("imr_fit function", {
  dat <- generate_simulated_data(500, 400, 2, 1, 1)
  # remove noise for exact recovery
  dat$Y <- dat$theta * dat$mask
  # create data object
  data <- imr_data(dat$Y, dat$X, dat$Z,val_prop =  0.3);data
  # fit model
  fit <- imr_fit(data, 2)
  # reconstruct response matrix
  rec <- reconstruct(fit, data)
  # get spearman correlation between estimates and original matrix. should be 1.
  estimated =
    testthat::expect_equal(round(evaluate(rec$estimates, dat$theta, "spearman")),
                           1,
                           info = "simple fit function test didn't result in 100% correlation")


})
