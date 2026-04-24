test_that("Test inversion scenarios", {
  dat <- generate_simulated_data(500, 400, 2, 1, 1, shared = TRUE, snr=0)
  # send an incomplete matrix
  expect_all_true(is.matrix(inv(as_incomplete(dat$Y))))
  #

})
