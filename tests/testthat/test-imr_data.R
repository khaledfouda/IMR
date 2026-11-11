sq_mat <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
non_sq_mat <- matrix(1:6, nrow = 2, ncol = 3)
d_mat <- as.matrix(dist(c(1, 2)))

test_that("imr_data keeps the 1.0.0 positional order: val_prop and seed before the masks", {
  expect_equal(names(formals(imr_data))[6:9], c("val_prop", "seed", "train_mask", "valid_mask"))

  dat <- generate_simulated_data(20, 30, 2, 1, 1, snr = 0, seed = 1)
  data <- imr_data(dat$Y, dat$X, dat$Z, NULL, NULL, 0.2, 1)
  expect_true(data$meta$split_data)
  expect_false(is.null(data$valid_mask))
  expect_equal(data$meta$val_prop, 0.2)
})

test_that("imr_data uses train_mask and valid_mask passed by name as the split", {
  dat <- generate_simulated_data(10, 12, 2, 0, 0, sparsity = 0.5, snr = 0, seed = 3)
  obs <- dat$Y != 0
  set.seed(3)
  valid <- matrix(0, 10, 12)
  valid[sample(which(obs), 15)] <- 1
  train <- obs * 1 - valid

  data <- imr_data(dat$Y, train_mask = train, valid_mask = valid)
  expect_true(data$meta$split_data)
  expect_equal(data$meta$n_train, sum(train))
  expect_equal(data$meta$n_valid, 15)
  # tolerance = 0: masking multiplies by exact 0/1, so the values must be bit-identical
  expect_equal(as.matrix(data$y_valid), dat$Y * valid, tolerance = 0)
  expect_equal(Matrix::which(data$y_valid != 0), which(valid == 1))
  expect_equal(as.matrix(data$y_train), dat$Y * train, tolerance = 0)
})

test_that("imr_similarity handles matrix inputs and invalid classes correctly", {
  # Valid square matrix
  res_mat <- imr_similarity(sq_mat)
  expect_output(print(res_mat))
  expect_s3_class(res_mat, "imr_similarity")
  expect_equal(res_mat$meta$source, "User Matrix (Normalized) (Inverted)")

  expect_error(imr_similarity(non_sq_mat), "must be square")

  # Input that is neither matrix nor character
  expect_error(imr_similarity(1:5), "must be a matrix or a character string")
})

test_that("imr_similarity computes Matern and RBF kernels and handles warnings", {
  expect_error(imr_similarity("matern"), "Distance matrix 'd' is required")
  expect_error(imr_similarity("gaussian", d = d_mat), "Unknown method")
expect_warning(
    res_matern <- imr_similarity("matern", d = d_mat, normalize = TRUE, invert = FALSE)
  )
  expect_equal(res_matern$meta$source, "Matern Kernel (Normalized)")
  expect_equal(res_matern$meta$params[["smoothness"]], 1.5)
  expect_output(print(res_matern))
  expect_warning(
    res_rbf <- imr_similarity("RBF", d = d_mat, rbf_ell = 2, invert = FALSE, normalize = FALSE)
  )
  expect_equal(res_rbf$meta$source, "RBF Kernel")
  expect_equal(res_rbf$meta$params[["ell"]], 2)
})

test_that("imr_similarity correctly applies jitter and invert parameters", {
  res_jitter <- imr_similarity(sq_mat, jitter = 0.5)
  expect_equal(res_jitter$meta$jitter, 0.5)
  expect_output(print(res_jitter))

  expect_warning(res_no_jitter <- imr_similarity(sq_mat, jitter = 0))
  expect_equal(res_no_jitter$meta$jitter, 0)
  expect_no_warning(
    res_inv <- imr_similarity("rbf", d = d_mat, invert = TRUE)
  )
  expect_true(res_inv$meta$inverted)
  expect_match(res_inv$meta$source, "\\(Inverted\\)")
})
