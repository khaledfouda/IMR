sq_mat <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
non_sq_mat <- matrix(1:6, nrow = 2, ncol = 3)
d_mat <- as.matrix(dist(c(1, 2)))

test_that("imr_similarity handles matrix inputs and invalid classes correctly", {
  # Valid square matrix
  res_mat <- imr_similarity(sq_mat)
  expect_output(print(res_mat))
  expect_s3_class(res_mat, "imr_similarity")
  expect_equal(res_mat$meta$source, "User Matrix")

  expect_error(imr_similarity(non_sq_mat), "must be square")

  # Input that is neither matrix nor character
  expect_error(imr_similarity(1:5), "must be a matrix or a character string")
})

test_that("imr_similarity computes Matern and RBF kernels and handles warnings", {
  expect_error(imr_similarity("matern"), "Distance matrix 'd' is required")
  expect_error(imr_similarity("gaussian", d = d_mat), "Unknown method")
expect_warning(
    res_matern <- imr_similarity("matern", d = d_mat),
    "Generated a raw Covariance matrix without inversion"
  )
  expect_equal(res_matern$meta$source, "Matern Kernel")
  expect_equal(res_matern$meta$params[["smoothness"]], 1.5)
  expect_output(print(res_matern))
  expect_warning(
    res_rbf <- imr_similarity("RBF", d = d_mat, rbf_ell = 2),
    "Generated a raw Covariance matrix without inversion"
  )
  expect_equal(res_rbf$meta$source, "RBF Kernel")
  expect_equal(res_rbf$meta$params[["ell"]], 2)
})

test_that("imr_similarity correctly applies jitter and invert parameters", {
  res_jitter <- imr_similarity(sq_mat, jitter = 0.5)
  expect_equal(res_jitter$meta$jitter, 0.5)
  expect_output(print(res_jitter))

  res_no_jitter <- imr_similarity(sq_mat, jitter = -1)
  expect_equal(res_no_jitter$meta$jitter, 0)
  expect_no_warning(
    res_inv <- imr_similarity("rbf", d = d_mat, invert = TRUE)
  )
  expect_true(res_inv$meta$inverted)
  expect_match(res_inv$meta$source, "\\(Inverted\\)")
})
