test_that("Test the creation of a sparse matrix", {
  # create sample data
  Y <- matrix(
    c(2, NA, 0,
      4, .5, NA,
      NA, NA, 0), 3, byrow= TRUE
  )
  # make it sparse with both NAs and 0s dropped
  Y <- as_incomplete(Y)
  # verify
  expect_equal(TRUE, is_incomplete(Y))
  expect_equal(length(Y@x), 3)
})
