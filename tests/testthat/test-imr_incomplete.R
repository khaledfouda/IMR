test_that("as_incomplete() drops NAs and 0s to create a sparse matrix", {
  # Create sample data
  y_input <- matrix(
    c(2,  NA,  0,
      4, 0.5, NA,
      NA,  NA,  0),
    nrow = 3,
    byrow = TRUE
  )

  # Make it sparse
  y_sparse <- as_incomplete(y_input)

  # Verify it was successfully converted to an incomplete object
  expect_true(is_incomplete(y_sparse))

  # Verify exactly 3 valid elements remain (2, 4, 0.5)
  expect_length(y_sparse@x, 3)
})
