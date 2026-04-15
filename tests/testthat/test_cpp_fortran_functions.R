library(testthat)


test_that("C++/Fortran functions", {
  library(Matrix)
  library(softImpute)
  # -------------------------------------------------------------------
  #  Data Preparations
  # -------------------------------------------------------------------
  #source("./tests/testthat/helper_generate_simulations.R")
  dat <- generate_simulated_data(100, 200, 5, 5, 5,sparsity_beta = 0)
  d1 <-  matrix(rbinom(10000,100,.2), 100, 100);d1 <- (d1 + t(d1)) / 2
  d2 <-  matrix(rbinom(20000,200,.2), 200, 200);d2 <- (d2 + t(d2)) / 2
  S <- imr_similarity(d1, invert = T, jitter = 1);S
  S2 <- imr_similarity(d2, invert = T, jitter = 1);S
  data <- imr_data(dat$y, dat$x, dat$z, S, S2, 0.2);data

  # -------------------------------------------------------------------
  # Test 1: Row Addition
  # -------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  add_per_row <- rnorm(nrow(mat))
  x_row_test <- mat@x + 0

  add_to_rows_inplace_cpp(x_row_test, mat@i, add_per_row)

  for(i in 1:ncol(mat)){
    mat[,i] <- mat[,i] + add_per_row
  }
  mat[dat$y == 0] = NA
  mat <- as.Incomplete(mat)
  expect_equal(x_row_test, mat@x,
               info = "Row additions failed to update correct indices.")

  # -------------------------------------------------------------------
  # Test 2: Column Addition
  # -------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  add_per_col <- rnorm(ncol(mat))
  x_row_test <- mat@x + 0 # this is very important

  add_to_cols_inplace_cpp(x_row_test, mat@p, add_per_col)

  for(i in 1:nrow(mat)){
    mat[i,] <- mat[i,] + add_per_col
  }
  mat[dat$y == 0] = NA
  mat <- as.Incomplete(mat)
  expect_equal(x_row_test,  mat@x,
               info = "Column additions failed to update correct indices.")
  #-----------------------------------------------------------------------
  # Test 3,4: row and column means
  #-----------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  nr <- nrow(mat)
  nc <- ncol(mat)

  expected_row_means <- rowMeans(dat$y)
  expected_col_means <- colMeans(dat$y)

  actual_row_means <- row_means_cpp(mat@x, mat@i, nr, nc)
  actual_col_means <- col_means_cpp(mat@x, mat@p, nr, nc)

  expect_equal(actual_row_means, expected_row_means,
               info = "row means do not match base R.")

  expect_equal(actual_col_means, expected_col_means,
               info = "col means do not match base R.")
  #-----------------------------------------------------------------------
  # Test 5,6: Soft-Thresholding
  #-----------------------------------------------------------------------
  lambda <- 0.2
  expected_beta <- matrix(
    sign(dat$beta) * pmax(abs(dat$beta) - lambda, 0),
    nrow = nrow(dat$beta),
    ncol = ncol(dat$beta)
  )
  actual_beta <- soft_threshold_cpp(dat$beta, lambda)
  expect_equal(actual_beta, expected_beta,
               info = "Fast soft-thresholding failed to match base R.")
  zeros_expected <- sum(expected_beta == 0)
  zeros_actual   <- sum(actual_beta == 0)
  expect_equal(zeros_actual, zeros_expected,
               info = "C++ function did not produce the correct number of exact zeros.")
  #-----------------------------------------------------------------------
  # Test 7: Frob ratio
  #-----------------------------------------------------------------------
  s1 <- svd(dat$y)
  s2 <- svd(dat$y + matrix(rnorm(2000),100,200))
  U1 <- s1$u; V1 <- s1$v; d1 <- s1$d
  U2 <- s2$u; V2 <- s2$v; d2 <- s2$d

  expected <- softImpute:::Frob(U1, d1, V1, U2, d2, V2)
  actual <- frob_ratio_cpp(U1, d1, V1, U2, d2, V2)
  expect_equal(actual, expected,
               info = "Frob ratio does not match.")
  #-----------------------------------------------------------------------
  # Test 8: Updating A without similarity
  #-----------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  lambda <- 1.2
  expected <-(mat %*% V1 + U1 %*% diag(d1,100,100)) %*%
    diag(1 / (1 + lambda / d1),100,100) %>%
    as.matrix()
  actual <- update_A_cpp(mat, U1, V1, d1, 1.2)
  expect_equal(actual, expected,
               info = "update_A does not match.")
  #-----------------------------------------------------------------------
  # Test 9: Updating B without similarity
  #-----------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  lambda <- 1.2
  expected <-( t(mat) %*% U1 + V1 %*% diag(d1,100,100)) %*%
    diag(1 / (1 + lambda / d1),100,100) %>%
    as.matrix()
  actual <- update_B_cpp(mat, U1, V1, d1, 1.2)
  expect_equal(actual, expected,
               info = "update_B does not match.")
  #-----------------------------------------------------------------------
  # Test 10: Updating A with similarity
  #-----------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  W = mat %*% V1 %*% diag(d1) + U1 %*% diag(d1^2)
  expected = matrix(NA, nrow(U1), ncol(U1))
  for(j in 1:ncol(expected))
    expected[,j] <- S$U %*% (diag(1,nrow(expected)) *
                               (1/(S$d+d1[j]))) %*%
    t(S$U) %*% W[,j]

  actual <- update_A_sim_cpp(mat, U1, V1, d1, S$U, S$d)
  expect_equal(actual, expected,
               info = "update_A_sim does not match.")
  #-----------------------------------------------------------------------
  # Test 11: Updating B with similarity
  #-----------------------------------------------------------------------
  mat <- as.Incomplete(dat$y)
  W = t(diag(d1) %*% t(U1) %*% mat + diag(d1^2) %*% t(V1))
  expected = matrix(NA, nrow(V1), ncol(V1))
  for(j in 1:ncol(expected))
    expected[,j] <- S2$U %*% (diag(1,nrow(expected)) *
                                (1/(S2$d+d1[j]))) %*% t(S2$U) %*% W[,j]

  actual <- update_B_sim_cpp(mat, U1, V1, d1, S2$U, S2$d)
  expect_equal(actual, expected,
               info = "update_A_sim does not match.")
  #-----------------------------------------------------------------------
  # Test 12,13: SVD small nr
  #-----------------------------------------------------------------------
  unsvd <- function(x) x$u %*% (t(x$v) * x$d)
  expected <- svd(t(dat$x))
  actual <- svd_small_nr_cpp(t(dat$x))
  expect_equal(actual$d, expected$d,
               info = "SVD small nr eigenvalues don't match")
  expect_equal(unsvd(actual), unsvd(expected),
               info = "SVD small nr reconstructions don't match")
  #-----------------------------------------------------------------------
  # Test 14,15: SVD small nc
  #-----------------------------------------------------------------------
  expected <- svd((dat$x))
  actual <- svd_small_nc_cpp((dat$x))
  expect_equal(actual$d, expected$d,
               info = "SVD small nr eigenvalues don't match")
  expect_equal(unsvd(actual), unsvd(expected),
               info = "SVD small nr reconstructions don't match")
  #-----------------------------------------------------------------------
  # Test 16: [Fortran]  crossprod.f90
  #-----------------------------------------------------------------------
  actual <- (dat$x %*% dat$beta)[as.matrix(data$valid_mask)==1]
  expected <- partial_crossprod(dat$x, dat$beta, data$valid_mask@i, data$valid_mask@p)
  expect_equal(actual, expected,
               info = "Crossprod don't match")
  #-----------------------------------------------------------------------
  # Test 17: [Fortran]  crossprodt.f90
  #-----------------------------------------------------------------------
  actual <- (dat$gamma %*% t(dat$z))[as.matrix(data$valid_mask)==1]
  expected <- partial_crossprod(dat$gamma, dat$z, data$valid_mask@i, data$valid_mask@p,T)
  expect_equal(actual, expected,
               info = "Crossprod don't match")


})
