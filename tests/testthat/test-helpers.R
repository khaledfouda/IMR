# Generate data once for the test file
dat <- generate_simulated_data(30, 50, 2, 1, 1, shared = TRUE, snr = 0)
# Pre-compute the incomplete matrix
inc_y <- as_incomplete(dat$Y)

test_that("inv() handles incomplete matrices", {
  expect_true(is.matrix(inv(inc_y)))
})

test_that("svd_opt() works correctly under cthin, rthin, and incomplete scenarios", {
  # cthin
  expect_type(svd_opt(dat$X, 2, 0.1)$d, "double")

  # rthin
  expect_type(svd_opt(t(dat$X), 2, 0.1)$d, "double")

  # Incomplete matrix: return only one eigenvalue
  base_svd_d1 <- svd(inc_y)$d[1]
  out <- svd_opt(inc_y, 5, base_svd_d1 - 1)$d
  expect_type(out, "double")
})

test_that("verify_warm_start() keeps a full-rank random start unchanged", {
  set.seed(6)
  n <- 20
  m <- 30
  r <- 3
  # a random initialization: orthonormal u and v, positive d, zero coefficients
  ws <- list(
    u = qr.Q(qr(matrix(rnorm(n * r), n, r))),
    d = c(3, 2, 1),
    v = qr.Q(qr(matrix(rnorm(m * r), m, r))),
    beta = matrix(0, 3, m),
    gamma = matrix(0, n, 2),
    beta0 = numeric(n),
    gamma0 = numeric(m)
  )
  # a rejected start gives a warning and NULL
  expect_no_warning(out <- verify_warm_start(ws, r))
  # no tolerance: with length(d) equal to the rank, nothing is truncated or padded
  expect_identical(out, ws)
})
