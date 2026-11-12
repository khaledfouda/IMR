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

test_that("get_metric('medae') is the median absolute error", {
  predicted <- c(1, 4, NA, 10, 3)
  actual <- c(2, 2, 5, 3, NA)
  # the complete pairs have absolute errors 1, 2 and 7: no rounding, no tolerance
  expect_identical(get_metric("medae")(predicted, actual), 2)
  expect_identical(evaluate(predicted, actual, metric = "medae"), 2)
  expect_identical(get_metric("medae")(predicted, actual, na.rm = FALSE), NA_real_)
})

# Data for the evaluate() tests: a dense truth, predictions, the observed set
set.seed(11)
ev_theta <- matrix(rnorm(20 * 15), 20, 15)
ev_pred <- ev_theta + matrix(rnorm(20 * 15, sd = 0.3), 20, 15)
ev_obs <- matrix(rbinom(20 * 15, 1, 0.5), 20, 15) == 1
ev_y <- as_incomplete(ev_theta * ev_obs)
ev_mask <- matrix(rbinom(20 * 15, 1, 0.4), 20, 15)

test_that("evaluate() uses only the observed entries of an incomplete actual", {
  # reference: the same metrics on the dense entries selected by hand; the
  # selected values are the same doubles, hence no tolerance
  expect_identical(evaluate(ev_pred, ev_y),
                   evaluate(ev_pred[ev_obs], ev_theta[ev_obs]))
  # the zero-filled comparison of the previous code gives another value
  expect_false(isTRUE(all.equal(
    evaluate(ev_pred, ev_y, "rmse"),
    sqrt(mean((ev_pred - ev_theta * ev_obs)^2))
  )))
})

test_that("evaluate() takes a vector aligned with the observed entries", {
  expect_identical(evaluate(ev_pred[ev_obs], ev_y),
                   evaluate(ev_pred[ev_obs], ev_theta[ev_obs]))
})

test_that("evaluate() uses only the entries where mask == 1", {
  sel <- ev_mask == 1
  ref <- evaluate(ev_pred[sel], ev_theta[sel])
  # a dense mask, the same mask stored sparse, and an aligned vector
  expect_identical(evaluate(ev_pred, ev_theta, mask = ev_mask), ref)
  expect_identical(evaluate(ev_pred, ev_theta, mask = as_incomplete(ev_mask)),
                   ref)
  expect_identical(evaluate(ev_pred[sel], ev_theta, mask = ev_mask), ref)
})

test_that("evaluate() intersects the mask with an incomplete actual", {
  # the mask also selects entries that ev_y does not store
  expect_true(any(ev_mask == 1 & !ev_obs))
  sel <- ev_mask == 1 & ev_obs
  expect_identical(evaluate(ev_pred, ev_y, mask = ev_mask),
                   evaluate(ev_pred[sel], ev_theta[sel]))
})

test_that("evaluate() rejects a mask of other dimensions", {
  expect_error(evaluate(ev_pred, ev_theta, mask = t(ev_mask)))
})

test_that("imr_tune_grid() accepts the median absolute error", {
  expect_identical(imr_tune_grid(metric = "medae")$metric, "MEDAE")
})
