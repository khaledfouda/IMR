# Independent tests of evaluate(): every reference is computed by hand on dense
# base-R objects, from the definitions of the metrics, never through evaluate()
# or get_metric().

# Data ----
# a dense truth without exact zeros, predictions, an observed set and a mask
set.seed(2026)
ei_n <- 12
ei_m <- 9
ei_theta <- matrix(rnorm(ei_n * ei_m, mean = 1), ei_n, ei_m)
ei_pred <- ei_theta + matrix(rnorm(ei_n * ei_m, sd = 0.4), ei_n, ei_m)
ei_obs <- matrix(rbinom(ei_n * ei_m, 1, 0.6), ei_n, ei_m) == 1
ei_mask <- matrix(rbinom(ei_n * ei_m, 1, 0.45), ei_n, ei_m)
ei_y <- as_incomplete(ei_theta * ei_obs)
ei_metrics <- c("rmse", "rrmse", "mae", "medae", "mape", "spearman")

# the definitions, on two complete vectors of the same length
ei_ref <- function(p, a) {
  e <- a - p
  list(
    rmse = sqrt(sum(e^2) / length(e)),
    rrmse = sqrt(sum(e^2) / sum(a^2)),
    mae = sum(abs(e)) / length(e),
    # the middle order statistic, or the mean of the two middle ones
    medae = sum(sort(abs(e))[c(floor((length(e) + 1) / 2),
                               ceiling((length(e) + 1) / 2))]) / 2,
    mape = 100 * sum(abs(e / a)) / length(e),
    # Spearman's rho is the Pearson correlation of the ranks
    spearman = cor(rank(a), rank(p))
  )
}

# Tolerance 1e-12 throughout: the reference and evaluate() use the same doubles
# and differ at most in the order of the floating-point operations (sum against
# mean, one square root against a ratio of two), a relative error of a few
# units of .Machine$double.eps on fewer than 110 terms.
ei_tol <- 1e-12

ei_expect_all <- function(predicted, actual, mask, p_ref, a_ref) {
  ref <- ei_ref(p_ref, a_ref)
  for (mt in ei_metrics) {
    expect_equal(evaluate(predicted, actual, metric = mt, mask = mask),
                 ref[[mt]], tolerance = ei_tol, info = mt)
  }
  tab <- evaluate(predicted, actual, mask = mask)
  expect_s3_class(tab, "data.frame")
  expect_identical(dim(tab), c(1L, 5L))
  # "medae" is not in the default table
  expect_identical(names(tab),
                   c("RMSE", "Rel_RMSE", "MAE", "MAPE", "Spearman_Rho"))
  expect_equal(unlist(tab, use.names = FALSE),
               unlist(ref[c("rmse", "rrmse", "mae", "mape", "spearman")],
                      use.names = FALSE),
               tolerance = ei_tol)
}

test_that("the hand-written median of the reference is the median", {
  # guards the reference itself: odd and even lengths
  # absolute errors 3, 1, 2 and 4, 1, 2, 3: small integers, exact
  expect_identical(ei_ref(c(1, 2, 3), c(4, 3, 5))$medae, 2)
  expect_identical(ei_ref(c(1, 2, 3, 4), c(5, 3, 5, 7))$medae, 2.5)
})

# Incomplete actual, no mask ----

test_that("incomplete actual, dense predicted: stored entries only", {
  expect_identical(length(ei_y@x), sum(ei_obs))
  ei_expect_all(ei_pred, ei_y, NULL, ei_pred[ei_obs], ei_theta[ei_obs])
})

test_that("incomplete actual, sparse dgCMatrix predicted: stored entries only", {
  # predictions stored at the observed entries only, and a fully stored matrix
  pred_sp <- Matrix::Matrix(ei_pred * ei_obs, sparse = TRUE)
  expect_s4_class(pred_sp, "dgCMatrix")
  ei_expect_all(pred_sp, ei_y, NULL, ei_pred[ei_obs], ei_theta[ei_obs])
  pred_full <- methods::as(ei_pred, "CsparseMatrix")
  ei_expect_all(pred_full, ei_y, NULL, ei_pred[ei_obs], ei_theta[ei_obs])
})

test_that("incomplete actual, vector predicted aligned with actual@x", {
  # column-major positions of the stored entries, built from the dense set
  pos <- which(ei_obs)
  ei_expect_all(ei_pred[pos], ei_y, NULL, ei_pred[pos], ei_theta[pos])
})

test_that("an incomplete actual is not compared on its zero-filled entries", {
  # the zero-filled RMSE over all n * m entries is a different number
  zero_filled <- sqrt(sum((ei_pred - ei_theta * ei_obs)^2) / (ei_n * ei_m))
  expect_gt(abs(evaluate(ei_pred, ei_y, "rmse") - zero_filled), 0.1)
})

# Mask, dense actual ----

test_that("dense, logical and sparse masks select the entries mask == 1", {
  sel <- ei_mask == 1
  ei_expect_all(ei_pred, ei_theta, ei_mask, ei_pred[sel], ei_theta[sel])
  ei_expect_all(ei_pred, ei_theta, sel, ei_pred[sel], ei_theta[sel])
  ei_expect_all(ei_pred, ei_theta, as_incomplete(ei_mask),
                ei_pred[sel], ei_theta[sel])
})

test_that("the three kinds of mask give the same value", {
  sel <- ei_mask == 1
  # the same entries in the same order: identical doubles, no tolerance
  dense <- evaluate(ei_pred, ei_theta, mask = ei_mask)
  expect_identical(evaluate(ei_pred, ei_theta, mask = sel), dense)
  expect_identical(evaluate(ei_pred, ei_theta, mask = as_incomplete(ei_mask)),
                   dense)
})

test_that("mask, vector predicted aligned with the masked entries", {
  pos <- which(ei_mask == 1)
  ei_expect_all(ei_pred[pos], ei_theta, ei_mask, ei_pred[pos], ei_theta[pos])
  ei_expect_all(ei_pred[pos], ei_theta, ei_mask == 1,
                ei_pred[pos], ei_theta[pos])
  ei_expect_all(ei_pred[pos], ei_theta, as_incomplete(ei_mask),
                ei_pred[pos], ei_theta[pos])
})

# Mask and incomplete actual ----

test_that("mask and incomplete actual: the intersection is evaluated", {
  sel <- ei_mask == 1 & ei_obs
  # the mask selects stored and unstored entries, so the intersection matters
  expect_true(any(ei_mask == 1 & !ei_obs))
  expect_gt(sum(sel), 10)
  ei_expect_all(ei_pred, ei_y, ei_mask, ei_pred[sel], ei_theta[sel])
  ei_expect_all(ei_pred, ei_y, ei_mask == 1, ei_pred[sel], ei_theta[sel])
  ei_expect_all(ei_pred, ei_y, as_incomplete(ei_mask),
                ei_pred[sel], ei_theta[sel])
})

test_that("mask and incomplete actual: a vector follows the masked entries", {
  # aligned with mask == 1, before the intersection with the stored entries
  pos <- which(ei_mask == 1)
  sel <- which(ei_mask == 1 & ei_obs)
  ei_expect_all(ei_pred[pos], ei_y, ei_mask, ei_pred[sel], ei_theta[sel])
  ei_expect_all(ei_pred[pos], ei_y, as_incomplete(ei_mask),
                ei_pred[sel], ei_theta[sel])
})

# No mask, actual not incomplete ----

test_that("vectors with NAs: pairwise removal, then the direct formulas", {
  set.seed(31)
  a <- rnorm(25, mean = 2)
  p <- a + rnorm(25, sd = 0.5)
  a[c(3, 17)] <- NA
  p[c(8, 17, 22)] <- NA
  ok <- !is.na(a) & !is.na(p)
  expect_identical(sum(ok), 21L)
  ei_expect_all(p, a, NULL, p[ok], a[ok])
})

test_that("dense matrices with NAs: pairwise removal over all entries", {
  a <- ei_theta
  p <- ei_pred
  a[c(2, 40, 77)] <- NA
  p[c(5, 40, 101)] <- NA
  ok <- !is.na(a) & !is.na(p)
  ei_expect_all(p, a, NULL, p[ok], a[ok])
})

test_that("complete dense matrices: all n * m entries are evaluated", {
  ei_expect_all(ei_pred, ei_theta, NULL,
                as.numeric(ei_pred), as.numeric(ei_theta))
})

# Errors ----

test_that("a mask with the wrong dimensions is an error", {
  expect_error(evaluate(ei_pred, ei_theta, mask = ei_mask[-1, ]))
  expect_error(evaluate(ei_pred, ei_theta, mask = ei_mask[, -1] == 1))
  expect_error(evaluate(ei_pred, ei_theta,
                        mask = as_incomplete(ei_mask[-1, ])))
  expect_error(evaluate(ei_pred, ei_y, mask = ei_mask[, -1]))
})

test_that("a vector predicted of the wrong length is an error", {
  # incomplete actual: one value short, and the full n * m vector (a vector is
  # never subset)
  expect_error(evaluate(ei_pred[ei_obs][-1], ei_y))
  expect_error(evaluate(as.numeric(ei_pred), ei_y))
  # mask
  pos <- which(ei_mask == 1)
  expect_error(evaluate(ei_pred[pos][-1], ei_theta, mask = ei_mask))
  expect_error(evaluate(as.numeric(ei_pred), ei_theta, mask = ei_mask))
  # mask and incomplete actual: a vector aligned with the intersection is
  # shorter than the masked entries
  sel <- which(ei_mask == 1 & ei_obs)
  expect_lt(length(sel), length(pos))
  expect_error(evaluate(ei_pred[sel], ei_y, mask = ei_mask))
})
