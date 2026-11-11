# imr_identify() is checked against the dense map T(Theta) of the paper
# (prop:identification), computed from the full fitted matrix. The function
# itself never forms that matrix, so the two routes are independent.
#
# Tolerance 1e-8 throughout: both routes are direct linear algebra on
# well-conditioned 20 x 30 problems, so they agree to about 1e-12. 1e-8 is the
# acceptance level of the task.

identify_setup <- function(row_covariates = TRUE, col_covariates = TRUE,
                           low_rank_component = TRUE, row_intercept = TRUE,
                           col_intercept = TRUE) {
  dat <- generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, snr = 2,
                                 seed = 11)
  data <- update(
    imr_data(dat$Y, dat$X, dat$Z, val_prop = 0),
    row_covariates = row_covariates, col_covariates = col_covariates,
    low_rank_component = low_rank_component, shared_beta = FALSE,
    shared_gamma = FALSE, row_similarity = FALSE, col_similarity = FALSE,
    row_intercept = row_intercept, col_intercept = col_intercept
  )
  fit <- suppressMessages(
    imr_fit(data, 3, lambda_m = 0.5, lambda_beta = 0.1, lambda_gamma = 0.1)
  )
  list(data = data, fit = fit)
}

# dense reference: T(theta) for the augmented designs xt and zt (NULL if absent)
dense_identify <- function(theta, xt, zt) {
  qx <- if (!is.null(xt)) qr(xt)
  qz <- if (!is.null(zt)) qr(zt)
  out_x <- if (is.null(qx)) theta else qr.resid(qx, theta)
  list(
    b = if (!is.null(qx)) qr.coef(qx, theta),
    g = if (!is.null(qz)) t(qr.coef(qz, t(out_x))),
    M = if (is.null(qz)) out_x else t(qr.resid(qz, t(out_x)))
  )
}

expect_identified <- function(fit, data) {
  tol <- 1e-8
  id <- imr_identify(fit, data)
  co <- id$coefficients
  theta <- reconstruct(fit, data, trace = FALSE)$estimates

  # the fitted matrix is unchanged
  expect_equal(reconstruct(id, data, trace = FALSE)$estimates, theta,
               tolerance = tol)

  xt <- cbind(if (!is.null(co$gamma0)) rep(1, nrow(theta)),
              if (!is.null(co$beta)) data$X)
  zt <- cbind(if (!is.null(co$beta0)) rep(1, ncol(theta)),
              if (!is.null(co$gamma)) data$Z)
  ref <- dense_identify(theta, xt, zt)
  b <- rbind(co$gamma0, co$beta)
  g <- cbind(co$beta0, co$gamma)
  has_m <- length(co$d) > 0
  M <- if (has_m) co$u %*% (co$d * t(co$v))

  # equality with the dense map
  if (!is.null(xt)) expect_equal(b, ref$b, tolerance = tol, ignore_attr = TRUE)
  if (!is.null(zt)) expect_equal(g, ref$g, tolerance = tol, ignore_attr = TRUE)
  if (has_m) expect_equal(M, ref$M, tolerance = tol, ignore_attr = TRUE)

  # (C1), (C2), (C3), relative to the size of the factors
  small <- function(prod, lhs, rhs)
    expect_lt(max(abs(prod)), tol * norm(lhs, "F") * norm(rhs, "F"))
  if (has_m && !is.null(xt)) small(crossprod(xt, M), xt, M)
  if (has_m && !is.null(zt)) small(M %*% zt, M, zt)
  if (!is.null(xt) && !is.null(zt)) small(crossprod(xt, g), xt, g)
  invisible(id)
}

test_that("imr_identify() gives T(Theta) and (C1)-(C3) for the full model", {
  s <- identify_setup()
  id <- expect_identified(s$fit, s$data)
  expect_true(id$meta$identified)
  # u and v stay orthonormal and d is positive and decreasing
  r <- length(id$coefficients$d)
  expect_equal(crossprod(id$coefficients$u), diag(r), tolerance = 1e-8)
  expect_equal(crossprod(id$coefficients$v), diag(r), tolerance = 1e-8)
  expect_true(all(diff(id$coefficients$d) <= 0) && all(id$coefficients$d > 0))
})

test_that("imr_identify() is idempotent", {
  s <- identify_setup()
  id <- imr_identify(s$fit, s$data)
  id2 <- imr_identify(id, s$data)
  co <- id$coefficients
  co2 <- id2$coefficients
  expect_equal(co2[c("beta", "gamma", "beta0", "gamma0", "d")],
               co[c("beta", "gamma", "beta0", "gamma0", "d")], tolerance = 1e-8)
  # u and v are compared through M: their signs are not unique
  expect_equal(co2$u %*% (co2$d * t(co2$v)), co$u %*% (co$d * t(co$v)),
               tolerance = 1e-8)
})

test_that("identified column intercepts are the column means of Theta", {
  # holds on the scale of the stored (centered) designs, fit$coefficients
  s <- identify_setup()
  id <- imr_identify(s$fit, s$data)
  theta <- reconstruct(s$fit, s$data, trace = FALSE)$estimates
  expect_equal(id$coefficients$gamma0, colMeans(theta), tolerance = 1e-8)
})

test_that("imr_identify() handles absent components", {
  # row side only
  s <- identify_setup(col_covariates = FALSE, row_intercept = FALSE)
  expect_identified(s$fit, s$data)
  # column side only
  s <- identify_setup(row_covariates = FALSE, col_intercept = FALSE)
  expect_identified(s$fit, s$data)
  # no low-rank component
  s <- identify_setup(low_rank_component = FALSE)
  expect_identified(s$fit, s$data)
  # intercepts and the low-rank component only
  s <- identify_setup(row_covariates = FALSE, col_covariates = FALSE)
  expect_identified(s$fit, s$data)
  # covariates without intercepts
  s <- identify_setup(row_intercept = FALSE, col_intercept = FALSE)
  expect_identified(s$fit, s$data)
})

test_that("the identified fit works with summary() and reconstruct_partial()", {
  s <- identify_setup()
  id <- imr_identify(s$fit, s$data)
  expect_output(print(id))
  expect_output(summary(id))
  y <- s$data$Y
  expect_equal(reconstruct_partial(id, s$data, y@i, y@p),
               reconstruct_partial(s$fit, s$data, y@i, y@p), tolerance = 1e-8)
  # the residuals are untouched, so the fit statistics are the same
  expect_identical(id$residuals, s$fit$residuals)
  expect_identical(id$meta$sum_squares$rss, s$fit$meta$sum_squares$rss)
})

test_that("imr_identify() stops on shared effects", {
  s <- identify_setup()
  # the flag alone triggers the error: no shared fit is needed
  s$fit$model$shared_beta <- TRUE
  expect_error(imr_identify(s$fit, s$data), "shared effects")
})
