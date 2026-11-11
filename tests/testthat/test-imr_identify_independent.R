# Independent checks of imr_identify() (prop:identification, main.tex). The
# references are the dense fitted matrix from reconstruct() and identities that
# the identified parts must satisfy; nothing here reuses the implementation.
#
# Tolerance 1e-8 throughout: every comparison is direct dense linear algebra on
# 20 x 30 problems with rank-3 factors, so routes agree to about 1e-12 and
# 1e-8 leaves headroom for conditioning without hiding a real error.

ident_sim <- function(seed) {
  generate_simulated_data(20, 30, 3, 3, 2, sparsity = 0.5, snr = 2,
                          seed = seed)
}

ident_fit <- function(dat, center = TRUE, scale = TRUE, val_prop = 0,
                      training = FALSE, X = dat$X, Z = dat$Z,
                      row_intercept = TRUE, col_intercept = TRUE) {
  data <- imr_data(dat$Y, X, Z, val_prop = val_prop, center = center,
                   scale = scale)
  data <- update(
    data, row_covariates = TRUE, col_covariates = TRUE,
    low_rank_component = TRUE, row_intercept = row_intercept,
    col_intercept = col_intercept, shared_beta = FALSE, shared_gamma = FALSE,
    row_similarity = FALSE, col_similarity = FALSE
  )
  fit <- suppressMessages(
    imr_fit(data, 3, lambda_m = 0.5, lambda_beta = 0.1, lambda_gamma = 0.1,
            training = training)
  )
  list(data = data, fit = fit)
}

test_that("the fitted matrix is unchanged without centering or scaling", {
  dat <- ident_sim(21)
  for (cs in list(c(FALSE, FALSE), c(TRUE, FALSE))) {
    s <- ident_fit(dat, center = cs[1], scale = cs[2])
    id <- imr_identify(s$fit, s$data)
    expect_equal(reconstruct(id, s$data, trace = FALSE)$estimates,
                 reconstruct(s$fit, s$data, trace = FALSE)$estimates,
                 tolerance = 1e-8)
  }
})

test_that("coef() of the identified fit reproduces the fit on the original X and Z", {
  dat <- ident_sim(22)
  # keep copies: the stored designs may be transformed in place
  X0 <- dat$X * 1
  Z0 <- dat$Z * 1
  s <- ident_fit(dat)
  id <- imr_identify(s$fit, s$data)
  co <- coef(id)
  M <- id$coefficients$u %*% (id$coefficients$d * t(id$coefficients$v))
  n <- nrow(X0)
  m <- nrow(Z0)
  theta <- outer(co$beta0, rep(1, m)) + outer(rep(1, n), co$gamma0) +
    X0 %*% co$beta + tcrossprod(co$gamma, Z0) + M
  expect_equal(theta, reconstruct(s$fit, s$data, trace = FALSE)$estimates,
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("the identified parts satisfy the Pythagoras identity", {
  s <- ident_fit(ident_sim(23))
  id <- imr_identify(s$fit, s$data)
  co <- id$coefficients
  theta <- reconstruct(s$fit, s$data, trace = FALSE)$estimates
  xt <- cbind(1, s$data$X)
  zt <- cbind(1, s$data$Z)
  b <- rbind(co$gamma0, co$beta)
  g <- cbind(co$beta0, co$gamma)
  M <- co$u %*% (co$d * t(co$v))
  xb <- xt %*% b
  gz <- tcrossprod(g, zt)
  # the three parts add up to Theta ...
  expect_equal(xb + gz + M, theta, tolerance = 1e-8, ignore_attr = TRUE)
  # ... and are mutually orthogonal in the Frobenius inner product
  expect_equal(norm(theta, "F")^2,
               norm(xb, "F")^2 + norm(gz, "F")^2 + norm(M, "F")^2,
               tolerance = 1e-8)
})

test_that("a duplicated column of X is rejected", {
  dat <- ident_sim(24)
  X <- cbind(dat$X, dat$X[, 1])
  s <- ident_fit(dat, center = FALSE, scale = FALSE, X = X)
  expect_error(imr_identify(s$fit, s$data), "rank deficient")
})

test_that("a constant column of X with a column intercept is rejected", {
  dat <- ident_sim(25)
  X <- cbind(dat$X, 1)
  s <- ident_fit(dat, center = FALSE, scale = FALSE, X = X,
                 col_intercept = TRUE)
  expect_error(imr_identify(s$fit, s$data), "rank deficient")
})

test_that("a constant column of Z with a row intercept is rejected", {
  dat <- ident_sim(26)
  Z <- cbind(dat$Z, 1)
  s <- ident_fit(dat, center = FALSE, scale = FALSE, Z = Z,
                 row_intercept = TRUE)
  expect_error(imr_identify(s$fit, s$data), "rank deficient")
})

test_that("a fit with training = TRUE is accepted and has no sum_squares", {
  dat <- ident_sim(27)
  s <- ident_fit(dat, val_prop = 0.3, training = TRUE)
  expect_null(s$fit$meta$sum_squares)
  id <- imr_identify(s$fit, s$data)
  expect_true(id$meta$identified)
  expect_null(id$meta$sum_squares)
  expect_equal(reconstruct(id, s$data, trace = FALSE)$estimates,
               reconstruct(s$fit, s$data, trace = FALSE)$estimates,
               tolerance = 1e-8)
})
