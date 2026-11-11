sim <- generate_simulated_data(30, 40, 3, 2, 2, sparsity = 0.5, snr = 1, seed = 31)
sd_x <- apply(sim$X, 2, sd)
sd_z <- apply(sim$Z, 2, sd)

test_that("imr_data scales the covariates and centers them only while the matching intercept is in the model", {
  data <- imr_data(sim$Y, sim$X, sim$Z)
  expect_equal(data$X_raw, sim$X)
  expect_equal(data$Z_raw, sim$Z)
  # no intercepts: scaled, not centered
  expect_false(data$meta$x_centered)
  expect_false(data$meta$z_centered)
  expect_equal(data$X, sweep(sim$X, 2, sd_x, "/"))
  expect_equal(data$Z, sweep(sim$Z, 2, sd_z, "/"))
  expect_equal(data$LX, 1 / svd(data$X)$d[1]^2)
  expect_equal(data$LZ, 1 / svd(data$Z)$d[1]^2)

  # gamma0 (column intercepts) absorbs the means of X, beta0 (row intercepts) those of Z
  with_col <- update(data, col_intercept = TRUE)
  expect_true(with_col$meta$x_centered)
  expect_false(with_col$meta$z_centered)
  expect_equal(with_col$X, scale(sim$X), ignore_attr = TRUE)
  expect_equal(with_col$LX, 1 / svd(with_col$X)$d[1]^2)
  expect_false(isTRUE(all.equal(with_col$LX, data$LX)))
  expect_equal(with_col$Z, data$Z)

  with_row <- update(data, row_intercept = TRUE)
  expect_true(with_row$meta$z_centered)
  expect_false(with_row$meta$x_centered)
  expect_equal(with_row$Z, scale(sim$Z), ignore_attr = TRUE)
  expect_equal(with_row$LZ, 1 / svd(with_row$Z)$d[1]^2)

  # switching the intercept off again restores the uncentered design
  back <- update(with_col, col_intercept = FALSE)
  expect_equal(back$X, data$X)
  expect_equal(back$LX, data$LX)

  # center = FALSE never centers; scale = FALSE keeps the original units
  no_center <- update(imr_data(sim$Y, sim$X, sim$Z, center = FALSE),
                      row_intercept = TRUE, col_intercept = TRUE)
  expect_false(no_center$meta$x_centered)
  expect_false(no_center$meta$z_centered)
  expect_equal(no_center$X, data$X)
  no_scale <- update(imr_data(sim$Y, sim$X, sim$Z, scale = FALSE), col_intercept = TRUE)
  expect_equal(no_scale$x_scale, c(1, 1))
  expect_equal(no_scale$X, sweep(sim$X, 2, colMeans(sim$X)))
})

test_that("imr_data stops on constant covariate columns when scaling", {
  expect_error(imr_data(sim$Y, cbind(sim$X, 0.1)), "Constant columns of X cannot be scaled: 3")
  expect_error(imr_data(sim$Y, Z = cbind(1, sim$Z)), "Constant columns of Z cannot be scaled: 1")
  expect_no_error(imr_data(sim$Y, cbind(sim$X, 0.1), scale = FALSE))
})

test_that("coef() and reconstruct() report the coefficients on the original scale", {
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  set.seed(1)
  fit <- suppressMessages(imr_fit(data, 3, lambda_m = 1, lambda_beta = 0.5, lambda_gamma = 0.5))
  s <- fit$coefficients
  n <- nrow(sim$Y)
  m <- ncol(sim$Y)

  beta <- s$beta / sd_x
  gamma <- sweep(s$gamma, 2, sd_z, "/")
  cf <- coef(fit)
  expect_equal(cf$beta, beta)
  expect_equal(cf$gamma, gamma)
  expect_equal(cf$gamma0, s$gamma0 - as.vector(t(beta) %*% colMeans(sim$X)))
  expect_equal(cf$beta0, s$beta0 - as.vector(gamma %*% colMeans(sim$Z)))
  expect_equal(cf[c("u", "d", "v")], s[c("u", "d", "v")])

  # the fitted matrix, built on the scale of the designs
  theta_s <- data$X %*% s$beta + s$gamma %*% t(data$Z) + s$u %*% (s$d * t(s$v)) +
    outer(s$beta0, rep(1, m)) + outer(rep(1, n), s$gamma0)
  rec <- reconstruct(fit, data, trace = FALSE)
  expect_equal(rec$estimates, theta_s, tolerance = 1e-10)
  expect_equal(rec$xbeta + rec$gammaz + rec$M + outer(rec$beta0, rep(1, m)) +
                 outer(rep(1, n), rec$gamma0), rec$estimates, tolerance = 1e-10)
  expect_equal(rec$beta, cf$beta)
  expect_equal(rec$gamma0, cf$gamma0)
  # reconstruct_partial() works on the designs and predicts the same entries
  expect_equal(reconstruct_partial(fit, data, data$Y@i, data$Y@p),
               rec$estimates[as.matrix(data$Y) != 0], tolerance = 1e-10)
  expect_output(summary(fit))
})

test_that("fits do not depend on the location and scale of the covariates", {
  a <- c(10, 0.5)
  b <- c(-3, 7)
  X2 <- sweep(sweep(sim$X, 2, a, "*"), 2, b, "+")
  data1 <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  data2 <- update(imr_data(sim$Y, X2, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  expect_equal(data2$X, data1$X)

  set.seed(1)
  fit1 <- suppressMessages(imr_fit(data1, 3, lambda_m = 1, lambda_beta = 0.5, lambda_gamma = 0.5))
  set.seed(1)
  fit2 <- suppressMessages(imr_fit(data2, 3, lambda_m = 1, lambda_beta = 0.5, lambda_gamma = 0.5))
  # tolerance 1e-8: the two designs differ by rounding only
  expect_equal(reconstruct(fit2, data2, trace = FALSE)$estimates,
               reconstruct(fit1, data1, trace = FALSE)$estimates, tolerance = 1e-8)
  expect_equal(coef(fit2)$beta, coef(fit1)$beta / a, tolerance = 1e-8)
})

test_that("imr_fit stops when the intercepts are changed without update()", {
  data <- imr_data(sim$Y, sim$X, sim$Z)
  data$model$col_intercept <- TRUE
  expect_error(imr_fit(data, 2), "does not match the intercepts")
})

test_that("the objective never increases on standardized designs with intercepts", {
  data <- update(imr_data(sim$Y, sim$X, sim$Z), row_intercept = TRUE, col_intercept = TRUE)
  Y <- as.matrix(data$Y)
  obs <- Y != 0
  lam <- c(m = 1, beta = 0.5, gamma = 0.5)
  # dense brute-force objective on the scale of the designs, where the penalties apply
  objective <- function(fit) {
    s <- fit$coefficients
    theta <- data$X %*% s$beta + s$gamma %*% t(data$Z) + s$u %*% (s$d * t(s$v)) +
      outer(s$beta0, rep(1, ncol(Y))) + outer(rep(1, nrow(Y)), s$gamma0)
    0.5 * sum((obs * (Y - theta))^2) + lam[["m"]] * sum(s$d) +
      lam[["beta"]] * sum(abs(s$beta)) + lam[["gamma"]] * sum(abs(s$gamma))
  }
  # with the seed fixed, the k-th iterate is the fit stopped at maxit = k
  obj <- vapply(1:20, function(k) {
    set.seed(1)
    objective(suppressMessages(imr_fit(
      data, rank = 4, lambda_m = lam[["m"]], lambda_beta = lam[["beta"]],
      lambda_gamma = lam[["gamma"]],
      convergence = imr_convergence(maxit = k, thresh = 1e-300, ls_initial = FALSE)
    )))
  }, numeric(1))
  expect_true(all(diff(obj) <= 1e-10 * obj[1]))
})
