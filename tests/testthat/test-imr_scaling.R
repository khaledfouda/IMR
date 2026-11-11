sim <- generate_simulated_data(30, 40, 3, 2, 2, sparsity = 0.5, snr = 1, seed = 31)
sd_x <- apply(sim$X, 2, sd)
sd_z <- apply(sim$Z, 2, sd)

test_that("imr_data centers and scales one stored copy of X and Z", {
  data <- imr_data(sim$Y, sim$X, sim$Z)
  expect_equal(data$X, scale(sim$X), ignore_attr = TRUE)
  expect_equal(data$Z, scale(sim$Z), ignore_attr = TRUE)
  expect_equal(data$x_center, colMeans(sim$X))
  expect_equal(data$x_scale, sd_x)
  expect_equal(data$z_center, colMeans(sim$Z))
  expect_equal(data$z_scale, sd_z)
  expect_equal(data$LX, 1 / svd(data$X)$d[1]^2)
  expect_equal(data$LZ, 1 / svd(data$Z)$d[1]^2)
  # the original matrix is the stored copy times the scales plus the centers
  expect_equal(sweep(sweep(data$X, 2, data$x_scale, "*"), 2, data$x_center, "+"), sim$X)
  expect_false("X_raw" %in% names(data))

  plain <- imr_data(sim$Y, sim$X, sim$Z, center = FALSE, scale = FALSE)
  expect_equal(plain$X, sim$X)
  expect_null(plain$x_center)
  expect_null(plain$x_scale)
  expect_equal(plain$LX, 1 / svd(sim$X)$d[1]^2)

  # update() switches centering and scaling on the stored copy, from any state to any other
  for (to_center in c(TRUE, FALSE)) for (to_scale in c(TRUE, FALSE)) {
    target <- imr_data(sim$Y, sim$X, sim$Z, center = to_center, scale = to_scale)
    for (start in list(data, plain)) {
      switched <- update(start, center = to_center, scale = to_scale)
      for (field in c("X", "x_center", "x_scale", "LX", "Z", "z_center", "z_scale", "LZ"))
        expect_equal(switched[[field]], target[[field]], info = field)
    }
  }
  # intercepts and other flags leave the stored copies alone
  expect_identical(update(data, row_intercept = TRUE, col_intercept = TRUE)$X, data$X)
})

test_that("imr_data stops on constant covariate columns when scaling", {
  expect_error(imr_data(sim$Y, cbind(sim$X, 0.1)), "Constant columns of X cannot be scaled: 3")
  expect_error(imr_data(sim$Y, Z = cbind(1, sim$Z)), "Constant columns of Z cannot be scaled: 1")
  unscaled <- imr_data(sim$Y, cbind(sim$X, 0.1), scale = FALSE)
  expect_error(update(unscaled, scale = TRUE), "Constant columns of X cannot be scaled: 3")
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

test_that("without intercepts, centering is part of the fitted covariate effects", {
  data <- imr_data(sim$Y, sim$X, sim$Z)
  set.seed(1)
  fit <- suppressMessages(imr_fit(data, 3, lambda_m = 1, lambda_beta = 0.5, lambda_gamma = 0.5))
  cf <- coef(fit)
  expect_null(cf$beta0)
  expect_null(cf$gamma0)
  expect_equal(cf$beta, fit$coefficients$beta / sd_x)
  rec <- reconstruct(fit, data, trace = FALSE)
  # (X - 1 xbar') beta and Gamma (Z - 1 zbar')' on the original scale
  expect_equal(rec$xbeta, sweep(sim$X, 2, colMeans(sim$X)) %*% cf$beta, tolerance = 1e-10)
  expect_equal(rec$gammaz, cf$gamma %*% t(sweep(sim$Z, 2, colMeans(sim$Z))), tolerance = 1e-10)
  expect_equal(rec$estimates, rec$xbeta + rec$gammaz + rec$M, tolerance = 1e-10)
})

test_that("imr_fit returns the objective at every iteration, and it never increases", {
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
  fit_at <- function(k, ...) {
    set.seed(1)
    suppressMessages(imr_fit(
      data, rank = 4, lambda_m = lam[["m"]], lambda_beta = lam[["beta"]],
      lambda_gamma = lam[["gamma"]],
      convergence = imr_convergence(maxit = k, thresh = 1e-300, ls_initial = FALSE), ...
    ))
  }

  expect_null(fit_at(3)$objective)
  obj <- fit_at(20, return_objective = TRUE)$objective
  expect_length(obj, 20)
  # with the seed fixed, the k-th iterate is the fit stopped at maxit = k
  for (k in c(1, 7, 20)) expect_equal(obj[k], objective(fit_at(k)), tolerance = 1e-10)
  expect_true(all(diff(obj) <= 1e-10 * obj[1]))
})

test_that("the returned objective includes the similarity penalties", {
  sr <- imr_similarity("rbf", d = as.matrix(dist(seq_len(nrow(sim$Y)))) / 10)
  sc <- imr_similarity("rbf", d = as.matrix(dist(seq_len(ncol(sim$Y)))) / 10)
  data <- imr_data(sim$Y, sim$X, sim$Z, similarity_rows = sr, similarity_cols = sc)
  set.seed(1)
  fit <- suppressMessages(imr_fit(data, 3, lambda_m = 2, lambda_beta = 0.5, lambda_gamma = 0.5,
                                  convergence = imr_convergence(maxit = 10),
                                  return_objective = TRUE))
  s <- fit$coefficients
  A <- s$u %*% diag(sqrt(s$d), length(s$d))
  B <- s$v %*% diag(sqrt(s$d), length(s$d))
  S_r <- sr$U %*% (sr$d * t(sr$U))
  S_c <- sc$U %*% (sc$d * t(sc$U))
  theta <- data$X %*% s$beta + s$gamma %*% t(data$Z) + A %*% t(B)
  expected <- 0.5 * sum(((sim$Y != 0) * (sim$Y - theta))^2) +
    2 / 2 * (sum(diag(t(A) %*% S_r %*% A)) + sum(diag(t(B) %*% S_c %*% B))) +
    0.5 * sum(abs(s$beta)) + 0.5 * sum(abs(s$gamma))
  expect_equal(fit$objective[fit$meta$n_iter], expected, tolerance = 1e-10)
})
