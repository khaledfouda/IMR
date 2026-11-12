# Independent tests of the KKT screening at the cold start of imr_solver().
# With lambda_m at or above sigma_1 of the residual of the model without M,
# M = 0 solves the problem and the fit is returned without iterating.
#
# Problem used throughout, unless a test says otherwise: both intercepts, two
# row covariates and one column covariate, about 60% of a 14 x 18 matrix observed.

screening_sim <- function() {
  generate_simulated_data(14, 18, 2, 2, 1, sparsity = 0.4, snr = 1, seed = 311)
}

screening_data <- function(sim, low_rank = TRUE, scale_y = 1) {
  update(imr_data(scale_y * sim$Y, sim$X, sim$Z), low_rank_component = low_rank,
         row_intercept = TRUE, col_intercept = TRUE)
}

# sigma_1 of the dense residual (zero at the unobserved entries) of a fit of the
# model without M: the smallest lambda_m at which M = 0 is a solution
sigma_1 <- function(fit) svd(as.matrix(fit$residuals))$d[1]

test_that("a screened fit with Lasso penalties equals the penalized model without M", {
  sim <- screening_sim()
  conv <- imr_convergence(maxit = 1e5, thresh = 1e-12)
  null_fit <- imr_fit(screening_data(sim, low_rank = FALSE), lambda_beta = 1,
                      lambda_gamma = 0.5, convergence = conv)
  # the penalties are active without removing every coefficient
  expect_true(any(null_fit$coefficients$beta == 0) && any(null_fit$coefficients$beta != 0))

  expect_no_message(
    fit <- imr_fit(screening_data(sim), rank = 3, lambda_m = 1.05 * sigma_1(null_fit),
                   lambda_beta = 1, lambda_gamma = 0.5, convergence = conv)
  )
  expect_equal(fit$meta$n_iter, 0)
  expect_true(fit$meta$converged)
  expect_true(all(fit$coefficients$d == 0))
  expect_equal(ncol(fit$coefficients$u), 1)
  expect_equal(ncol(fit$coefficients$v), 1)
  # the same deterministic fit of the model without M on both sides: default tolerance
  expect_equal(fit$residuals, null_fit$residuals)
  covariate_terms <- c("beta", "gamma", "beta0", "gamma0")
  expect_equal(fit$coefficients[covariate_terms], null_fit$coefficients[covariate_terms])
})

test_that("a screened fit with Lasso penalties satisfies the KKT conditions of the full problem", {
  sim <- screening_sim()
  lambda_m <- 6 # above sigma_1 of the residual of the penalized model without M (about 5.4)
  lambda_beta <- 1
  lambda_gamma <- 0.5
  fit <- imr_fit(screening_data(sim), rank = 3, lambda_m = lambda_m, lambda_beta = lambda_beta,
                 lambda_gamma = lambda_gamma,
                 convergence = imr_convergence(maxit = 1e5, thresh = 1e-12))
  expect_equal(fit$meta$n_iter, 0)

  # dense residual, and the designs on the scale where the penalties apply
  resid <- as.matrix(fit$residuals)
  x_s <- sweep(sweep(sim$X, 2, fit$scaling$x_center), 2, fit$scaling$x_scale, "/")
  z_s <- sweep(sweep(sim$Z, 2, fit$scaling$z_center), 2, fit$scaling$z_scale, "/")
  beta <- fit$coefficients$beta
  gamma <- fit$coefficients$gamma
  grad_beta <- crossprod(x_s, resid)
  grad_gamma <- resid %*% z_s

  # tolerance 1e-3 (absolute, against gradients of order 1): the sub-fit stops at
  # a relative change of 1e-12, which leaves the stationarity conditions
  # satisfied to about 1e-5 on this problem
  tol <- 1e-3
  # M = 0: the subgradient condition of the nuclear norm
  expect_lte(svd(resid)$d[1], lambda_m)
  # unpenalized intercepts
  expect_lt(max(abs(rowSums(resid))), tol)
  expect_lt(max(abs(colSums(resid))), tol)
  # Lasso: X_s' E = lambda_beta * sign(beta) on the support, |.| <= lambda_beta off it
  expect_lt(max(abs(grad_beta[beta != 0] - lambda_beta * sign(beta[beta != 0]))), tol)
  expect_lt(max(abs(grad_beta[beta == 0])), lambda_beta + tol)
  expect_lt(max(abs(grad_gamma[gamma != 0] - lambda_gamma * sign(gamma[gamma != 0]))), tol)
  expect_lt(max(abs(grad_gamma[gamma == 0])), lambda_gamma + tol)
})

test_that("just below the critical lambda_m the penalized fit iterates to a nonzero M", {
  sim <- screening_sim()
  null_fit <- imr_fit(screening_data(sim, low_rank = FALSE), lambda_beta = 1, lambda_gamma = 0.5,
                      convergence = imr_convergence(maxit = 1e5, thresh = 1e-12))
  fit <- suppressMessages(
    imr_fit(screening_data(sim), rank = 3, lambda_m = 0.95 * sigma_1(null_fit),
            lambda_beta = 1, lambda_gamma = 0.5,
            convergence = imr_convergence(maxit = 1e4, thresh = 1e-8))
  )
  expect_gt(fit$meta$n_iter, 0)
  expect_gt(fit$coefficients$d[1], 0)
})

# In the four tests below, lambda_m = 10 * ||Y||_F. The residual of any model
# without M that does not increase the least-squares loss from the zero start
# has sigma_1 <= ||Y||_F, so this value is above the critical one in every case.

test_that("a fit with ls_initial = FALSE is not screened above the critical lambda_m", {
  sim <- screening_sim()
  data <- screening_data(sim)
  lambda_m <- 10 * sqrt(sum(sim$Y^2))
  # maxit = 1e4 lets the sub-fit without M stop by the rule, which screening requires
  screened <- imr_fit(data, rank = 3, lambda_m = lambda_m,
                      convergence = imr_convergence(maxit = 1e4, thresh = 1e-8))
  expect_equal(screened$meta$n_iter, 0)
  fit <- suppressMessages(
    imr_fit(data, rank = 3, lambda_m = lambda_m,
            convergence = imr_convergence(maxit = 1e4, thresh = 1e-8, ls_initial = FALSE))
  )
  expect_gt(fit$meta$n_iter, 0)
})

test_that("a warm start with nonzero d is not screened above the critical lambda_m", {
  sim <- screening_sim()
  data <- screening_data(sim)
  lambda_m <- 10 * sqrt(sum(sim$Y^2))
  start <- suppressMessages(
    imr_fit(data, rank = 3, lambda_m = lambda_m / 1e3,
            convergence = imr_convergence(maxit = 1e4, thresh = 1e-8))
  )
  expect_gt(start$coefficients$d[1], 0)
  fit <- suppressMessages(
    imr_fit(data, rank = 3, lambda_m = lambda_m, warm_start = start,
            convergence = imr_convergence(maxit = 5, thresh = 1e-8))
  )
  expect_gt(fit$meta$n_iter, 0)
})

test_that("a model with similarity matrices is not screened, and is once they are switched off", {
  sim <- screening_sim()
  sr <- imr_similarity("rbf", d = as.matrix(dist(seq_len(nrow(sim$Y)))) / 5)
  sc <- imr_similarity("rbf", d = as.matrix(dist(seq_len(ncol(sim$Y)))) / 5)
  data <- imr_data(sim$Y, similarity_rows = sr, similarity_cols = sc)
  lambda_m <- 10 * sqrt(sum(sim$Y^2))
  conv <- imr_convergence(maxit = 5, thresh = 1e-8)

  both <- suppressMessages(imr_fit(data, rank = 3, lambda_m = lambda_m, convergence = conv))
  expect_gt(both$meta$n_iter, 0)
  # one similarity matrix is enough to leave the nuclear-norm case
  rows_only <- suppressMessages(
    imr_fit(update(data, col_similarity = FALSE), rank = 3, lambda_m = lambda_m, convergence = conv)
  )
  expect_gt(rows_only$meta$n_iter, 0)
  none <- imr_fit(update(data, row_similarity = FALSE, col_similarity = FALSE), rank = 3,
                  lambda_m = lambda_m, convergence = conv)
  expect_equal(none$meta$n_iter, 0)
})

test_that("a fit whose sub-fit without M reaches maxit is not screened", {
  sim <- screening_sim()
  data <- screening_data(sim)
  lambda_m <- 10 * sqrt(sum(sim$Y^2))
  # one iteration cannot satisfy the stopping rule of the model with covariates
  fit <- suppressMessages(
    imr_fit(data, rank = 3, lambda_m = lambda_m,
            convergence = imr_convergence(maxit = 1, thresh = 1e-8))
  )
  expect_gt(fit$meta$n_iter, 0)
  # the same call with enough iterations for the sub-fit is screened
  screened <- imr_fit(data, rank = 3, lambda_m = lambda_m,
                      convergence = imr_convergence(maxit = 1e4, thresh = 1e-8))
  expect_equal(screened$meta$n_iter, 0)
})

test_that("a warm start with all d equal to zero gives the cold fit, without a warning", {
  sim <- screening_sim()
  data <- screening_data(sim)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-8)
  screened <- imr_fit(data, rank = 3, lambda_m = 10 * sqrt(sum(sim$Y^2)), lambda_beta = 1,
                      lambda_gamma = 0.5, convergence = conv)
  expect_true(all(screened$coefficients$d == 0))

  # lambda_m = 2 is well below the critical value (about 5.4): both fits iterate
  expect_no_warning(
    warm <- suppressMessages(imr_fit(data, rank = 3, lambda_m = 2, lambda_beta = 1,
                                     lambda_gamma = 0.5, convergence = conv,
                                     warm_start = screened))
  )
  cold <- suppressMessages(imr_fit(data, rank = 3, lambda_m = 2, lambda_beta = 1,
                                   lambda_gamma = 0.5, convergence = conv))
  expect_gt(cold$meta$n_iter, 0)
  # the warm start is discarded, so both fits run the same computation: default tolerance
  expect_equal(warm$coefficients, cold$coefficients)
  expect_equal(warm$meta$n_iter, cold$meta$n_iter)
})

test_that("a warm start without d still raises a warning", {
  sim <- screening_sim()
  data <- screening_data(sim)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-8)
  start <- imr_fit(data, rank = 3, lambda_m = 10 * sqrt(sum(sim$Y^2)), convergence = conv)$coefficients
  start$d <- NULL
  expect_warning(
    suppressMessages(imr_fit(data, rank = 3, lambda_m = 2, convergence = conv, warm_start = start))
  )
})

test_that("imr_tune() starts the lambda_m path at a fit with M = 0, without a warning", {
  # row covariates and column intercepts only, on a 20 x 24 matrix
  sim <- generate_simulated_data(20, 24, 2, 2, 0, sparsity = 0.4, snr = 1, seed = 317)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-7)
  data <- update(imr_data(sim$Y, sim$X, val_prop = 0.2, seed = 5), col_intercept = TRUE)
  grid <- imr_tune_grid(beta = 0, gamma = 0, nuclear = c(0, NA, 5, 2), rank = c(2, 4, 2, 2))
  grid <- imr_set_grid_limits(data, grid, convergence = conv, verify_iter = 0)

  for (fast in c(TRUE, FALSE)) {
    expect_no_warning(
      out <- suppressMessages(imr_tune(data, grid, final_fit = FALSE, fast_nuclear = fast,
                                       convergence = conv, verbose = 0, n_cores = 1, seed = 1))
    )
    expect_equal(out$history$rank_out[1], 0)
    # the first point of the path is the grid maximum: the same number on both sides
    expect_equal(out$history$lambda_m[1], grid$nuclear$max)
    # the path leaves M = 0 afterwards: the cold restart did not stall the tuning
    expect_true(any(out$history$rank_out > 0))
  }
})

test_that("return_objective = TRUE on a screened fit returns an empty objective", {
  sim <- screening_sim()
  fit <- imr_fit(screening_data(sim), rank = 3, lambda_m = 10 * sqrt(sum(sim$Y^2)),
                 convergence = imr_convergence(maxit = 1e4, thresh = 1e-8),
                 return_objective = TRUE)
  expect_equal(fit$meta$n_iter, 0)
  expect_length(fit$objective, 0)
})

test_that("screening is equivariant to the scale of Y and lambda_m", {
  sim <- screening_sim()
  null_fit <- imr_fit(screening_data(sim, low_rank = FALSE),
                      convergence = imr_convergence(maxit = 1e5, thresh = 1e-12))
  critical <- sigma_1(null_fit)
  conv <- imr_convergence(maxit = 1e4, thresh = 1e-10)
  n_iter <- function(scale_y, factor) {
    suppressMessages(
      imr_fit(screening_data(sim, scale_y = scale_y), rank = 3,
              lambda_m = scale_y * factor * critical, convergence = conv)
    )$meta$n_iter
  }
  # 3% on either side of the critical value, far more than the error of the
  # sub-fit, so the comparison does not depend on rounding
  for (scale_y in c(1, 1e3, 1e-3)) {
    expect_equal(n_iter(scale_y, 1.03), 0)
    expect_gt(n_iter(scale_y, 0.97), 0)
  }
})
