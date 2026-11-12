sim <- function(seed) {
  generate_simulated_data(25, 30, 2, 1, 1, sparsity = 0.5, snr = 1,
                          seed = seed)
}
lr_grid <- function(...) {
  imr_tune_grid(nuclear = c(0, NA, 4, 2), rank = c(2, 6, 2, 2), ...)
}
collect_messages <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  return(msgs)
}

test_that("imr() tuned case returns fit, data and a computed grid", {
  d <- sim(1)
  out <- suppressMessages(imr(d$Y, grid = lr_grid(), row_intercept = TRUE,
                              col_intercept = FALSE, n_cores = 1, seed = 1,
                              verbose = 0))
  expect_s3_class(out$fit, "imr_fit")
  expect_s3_class(out$data, "imr_data")
  expect_true(out$data$meta$split_data)
  expect_true(out$data$model$row_intercept)
  expect_false(out$data$model$col_intercept)
  expect_true(is.numeric(out$grid$nuclear$max))
  expect_true(all(c("fit", "params", "history", "metric", "time_secs",
                    "data", "grid") %in% names(out)))
})

test_that("imr() is reproducible with the same seed", {
  d <- sim(2)
  run <- function() suppressMessages(
    imr(d$Y, grid = lr_grid(), n_cores = 1, seed = 7, verbose = 0)
  )
  a <- run()
  b <- run()
  # same seed, same split and same deterministic solver: default tolerance
  expect_equal(a$params, b$params)
  expect_equal(a$fit$coefficients, b$fit$coefficients)
})

test_that("imr() fits once when the grid fixes every parameter", {
  d <- sim(3)
  grid <- imr_tune_grid(nuclear = 0.1, rank = 2)
  out <- suppressMessages(imr(d$Y, grid = grid, val_prop = 0, verbose = 0))
  expect_named(out, c("fit", "data", "grid"))
  expect_false(out$data$meta$split_data)
  expect_false(out$fit$meta$training)
  expect_equal(out$fit$meta$rank, 2)
  expect_true(0.1 %in% unlist(out$fit$meta$lambdas))
})

test_that("the fixed rule depends on the covariates given", {
  d <- sim(4)
  grid <- imr_tune_grid(beta = c(0, NA, 3), nuclear = 0.1, rank = 2)
  with_x <- suppressMessages(imr(d$Y, X = d$X, grid = grid, n_cores = 1,
                                 seed = 1, verbose = 0))
  expect_true("params" %in% names(with_x))
  no_x <- suppressMessages(imr(d$Y, grid = grid, n_cores = 1, seed = 1,
                               verbose = 0))
  expect_false("params" %in% names(no_x))
})

test_that("imr() metric defaults to the grid and is overridden by the argument", {
  d <- sim(5)
  grid <- lr_grid(metric = "mae")
  run <- function(...) suppressMessages(
    imr(d$Y, grid = grid, n_cores = 1, seed = 1, verbose = 0, ...)
  )
  expect_equal(run()$metric, "MAE")
  expect_equal(run(metric = "medae")$metric, "MEDAE")
  expect_error(run(metric = "not_a_metric"))
})

test_that("medae equals the median absolute error with pairwise NA removal", {
  set.seed(11)
  p <- rnorm(40)
  a <- rnorm(40)
  p[3] <- NA
  a[10] <- NA
  ref <- median(abs(a - p), na.rm = TRUE)
  expect_equal(get_metric("medae")(p, a), ref)
  expect_equal(get_metric("MedAE")(p, a), ref)
  expect_equal(unname(evaluate(p, a, metric = "medae")), ref)
  pm <- matrix(p, 8, 5)
  am <- matrix(a, 8, 5)
  expect_equal(get_metric("medae")(pm, am), ref)
  expect_equal(unname(evaluate(pm, am, metric = "MEDAE")), ref)
  expect_equal(names(evaluate(p, a)),
               c("RMSE", "Rel_RMSE", "MAE", "MAPE", "Spearman_Rho"))
})

test_that("imr() verbosity controls the messages", {
  d <- sim(6)
  tuned <- function(v) imr(d$Y, grid = lr_grid(), thresh = 1e-3, n_cores = 1,
                           seed = 1, verbose = v)
  fixed <- function(v) imr(d$Y, grid = imr_tune_grid(nuclear = 0.1, rank = 2),
                           thresh = 1e-3, verbose = v)
  expect_silent(tuned(0))
  m1 <- collect_messages(tuned(1))
  expect_length(m1, 6)
  expect_true(all(startsWith(m1, "imr: ")))
  mf <- collect_messages(fixed(1))
  expect_length(mf, 4)
  expect_true(all(startsWith(mf, "imr: ")))
  expect_gt(length(collect_messages(tuned(2))), 6)
})

test_that("imr() stops on invalid arguments in the tuned case", {
  d <- sim(7)
  run <- function(...) imr(d$Y, grid = lr_grid(), n_cores = 1, verbose = 0,
                           ...)
  expect_error(run(val_prop = 0), "val_prop")
  expect_error(run(val_prop = 1), "val_prop")
  expect_error(imr(d$Y, grid = list(nuclear = 1), verbose = 0))
  expect_error(suppressMessages(run(thresh = 0)))
})

test_that("final_fit = FALSE leaves the model fitted on the training split", {
  d <- sim(8)
  run <- function(...) suppressMessages(
    imr(d$Y, grid = lr_grid(), n_cores = 1, seed = 1, verbose = 0, ...)
  )
  expect_false(run()$fit$meta$training)
  expect_true(run(final_fit = FALSE)$fit$meta$training)
})
