#-----------------------------------------------------------
#' @export
nlrr.cv <- function(
  inp.dat,
  row_intercept = FALSE,
  col_intercept = FALSE,
  lambda_beta = NULL,
  lambda_gamma = NULL,
  hpar = get_imr_default_hparams(),
  error_function = error_metric$rmse,
  thresh = 1e-6,
  maxit = 300,
  verbose = 0,
  seed = NULL,
  ls_initial = FALSE
) {
  #-------------------
  stopifnot(is.Incomplete(inp.dat$Y))
  stopifnot(is.Incomplete(inp.dat$y_train))
  stopifnot(is.Incomplete(inp.dat$y_valid))
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  #-------------------------------
  # set flags
  beta_flag <- !(is.null(inp.dat$Xq))
  gamma_flag <- !(is.null(inp.dat$Zq))
  # if neither beta or gamma are provided then send to cv_M
  if (!(beta_flag | gamma_flag)) {
    return(IMR::imr.cv_M(
      y_train = inp.dat$y_train,
      y_valid = inp.dat$y_valid,
      Y_full = inp.dat$Y,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      hpar = hpar,
      error_function = error_function,
      thresh = thresh,
      trace = verbose > 0,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed
    ))
  }

  # obtain upperbounds to the lambda hyperparameters
  if (beta_flag & is.null(hpar$beta$lambda_max) & is.null(lambda_beta)) {
    hpar$beta$lambda_max <- get_lambda_lasso_max(
      y_train = inp.dat$y_train,
      X = inp.dat$Xq,
      y_valid = inp.dat$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 100,
      verbose = verbose
    )
  }
  if (gamma_flag & is.null(hpar$gamma$lambda_max) & is.null(lambda_gamma)) {
    hpar$gamma$lambda_max <- get_lambda_lasso_max(
      y_train = inp.dat$y_train,
      Z = inp.dat$Zq,
      y_valid = inp.dat$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 100,
      verbose = verbose
    )
  }

  if (beta_flag & is.null(lambda_beta)) {
    lambda_beta_grid <- seq(
      from = hpar$beta$lambda_max,
      to = 0,
      length.out = hpar$beta$n.lambda
    )
  } else {
    lambda_beta_grid <- c(if (is.null(lambda_beta)) 0 else lambda_beta)
  }

  if (gamma_flag & is.null(lambda_gamma)) {
    lambda_gamma_grid <- seq(
      from = hpar$gamma$lambda_max,
      to = 0,
      length.out = hpar$gamma$n.lambda
    )
  } else {
    lambda_gamma_grid <- c(if (is.null(lambda_gamma)) 0 else lambda_gamma)
  }

  #---------------------------
  # parallel setup
  grid <- list(
    lambda_beta = lambda_beta_grid,
    lambda_gamma = lambda_gamma_grid
  )

  nllr.fit.step <- function(lambda_beta, lambda_gamma, ...) {
    if (!is.null(seed)) set.seed(seed)
    fit <- IMR::imr.fit_no_low_rank(
      Y = inp.dat$y_train,
      X = inp.dat$Xq,
      Z = inp.dat$Zq,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      thresh = thresh,
      maxit = maxit,
      trace = FALSE
    )
    estim <- 0
    if (!is.null(fit$beta)) {
      estim <- partial_crossprod(inp.dat$Xq, fit$beta, inp.dat$y_valid@i, inp.dat$y_valid@p)
    }
    if (!is.null(fit$gamma)) {
      estim <- estim + partial_crossprod(fit$gamma, inp.dat$Zq, inp.dat$y_valid@i, inp.dat$y_valid@p, TRUE)
    }
    if (!is.null(fit$beta0)) {
      estim <- estim + partial_crossprod(
        matrix(fit$beta0, ncol = 1),
        matrix(1, 1, ncol(inp.dat$y_train)),
        inp.dat$y_valid@i, inp.dat$y_valid@p
      )
    }
    # add_to_rows_inplace_cpp(y_valid@x, y_valid@i, old_fit$beta0)
    if (!is.null(fit$gamma0)) {
      estim <- estim + partial_crossprod(
        matrix(1, nrow(inp.dat$y_train), 1),
        matrix(fit$gamma0, nrow = 1),
        inp.dat$y_valid@i, inp.dat$y_valid@p
      )
    }
    # add_to_cols_inplace_cpp(y_valid@x, y_valid@p, old_fit$gamma0)

    fit$error <- error_function(estim, inp.dat$y_valid@x)
    fit$lambda_beta <- lambda_beta
    fit$lambda_gamma <- lambda_gamma
    fit
  }

  results <- parallel_grid(grid, nllr.fit.step,
    "list",
    .packages = "IMR",
    .progress = TRUE,
    .seed = seed,
    # y_train = inp.dat$y_train,
    # y_valid = inp.dat$y_valid,
    # X = inp.dat$Xq,
    # Z = inp.dat$Zq,
    row_intercept = row_intercept,
    col_intercept = col_intercept,
    thresh = thresh,
    maxit = maxit,
    trace = verbose >= 2,
    seed = seed,
    error_function = error_function
  )


  # Select the best fit
  errors <- vapply(results, `[[`, numeric(1), "error")
  best_idx <- which.min(errors)
  best_fit <- results[[best_idx]]
  #--------------------
  # message >>
  if (verbose >= 1) {
    for (res in results) {
      message(sprintf(
        "<< lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f  >>",
        res$lambda_beta,
        sum(res$beta == 0) / length(res$beta),
        res$lambda_gamma,
        sum(res$gamma == 0) / length(res$gamma),
        res$error
      ))
    }
    message(sprintf(
      "<< Best fit >> lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f >>",
      best_fit$lambda_beta,
      sum(best_fit$beta == 0) / length(best_fit$beta),
      best_fit$lambda_gamma,
      sum(best_fit$gamma == 0) / length(best_fit$gamma),
      best_fit$error
    ))
    best_fit$init_hparams <- hpar
  }
  rm(results)

  nfit <- IMR::imr.fit_no_low_rank(
    Y = inp.dat$Y,
    X = inp.dat$Xq,
    Z = inp.dat$Zq,
    lambda_beta = best_fit$lambda_beta,
    lambda_gamma = best_fit$lambda_gamma,
    row_intercept = row_intercept,
    col_intercept = col_intercept,
    thresh = thresh,
    maxit = maxit,
    warm_start = best_fit,
    trace = FALSE
  )

  inp.dat$y_train <- IMR::as.Incomplete(nfit$resid * inp.dat$train_mask)
  inp.dat$y_valid <- IMR::as.Incomplete(nfit$resid * inp.dat$valid_mask)
  # inp.dat$y_train <- as(nfit$resid * (1 - inp.dat$valid_mask), "Incomplete")
  # inp.dat$y_valid <- as(nfit$resid * (inp.dat$valid_mask), "Incomplete")
  mfit <- IMR::imr.cv_M(inp.dat$y_train, inp.dat$y_valid,
    Y_full = nfit$resid, hpar = hpar,
    error_function = error_function, trace = verbose > 0,
    ls_initial = FALSE
  )

  nfit$u <- mfit$fit$u
  nfit$v <- mfit$fit$v
  nfit$d <- mfit$fit$d
  # best_fit$fit <- nfit
  best_fit$fit <- IMR::imr.fit(
    Y = inp.dat$Y,
    X = inp.dat$Xq,
    Z = inp.dat$Zq,
    r = mfit$rank_M,
    lambda_m = mfit$lambda_m,
    lambda_beta = best_fit$lambda_beta,
    lambda_gamma = best_fit$lambda_gamma,
    row_intercept = row_intercept,
    col_intercept = col_intercept,
    thresh = thresh,
    maxit = maxit,
    trace = verbose > 1,
    ls_initial = F,
    warm_start = nfit,
  )
  mfit$rank_M -> best_fit$rank_M
  mfit$lambda_m -> best_fit$lambda_m

  return(best_fit)
}
