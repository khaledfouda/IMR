#' @export
imr_tune_grid <- function(beta = c(0, NA, 20), # min, max, length
                          gamma = c(0, NA, 20),
                          laplace = c(0, NA, 10, 3), # min, max, length, streak
                          rank = c(2, 30, 2) # min, max, step
) {
  parse_param <- function(p, is_rank = FALSE, is_laplace = FALSE) {
    len <- length(p)
    stopifnot(len >= 1)
    stopifnot(!(is_rank && len < 2))
    if (is_rank) {
      list(
        min  = p[1],
        max  = if (len == 1) p[1] else p[2],
        step = if (len >= 3) p[3] else 2
      )
    } else {
      out <- list(
        min     = p[1],
        max     = if (len == 1) p[1] else if (is.na(p[2])) "auto" else p[2],
        length  = if (len == 1) 1 else if (len >= 3) p[3] else 20
      )
      if (is_laplace) {
        out$streaks <- if (len >= 5) p[5] else 1
      }
      out
    }
  }

  structure(
    list(
      beta    = parse_param(beta),
      gamma   = parse_param(gamma),
      rank    = parse_param(rank, is_rank = TRUE),
      laplace = parse_param(laplace, is_laplace = TRUE)
    ),
    class = "imr_tune_grid"
  )
}

#' @export
print.imr_tune_grid <- function(x, ...) {
  cat("\n== IMR Hyperparameter Configuration ==\n")

  fmt_range <- function(param) {
    if (param$length == 1) {
      return(paste0("Fixed at ", param$min))
    }

    # max_val <- if (identical(param$max, "auto")) "auto" else param$max

    sprintf(
      "Range: %s -> %s (Grid: %d points)",
      param$min, if (is.numeric(param$max)) round(param$max, 6) else param$max, param$length
    )
  }

  cat(sprintf("%-18s %s\n", "Beta:", fmt_range(x$beta)))
  cat(sprintf("%-18s %s\n", "Gamma:", fmt_range(x$gamma)))
  steps_str <- paste(x$laplace$steps, collapse = ", ")
  cat(sprintf(
    "%-18s Range: %s -> %s (Length: %s, Streaks: %s)\n",
    "Laplace:",
    x$laplace$min,
    if (is.numeric(x$laplace$max)) round(x$laplace$max, 6) else x$laplace$max,
    x$laplace$length,
    x$laplace$streaks
  ))
  cat(sprintf(
    "%-18s Range: %d -> %d (Step: %s)\n",
    "Rank:",
    x$rank$min,
    x$rank$max,
    x$rank$step
  ))

  cat("===========================================================\n\n")

  invisible(x)
}

#-----------------------------------------------------
#' @export
imr_set_grid_limits <- function(data,
                                grid,
                                default_rank = 2,
                                default_lambda_m = 0,
                                default_lambda_beta = 0,
                                default_lambda_gamma = 0,
                                convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE),
                                bisection_iter = 15,
                                verbose = 0) {
  if (data$meta$has_X && grid$beta$max == "auto") {
    grid$beta$max <- IMR::get_lambda_lasso_max(data, "beta",
      rank = default_rank,
      lambda_m = default_lambda_m,
      convergence = convergence,
      bisection_iter = bisection_iter,
      verbose = verbose
    )
  }
  if (data$meta$has_Z && grid$gamma$max == "auto") {
    grid$gamma$max <- IMR::get_lambda_lasso_max(data, "gamma",
      rank = default_rank,
      lambda_m = default_lambda_m,
      convergence = convergence,
      bisection_iter = bisection_iter,
      verbose = verbose
    )
  }
  if (grid$laplace$max == "auto") {
    grid$laplace$max <- IMR::get_lambda_m_max(data,
      rank = default_rank,
      lambda_beta = default_lambda_beta,
      lambda_gamma = default_lambda_gamma,
      convergence = convergence,
      bisection_iter = bisection_iter,
      verbose = verbose
    )
  }
  return(grid)
}

#------------------------------------------------------
#' @export
get_lambda_m_max <-
  function(data,
           lambda_beta = 0,
           lambda_gamma = 0,
           rank = 2,
           bisection_iter = 15,
           convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE),
           verbose = 0) {
    need_fit <- any(!is.null(data$Xq), !is.null(data$Zq), data$model$intercept_row, data$model$intercept_col)


    if (!need_fit) {
      lambda_kkt <- IMR::svd_opt(data$Y, 1)$d[1]
    } else {
      fit <- IMR::imr_fit(data,
        rank = 0,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        convergence = convergence
      )
      # return largest singular value
      lambda_kkt <- IMR::svd_opt(fit$residuals, 1)$d[1]
    }
    lower <- 0
    upper <- lambda_kkt
    baseline_fit <- NULL

    # helper, fits a single model
    fit_test <- function(lam) {
      IMR::imr_fit(data,
        rank = rank,
        lambda_m = lam,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        convergence = convergence,
        warm_start = baseline_fit
      )
    }
    baseline_fit <- fit_test(0)

    # we begin by adjusting the upperbound (in case the KKT bound isn't enough)
    for (i in seq_len(bisection_iter)) {
      coefs <- fit_test(upper)$coefficients$d
      zero_ratio <- mean(abs(coefs) < 1e-6)

      if (zero_ratio == 1) {
        # we achieved an upper-bound. Let's return.
        break
      } else {
        # It is not fully sparse. We must try a HIGHER lambda.
        upper <- lambda_kkt * (i + 1)
      }
      if (verbose >= 2) {
        message(sprintf("lambda = %.4f, zero ratio = %.2f", upper, zero_ratio))
      }
    }

    for (i in seq_len(bisection_iter)) {
      mid <- (lower + upper) / 2
      coefs <- fit_test(mid)$coefficients$d
      zero_ratio <- mean(abs(coefs) < 1e-6)

      if (zero_ratio == 1) {
        # It is fully sparse. We can try a lower lambda.
        upper <- mid
      } else {
        # It is NOT fully sparse. We must try a higher lambda.
        lower <- mid
      }
      if (verbose >= 2) {
        message(sprintf("lambda = %.4f, zero ratio = %.2f", mid, zero_ratio))
      }
    }

    # The upper bound is guaranteed to be the minimum lambda that achieves 100% sparsity
    lambda_sup <- upper

    if (verbose > 0) {
      message(sprintf(
        "Target: Laplace | KKT max: %.3f | Empiric Sup: %.3f (%.1f%% of KKT)",
        lambda_kkt, lambda_sup, 100 * lambda_sup / lambda_kkt
      ))
    }

    return(lambda_sup)
  }
#------------------------------------------------
#' Find the minimum Lasso lambda that forces all covariates to zero
#' @export
get_lambda_lasso_max <- function(
  data, # Must be an 'imr_data' S3 object
  target = c("beta", "gamma"),
  rank = 2,
  lambda_m = 0,
  bisection_iter = 15,
  convergence = imr_convergence(),
  verbose = 0
) {
  stopifnot(inherits(data, "imr_data"))
  target <- match.arg(target, c("beta", "gamma"))
  is_beta <- target == "beta"
  # remove the other set of covariates from the data.
  if (is_beta) data$Zq <- NULL else data$Xq <- NULL

  if (is_beta && !data$meta$has_X) stop("Target is 'beta' but no X matrix found in data.")
  if (!is_beta && !data$meta$has_Z) stop("Target is 'gamma' but no Z matrix found in data.")


  # ---  Baseline Fit  ---
  # Fit the model but keep our target penalized to 0
  baseline_fit <- imr_fit(
    data = data,
    rank = rank,
    lambda_m = lambda_m,
    lambda_beta = 0,
    lambda_gamma = 0,
    convergence = convergence
  )

  # --- Calculate training Residuals for KKT Conditions ---
  resids <- baseline_fit$residuals

  # ---  KKT Max  ---
  # The theoretical minimum lambda that forces all coefs to 0 is max(abs(X^T * Residuals))
  if (is_beta) {
    lambda_kkt <- max(abs(crossprod(data$Xq, resids)))
  } else {
    lambda_kkt <- max(abs(resids %*% data$Zq))
  }

  # --- Bisection Search for Supremum ---
  # Supermum is typicaly smaller than the KKT max due to missing data.

  lower <- 0
  upper <- lambda_kkt


  # Helper to fit a single model
  fit_test <- function(lam) {
    imr_fit(
      data = data,
      rank = rank,
      lambda_m = lambda_m,
      lambda_beta = if (is_beta) lam else 0,
      lambda_gamma = if (!is_beta) lam else 0,
      warm_start = baseline_fit,
      convergence = convergence
    )
  }
  # we begin by adjusting the upperbound (in case the KKT bound isn't enough)
  for (i in seq_len(bisection_iter)) {
    test_model <- fit_test(upper)
    # Check if all coefficients are effectively zero
    coefs <- if (is_beta) test_model$coefficients$beta else test_model$coefficients$gamma
    zero_ratio <- mean(abs(coefs) < 1e-6)

    if (zero_ratio == 1) {
      # we achieved an upper-bound. Let's return.
      break
    } else {
      # It is not fully sparse. We must try a HIGHER lambda.
      upper <- lambda_kkt * (i + 1)
    }
    if (verbose >= 2) {
      message(sprintf("lambda = %.4f, zero ratio = %.2f", upper, zero_ratio))
    }
  }

  for (i in seq_len(bisection_iter)) {
    mid <- (lower + upper) / 2
    test_model <- fit_test(mid)

    # Check if all coefficients are effectively zero
    coefs <- if (is_beta) test_model$coefficients$beta else test_model$coefficients$gamma
    zero_ratio <- mean(abs(coefs) < 1e-6)

    if (zero_ratio == 1) {
      # It is fully sparse. We can try a lower lambda.
      upper <- mid
    } else {
      # It is NOT fully sparse. We must try a higher lambda.
      lower <- mid
    }
    if (verbose >= 2) {
      message(sprintf("lambda = %.4f, zero ratio = %.2f", mid, zero_ratio))
    }
  }

  # The upper bound is guaranteed to be the minimum lambda that achieves 100% sparsity
  lambda_sup <- upper

  if (verbose > 0) {
    message(sprintf(
      "Target: %s | KKT max: %.3f | Empiric Sup: %.3f (%.1f%% of KKT)",
      ifelse(is_beta, "Beta", "Gamma"), lambda_kkt, lambda_sup, 100 * lambda_sup / lambda_kkt
    ))
  }

  return(lambda_sup)
}

#-----------------------------------------------------------------------
