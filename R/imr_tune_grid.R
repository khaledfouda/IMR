#' Define Hyperparameter Tuning Grid for IMR
#'
#' @description
#' Creates a configuration object defining the search space for hyperparameter
#' tuning in an IMR model. This object
#' is used by the `imr_tune` function to tuning and finding the optimal hyperparamters.
#'
#' @param beta Numeric vector defining the grid for the row covariate penalty
#'   (`lambda_beta`). Format: `c(min, max, length)`. Defaults to `c(0, NA, 20)`.
#' @param gamma Numeric vector defining the grid for the column covariate penalty
#'   (`lambda_gamma`). Format: `c(min, max, length)`. Defaults to `c(0, NA, 20)`.
#' @param laplace Numeric vector defining the grid for the low-rank component
#'   penalty (`lambda_m`). Format: `c(min, max, length, streaks)`. The `streaks`
#'   parameter dictates the patience for early stopping. Defaults to `c(0, NA, 20, 2)`.
#' @param rank Numeric vector defining the grid for the rank of the low-rank
#'   component. Format: `c(min, max, step, streaks)`. The `streaks` parameter
#'   dictates the patience for early stopping. Defaults to `c(2, 30, 2, 2)`.
#'
#' @details
#' Regarding the inputs for `beta`, `gamma`, and `laplace`:
#' * **Automatic Maximum:** If the second element (`max`) is `NA`, the upper limit
#'     can be automatically determined using the function `imr_set_grid_limits`.
#' * **Fixed Value:** If a vector of length 1 is provided (e.g., `beta = 0.5`),
#'     the parameter is held constant during tuning, bypassing the grid search
#'     for that specific hyperparameter. This also holds for the `rank`.
#' Generally, if the input vector is of length 1, then this value is used as a fixed
#'    value. If the input vector is of length 2, then they are used as the minimum and
#'    maximum and values and the rest are set to the default values.
#'
#' @return An object of class `"imr_tune_grid"`, which is a list containing the
#'   parsed grid boundaries, lengths, step sizes, and early stopping criteria
#'   for each parameter.
#'
#' @export
#'
#' @examples
#' # Default grid (max limits to be calculated automatically)
#' default_grid <- imr_tune_grid()
#'
#' # Fix beta to 0, tune gamma and laplace, and ranks 2 through 10
#' custom_grid <- imr_tune_grid(
#'   beta = 0,
#'   gamma = c(0.01, 1, 10),
#'   laplace = c(0, NA, 15, 3), # 15 points, patience of 3
#'   rank = c(2, 10, 2, 2)
#' )
imr_tune_grid <- function(beta = c(0, NA, 20), # min, max, length
                          gamma = c(0, NA, 20),
                          laplace = c(0, NA, 20, 2), # min, max, length, streak
                          rank = c(2, 30, 2, 2) # min, max, step, streak
) {
  parse_param <- function(p, is_rank = FALSE, is_laplace = FALSE) {
    len <- length(p)
    stopifnot(len >= 1)
    stopifnot(!(is_rank && len < 2))
    if (is_rank) {
      list(
        min = p[1],
        max = if (len == 1) p[1] else p[2],
        step = if (len >= 3) p[3] else 2,
        streaks = if (len >= 4) p[4] else 2
      )
    } else {
      out <- list(
        min     = p[1],
        max     = if (len == 1) p[1] else if (is.na(p[2])) "auto" else p[2],
        length  = if (len == 1) 1 else if (len >= 3) p[3] else 20
      )
      if (is_laplace) {
        out$streaks <- if (len >= 4) p[4] else 2
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
    "%-18s Range: %d -> %d   (Step: %s, Streaks: %d)\n",
    "Rank:",
    x$rank$min,
    x$rank$max,
    x$rank$step,
    x$rank$streaks
  ))

  cat("===========================================================\n\n")

  invisible(x)
}

#-----------------------------------------------------
#' Automatically Determine Hyperparameter Grid Limits
#'
#' @description
#' Calculates the appropriate upper limits (`max`) for the hyperparameter tuning
#' grid when they are specified as `"auto"`. It finds the smallest regularization
#' parameter that forces all corresponding coefficients to zero.
#'
#' @param data An object of class `"imr_data"` containing the model structure
#'   and matrices.
#' @param grid An object of class `"imr_tune_grid"`, typically created by
#'   `imr_tune_grid()`.
#' @param default_rank Integer. The default rank to use when searching for the
#'   maximum `lambda_beta` or `lambda_gamma`. Defaults to `2`.
#' @param default_lambda_m Numeric. The default penalty for the low-rank component
#'   when searching for covariate penalty limits. Defaults to `0`.
#' @param default_lambda_beta Numeric. The default penalty for row covariates
#'   when searching for the maximum `lambda_m`. Defaults to `0`.
#' @param default_lambda_gamma Numeric. The default penalty for column covariates
#'   when searching for the maximum `lambda_m`. Defaults to `0`.
#' @param convergence An `"imr_convergence"` object controlling the internal model
#'   fits during the search. Defaults to `IMR::imr_convergence(trace = FALSE, ls_initial = FALSE)`.
#' @param bisection_iter Integer. The number of iterations to use in the bisection
#'   algorithm to refine the upper limit. Defaults to `15`.
#' @param verbose Integer. Controls the level of printed output during the search.
#'   Defaults to `0` (silent). Set it `1` to output the final parameters and `2` for
#'   an output at every iteration.
#'
#' @details
#' The function isolates each hyperparameter to find its maximum using a
#' two-step process:
#' \enumerate{
#'   \item **KKT Initialization:** It calculates an initial theoretical upper limit
#'     based on the Karush-Kuhn-Tucker (KKT) conditions.
#'   \item **Refinement:** Because the theoretical KKT value can sometimes
#'     be slightly over or underestimated in practice, the function runs a bisection
#'     search for `bisection_iter` iterations to empirically find the exact smallest
#'     penalty value that shrinks all targeted parameters to zero.
#' }
#'
#' When searching for the maximum `lambda_beta` or `lambda_gamma`, the other
#' covariate penalty is effectively set to infinity (ignored), while the low-rank
#' component is fixed at `default_rank` and `default_lambda_m`. Conversely, when
#' searching for the maximum `lambda_m` (laplace), the covariate penalties are fixed
#' at `default_lambda_beta` and `default_lambda_gamma`.
#'
#' @return The modified `"imr_tune_grid"` object with all `"auto"` maximums
#'   replaced by their calculated numeric values.
#'
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
  if (data$model$row_covariates && data$meta$has_X && grid$beta$max == "auto") {
    grid$beta$max <- imr_get_lambda_lasso_max(data, "beta",
      rank = default_rank,
      lambda_m = default_lambda_m,
      convergence = convergence,
      bisection_iter = bisection_iter,
      verbose = verbose
    )
  }
  if (data$model$col_covariates && data$meta$has_Z && grid$gamma$max == "auto") {
    grid$gamma$max <- imr_get_lambda_lasso_max(data, "gamma",
      rank = default_rank,
      lambda_m = default_lambda_m,
      convergence = convergence,
      bisection_iter = bisection_iter,
      verbose = verbose
    )
  }
  if (data$model$low_rank_component && grid$laplace$max == "auto") {
    grid$laplace$max <- imr_get_lambda_m_max(data,
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
imr_get_lambda_m_max <-
  function(data,
           lambda_beta = 0,
           lambda_gamma = 0,
           rank = 2,
           bisection_iter = 15,
           convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE),
           verbose = 0) {
    need_fit <- any(!is.null(data$Xq), !is.null(data$Zq), data$model$row_intercept, data$model$col_intercept)


    if (!need_fit) {
      lambda_kkt <- svd_opt(data$Y, 1)$d[1]
    } else {
      data$model$low_rank_component <- FALSE
      fit <- imr_fit(data,
        rank = 0,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        convergence = convergence
      )
      data$model$low_rank_component <- TRUE
      # return largest singular value
      lambda_kkt <- IMR::svd_opt(fit$residuals, 1)$d[1]
    }
    lower <- 0
    upper <- lambda_kkt

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
    baseline_fit <- NULL
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
        # upper <- lambda_kkt * (i + 1)
        lower <- upper
        upper <- upper * 1.5
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
imr_get_lambda_lasso_max <- function(
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
  if (is_beta) {
    data$model$row_covariates <- FALSE
  } else {
    data$model$col_covariates <- FALSE
  }

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

  #-- put them back in the model
  if (is_beta) {
    data$model$row_covariates <- TRUE
  } else {
    data$model$col_covariates <- TRUE
  }

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
  baseline_fit <- NULL
  baseline_fit <- fit_test(0)
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
      # upper <- lambda_kkt * (i + 1)
      lower <- upper
      upper <- upper * 1.5
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
