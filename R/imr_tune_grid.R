#' @title Define Hyperparameter Tuning Grid for IMR
#'
#' @description
#' Constructs a configuration object specifying the search space for hyperparameter
#' optimization within an Incomplete Matrix Regression (IMR) framework. This object
#' facilitates the execution of the `imr_tune` function by defining the domain for identifying
#'  optimal hyperparameters.
#'
#' @param beta A numeric vector defining the search grid for the row-wise covariate
#'   regularization parameter (\eqn{\lambda_\beta}). Expected format: `c(min, max, length)`.
#'   Defaults to `c(0, NA, 20)`.
#' @param gamma A numeric vector defining the search grid for the column-wise covariate
#'   regularization parameter (\eqn{\lambda_\Gamma}). Expected format: `c(min, max, length)`.
#'   Defaults to `c(0, NA, 20)`.
#' @param nuclear A numeric vector defining the search grid for the nuclear norm
#'   penalty (\eqn{\lambda_M}) applied to the low-rank component. Expected format:
#'   `c(min, max, length, streaks)`, where `streaks` specifies the patience threshold
#'   for early stopping. Defaults to `c(0, NA, 20, 2)`.
#' @param rank A numeric vector defining the search grid for the rank of the
#'   low-rank component. Expected format: `c(min, max, step, streaks)`, where
#'   `streaks` specifies the patience threshold for early stopping.
#'   Defaults to `c(2, 30, 2, 2)`.
#'
#' @details
#' Specifications for `beta`, `gamma`, and `nuclear` parameters are subject to the
#' following conventions:
#' * If the second element (`max`) is
#'     specified as `NA`, the upper bound of the grid is determined internally
#'     via `imr_set_grid_limits`. You need to call that function before calling `imr_tune`.
#' *  Providing a vector of length 1 (e.g., `beta = 0.5`)
#'     constrains the parameter to that value throughout the tuning process.
#' *  If a vector of length 2 is provided, the elements are
#'     interpreted as the minimum and maximum bounds, respectively, with remaining
#'     grid attributes reverting to their default values.
#'
#' @return An object of class `"imr_tune_grid"`, a list of the parameters and their values.
#'
#' @export
#' @seealso [print.imr_tune_grid()], [imr_set_grid_limits()]
#' @examples
#' # Initialize default grid with automated limit calculation
#' default_grid <- imr_tune_grid()
#'
#' # Fix beta at 0 and define custom search spaces for gamma, nuclear, and rank
#' custom_grid <- imr_tune_grid(
#'   beta = 0,
#'   gamma = c(0.01, 1, 10),
#'   nuclear = c(0, NA, 15, 3), # 15 points, patience threshold of 3
#'   rank = c(2, 10, 2, 2)
#' )
#'
#' # print the grid's information
#' print(custom_grid)
imr_tune_grid <- function(beta = c(0, NA, 20), # min, max, length
                          gamma = c(0, NA, 20),
                          nuclear = c(0, NA, 20, 2), # min, max, length, streak
                          rank = c(2, 30, 2, 2) # min, max, step, streak
) {
  parse_param <- function(p, is_rank = FALSE, is_nuclear = FALSE) {
    len <- length(p)
    stopifnot(len >= 1)
    #stopifnot(!(is_rank && len < 2))
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
      if (is_nuclear) {
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
      nuclear = parse_param(nuclear, is_nuclear = TRUE)
    ),
    class = "imr_tune_grid"
  )
}

#' @title Summary of the hyperparameter grid
#' @param x An `imr_tune_grid` object
#' @param ... Additional arguments to comply with generic function
#' @return The input `imr_tune_grid` object `x`, invisibly. Called for its side
#'   effect of printing the configured search ranges to the console.
#' @seealso [imr_tune_grid()]
#' @inherit imr_tune_grid examples
#' @export
#' @method print imr_tune_grid
print.imr_tune_grid <- function(x, ...) {
  cat("\n== IMR Hyperparameter Configuration ==\n")

  fmt_range <- function(param, nuclear=FALSE, rank = FALSE) {
    if (param$min == param$max) {
      return(paste0("Fixed at ", param$min))
    }

    # max_val <- if (identical(param$max, "auto")) "auto" else param$max

    if(nuclear){
      text_partial <- sprintf("(Length: %s, Streaks: %s)", param$length, param$streaks)
    }else if(rank){
      text_partial <- sprintf("(Step: %s, Streaks: %d)", param$step, param$streaks)
    }else
      text_partial <- sprintf("(Grid: %d points)", param$length)

    paste(
    sprintf(
      "Range: %s -> %s ",
      param$min, if (is.numeric(param$max)) round(param$max, 6) else param$max
    ),
    text_partial)
  }

  cat(sprintf("%-18s %s\n", "Beta:", fmt_range(x$beta)))
  cat(sprintf("%-18s %s\n", "Gamma:", fmt_range(x$gamma)))
  cat(sprintf("%-18s %s\n", "Nuclear:", fmt_range(x$nuclear, nuclear=TRUE)))
  cat(sprintf("%-18s %s\n", "Rank:", fmt_range(x$rank, rank=TRUE)))

  cat("===========================================================\n\n")

  invisible(x)
}

#-----------------------------------------------------
#' @title Automatically Determine Hyperparameter Grid Maximum Values
#'
#' @description
#' Computes the optimal upper bounds for the hyperparameter tuning
#' grid when specifications are set to `"auto"`. The function identifies the
#' minimal regularization parameter required to shrink all corresponding
#' coefficients to zero.
#'
#' @param data An object of class `"imr_data"` containing the response matrices
#'   and covariate structures.
#' @param grid An object of class `"imr_tune_grid"`, initialized via
#'   `imr_tune_grid()`.
#' @param default_rank,default_lambda_m Integer,Numeric. The fixed rank and low-rank component
#'   penalty used as reference when
#'   estimating the maximum \eqn{\lambda_\beta} or \eqn{\lambda_\Gamma}.
#'   Default to `2` and `0`.
#' @param default_lambda_beta,default_lambda_gamma Numeric. The row and column
#'   covariate penalties used as
#'   a reference when estimating the maximum \eqn{\lambda_M}. Defaults to `0`.
#' @param convergence An `"imr_convergence"` object specifying the numerical
#'   tolerances and optimization parameters for internal model fits. Defaults
#'   to `imr_convergence(trace = FALSE, ls_initial = FALSE)`.
#' @param bisection_iter Integer. The number of iterations for the bisection
#'   algorithm employed to refine the estimated upper bound. Defaults to `5`.
#' @param verbose Integer. Level of diagnostic output. `0` (default) is silent,
#'   `1` reports final parameters, and `2` provides per-iteration updates.
#'
#' @details
#' The procedure isolates each hyperparameter to determine its saturation
#' point through a two-stage estimation process:
#' \enumerate{
#'   \item An initial theoretical upper bound is
#'     derived based on the Karush-Kuhn-Tucker (KKT)  conditions.
#'   \item The function executes a bisection search over
#'     `bisection_iter` iterations using KKT estimates as initial values. This identifies the minimum
#'     of the penalty values that result in a zero-solution for the
#'     targeted parameters.
#' }
#'
#' During the estimation of the maximum \eqn{\lambda_\beta} or \eqn{\lambda_\Gamma},
#' the other covariate penalty is treated as infinite, while the
#' low-rank component is set to `default_rank` and `default_lambda_m`.
#' Conversely, estimation of the maximum \eqn{\lambda_M} (nuclear norm penalty)
#' assumes the covariate penalties are fixed at `default_lambda_beta` and
#' `default_lambda_gamma`.
#'
#' @return A modified `"imr_tune_grid"` object where all `"auto"` placeholders
#'   are replaced by the numerically determined maximum values.
#'
#' @seealso [imr_tune_grid()], [imr_tune()], [imr_data()]
#' @inherit imr_tune examples
#'
#' @export
imr_set_grid_limits <- function(data,
                                grid,
                                default_rank = 2,
                                default_lambda_m = 0,
                                default_lambda_beta = 0,
                                default_lambda_gamma = 0,
                                convergence = imr_convergence(trace = FALSE, ls_initial = FALSE),
                                bisection_iter = 5,
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
  if (data$model$low_rank_component && grid$nuclear$max == "auto") {
    grid$nuclear$max <- imr_get_lambda_m_max(data,
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
#' @noRd
imr_get_lambda_m_max <-
  function(data,
           lambda_beta = 0,
           lambda_gamma = 0,
           rank = 2,
           bisection_iter = 15,
           convergence = imr_convergence(trace = FALSE, ls_initial = FALSE),
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
      lambda_kkt <- svd_opt(fit$residuals, 1)$d[1]
    }
    lower <- 0
    upper <- lambda_kkt

    # helper, fits a single model
    fit_test <- function(lam) {
      imr_fit(data,
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

    # we begin by adjusting the upper bound (in case the KKT bound isn't enough)
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
        message(sprintf("parameter = %.4f, zero ratio = %.2f", upper, zero_ratio))
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
        message(sprintf("parameter = %.4f, zero ratio = %.2f", mid, zero_ratio))
      }
    }

    # The upper bound is guaranteed to be the minimum lambda that achieves 100% sparsity
    lambda_sup <- upper

    if (verbose > 0) {
      message(sprintf(
        "Target: nuclear | KKT max: %.3f | Empiric Sup: %.3f (%.1f%% of KKT)",
        lambda_kkt, lambda_sup, 100 * lambda_sup / lambda_kkt
      ))
    }

    return(lambda_sup)
  }
#------------------------------------------------
#' Find the minimum Lasso lambda that forces all covariates to zero
#' @noRd
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
  # we begin by adjusting the upper bound (in case the KKT bound isn't enough)
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
      message(sprintf("parameter = %.4f, zero ratio = %.2f", upper, zero_ratio))
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
      message(sprintf("parameter = %.4f, zero ratio = %.2f", mid, zero_ratio))
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
