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

# -----------------------------------------------------
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


#' Dual norm of the loss gradient at the null solution
#'
#'  value of the penalty at which the target component is zero.
#' For convex problems (and this is a convex problem), the KKT value is the supremum
#'
#' @param resid An `Incomplete` residual matrix from the model fitted WITHOUT
#'   the target component, all other parameters at their conditional optimum.
#' @param target One of "nuclear", "beta", "gamma".
#' @param data The `imr_data` object (for `Xq`, `Zq`, and the shared flags).
#' @param huber_c Transition point of the Huber loss, or `Inf` for least
#'   squares. The gradient involves psi_c(r), not r.
#' @noRd
imr_lambda_max_kkt <- function(resid, target, data, huber_c = Inf) {

  R <- resid
  if (is.finite(huber_c)) R@x <- pmin(pmax(R@x, -huber_c), huber_c)

  switch(
    target,
    # grad = -P_Omega*(psi_c(r)); dual of the nuclear norm is the operator norm
    nuclear = svd_opt(R, 1)$d[1],

    # Shared effects sum the gradient across columns/rows
    beta = if (isTRUE(data$model$shared_beta)) {
      max(abs(crossprod(data$Xq, Matrix::rowSums(R))))
    } else {
      max(abs(crossprod(data$Xq, R)))
    },

    gamma = if (isTRUE(data$model$shared_gamma)) {
      max(abs(Matrix::colSums(R) %*% data$Zq))
    } else {
      max(abs(R %*% data$Zq))
    },

    stop("unknown target")
  )
}











#' Penalty value at which a model is exactly zero
#'
#' Computes the KKT bound, then optionally spends a small number of
#' model fits verifying it and refining downward on the log scale. Refinement
#' is off by default.
#'
#' @param verify_iter Integer. Fits used to confirm the bound does zero the
#'   target, expanding geometrically if not. `0` trusts the KKT bound.
#' @param refine_iter Integer. Log-scale bisection steps used to shrink the
#'   bound afterwards. `0` (the default) keeps the exact value.
#' @param zero_tol Relative tolerance for "all coefficients zero", measured
#'   against the un-penalized fit.
#' @noRd
imr_lambda_max <- function(data,
                           target = c("nuclear", "beta", "gamma"),
                           rank = 2,
                           lambda_m = 0,
                           lambda_beta = 0,
                           lambda_gamma = 0,
                           huber_shift = 0,
                           convergence = imr_convergence(trace = FALSE),
                           verify_iter = 1L,
                           refine_iter = 0L,
                           zero_tol = 1e-4,
                           verbose = 0) {

  target <- match.arg(target)
  training <- isTRUE(data$meta$split_data)

  # --- Null model: turn off the target ---------
  null_data <- data
  switch(target,
         nuclear = null_data$model$low_rank_component <- FALSE,
         beta    = null_data$model$row_covariates     <- FALSE,
         gamma   = null_data$model$col_covariates     <- FALSE
  )
  null_fit <- imr_fit(null_data,
                      rank = if (target == "nuclear") 0 else rank,
                      lambda_m = lambda_m, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma,
                      huber_shift = huber_shift, convergence = convergence, training = training
  )

  huber_c <- if (huber_shift > 0) null_fit$meta$huber$huber_c else Inf
  lambda_kkt <- imr_lambda_max_kkt(null_fit$residuals, target, data, huber_c)

  if (verify_iter <= 0L && refine_iter <= 0L) return(lambda_kkt)

  # ---  in case of fine-tuning ----------------------------------
  extract <- function(f) switch(target,
                                nuclear = f$coefficients$d,
                                beta    = f$coefficients$beta,
                                gamma   = f$coefficients$gamma
  )
  fit_at <- function(lam, warm) imr_fit(data,
                                        rank = rank,
                                        lambda_m     = if (target == "nuclear") lam else lambda_m,
                                        lambda_beta  = if (target == "beta")    lam else lambda_beta,
                                        lambda_gamma = if (target == "gamma")   lam else lambda_gamma,
                                        huber_shift = huber_shift, convergence = convergence,
                                        training = training, warm_start = warm
  )

  # Dense fit helps refine the zero_tol. The reason is in case the values are too large
  # then the zero tol shouldn't be too small, to save time
  dense_fit <- fit_at(0, NULL)
  scale_ref <- max(abs(extract(dense_fit)), na.rm = TRUE)
  # if it's already zero, no reason to continue. We also guard against NA values, but
  # that shouldn't happen so the first condition should never be true.
  if (!is.finite(scale_ref) || scale_ref <= 0) return(lambda_kkt)

  # --- see if the limit is low, then test larger values (doubling the value at each iteration)
  is_zero <- function(lam) max(abs(extract(fit_at(lam, dense_fit)))) < zero_tol * scale_ref
  upper <- lambda_kkt
  lower <- 0
  ok <- FALSE
  for (i in seq_len(max(verify_iter, 1L))) {
    if (is_zero(upper)) { ok <- TRUE; break }
    lower <- upper
    upper <- upper * 2
    if (verbose >= 2) message(sprintf("  expand: lambda = %.5g", upper))
  }
  # in case we never arrive at that value
  if (!ok && verbose > 0)
    warning(sprintf(
      "%s: the target is not zero at lambda = %.5g. try a smaller `thresh` or a larger `maxit`.",
      target, upper), call. = FALSE)

  # --- Optional refinement on the LOG scale ---------------------------
  # this is mid-point search.
  # The grid is exp(seq(log(max), log(min), ...))
  if (refine_iter > 0L && ok) {
    lo <- log(max(lower, upper * 1e-6))
    hi <- log(upper)
    for (i in seq_len(refine_iter)) {
      mid <- (lo + hi) / 2
      if (is_zero(exp(mid))) hi <- mid else lo <- mid
      if (verbose >= 2) message(sprintf("  refine: lambda = %.5g", exp(mid)))
    }
    upper <- exp(hi)
  }

  if (verbose > 0)
    message(sprintf("Target: %-7s | KKT: %.5g | returned: %.5g (%.1f%% of KKT)%s",
                    target, lambda_kkt, upper, 100 * upper / lambda_kkt,
                    if (huber_shift > 0) sprintf(" | Huber c = %.4g", huber_c) else ""))
  upper
}







