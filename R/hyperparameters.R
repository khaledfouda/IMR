#' @export
imr_tune_grid <- function(beta = c(0, NA, 20, 0), # min, max, length, default
                            gamma = c(0, NA, 20, 0),
                            laplace = c(0, NA, 10, 0, 3), # min, max, length, default, streak
                            rank = c(2, 30, 2, 2) # min, max, step, default
                            ) {

  parse_param <- function(p, is_rank = FALSE, is_laplace = FALSE) {
    len <- length(p)
    stopifnot(len >= 1)
    stopifnot(!(is_rank && len < 2))
    if (is_rank) {
      list(
        min  = p[1],
        max  = if(len == 1) p[1] else p[2],
        step = if (len >= 3) p[3] else 2,
        default = if(len == 1) p[1] else if (len >= 4) p[4] else 2
      )
    } else {
      out = list(
        min     = p[1],
        max     = if (len == 1) p[1] else if(is.na(p[2])) "auto" else p[2],
        length  = if (len == 1) 1 else if (len >= 3) p[3] else 20,
        default = if (len == 1) p[1] else if (len >= 4) p[4] else 0
      )
      if(is_laplace)
        out$streaks = if(len >= 5) p[5] else 1
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

    #max_val <- if (identical(param$max, "auto")) "auto" else param$max

    sprintf(
      "Range: %s -> %s (Grid: %d points) | Default: %s",
      param$min, if(is.numeric(param$max)) round(param$max,6) else param$max, param$length, param$default
    )
  }

  cat(sprintf("%-18s %s\n", "Beta:", fmt_range(x$beta)))
  cat(sprintf("%-18s %s\n", "Gamma:", fmt_range(x$gamma)))
  steps_str <- paste(x$laplace$steps, collapse = ", ")
  cat(sprintf(
    "%-18s Range: %s -> %s (Length: %s, Streaks: %s) | Default: %s\n",
    "Laplace:",
    x$laplace$min,
    if(is.numeric(x$laplace$max)) round(x$laplace$max,6) else x$laplace$max,
    x$laplace$length,
    x$laplace$streaks,
    x$laplace$default
  ))
  cat(sprintf(
    "%-18s Range: %d -> %d (Step: %s) | Default: %s\n",
    "Rank:",
    x$rank$min,
    x$rank$max,
    x$rank$step,
    x$rank$default
  ))

  cat("===========================================================\n\n")

  invisible(x)
}

#-----------------------------------------------------
#' @export
imr_set_grid_limits <- function(data,
                                grid,
                                intercept_row = FALSE,
                                intercept_col = FALSE,
                                shared_beta = FALSE,
                                shared_gamma = FALSE,
                                convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE),
                                bisection_iter = 15,
                                verbose = 0){

  if(data$meta$has_X && grid$beta$max == "auto")
    grid$beta$max <- IMR::get_lambda_lasso_max(data, "beta",
                                               rank = grid$rank$default,
                                               lambda_m = grid$laplace$default,
                                               intercept_row = intercept_row,
                                               intercept_col = intercept_col,
                                               shared_effects = shared_beta,
                                               convergence = convergence,
                                               bisection_iter = bisection_iter,
                                               verbose = verbose)
  if(data$meta$has_Z && grid$gamma$max == "auto")
    grid$gamma$max <- IMR::get_lambda_lasso_max(data, "gamma",
                                               rank = grid$rank$default,
                                               lambda_m = grid$laplace$default,
                                               intercept_row = intercept_row,
                                               intercept_col = intercept_col,
                                               shared_effects = shared_gamma,
                                               convergence = convergence,
                                               bisection_iter = bisection_iter,
                                               verbose = verbose)
  if(grid$laplace$max == "auto")
    grid$laplace$max <- IMR::get_lambda_m_max(data,
                                              intercept_row = intercept_row,
                                              intercept_col = intercept_col,
                                              rank = grid$rank$default,
                                              lambda_beta = grid$beta$default,
                                              lambda_gamma = grid$gamma$default,
                                              shared_beta = shared_beta,
                                              shared_gamma = shared_gamma,
                                              convergence = convergence,
                                              bisection_iter = bisection_iter,
                                              verbose = verbose)
  return(grid)
}

#------------------------------------------------------
#' @export
get_lambda_m_max <-
  function(data,
           intercept_row = FALSE,
           intercept_col = FALSE,
           lambda_beta = 0,
           lambda_gamma = 0,
           rank        = 2,
           shared_beta = FALSE,
           shared_gamma = FALSE,
           bisection_iter = 15,
           convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE),
           verbose = 0) {
    need_fit <- any(!is.null(data$Xq), !is.null(data$Zq), intercept_row, intercept_col)


    if (!need_fit) {
      lambda_kkt = IMR::svd_opt(data$Y, 1)$d[1]
    }else{
      fit <- IMR::imr_fit(data,
        rank = 0,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        intercept_row = intercept_row,
        intercept_col = intercept_col,
        shared_beta = shared_beta,
        shared_gamma = shared_gamma,
        convergence = convergence
      )
      # return largest singular value
      lambda_kkt =IMR::svd_opt(fit$residuals, 1)$d[1]
    }
    lower <- 0
    upper <- lambda_kkt
    baseline_fit <- NULL

    #helper, fits a single model
    fit_test <- function(lam){
      IMR::imr_fit(data,
                   rank = rank,
                   lambda_m =  lam,
                   lambda_beta = lambda_beta,
                   lambda_gamma = lambda_gamma,
                   intercept_row = intercept_row,
                   intercept_col = intercept_col,
                   shared_beta = shared_beta,
                   shared_gamma = shared_gamma,
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
  intercept_row = TRUE,
  intercept_col = TRUE,
  shared_effects = FALSE, # Maps shared_beta/shared_gamma flags
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
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    shared_beta = shared_effects,
    shared_gamma = shared_effects,
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
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_beta = shared_effects,
      shared_gamma = shared_effects,
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
#' @export
adaptive_tuner <- function(
  eval_fun,
  step_sizes = c(1, 0.1, 0.01),
  start_value = 0,
  end_value = 20, # if start < end then it's ascending.
  inc_streak_to_stop = 2,
  .warm_start = NULL,
  ... # all the other parameters being passed to eval_fun
) {
  results <- data.frame()
  best_overall <- list(parameter = NA, error = Inf, fit = NULL)
  ascending <- start_value < end_value
  direction <- if (ascending) 1 else -1
  current_start <- start_value
  fit <- .warm_start
  for (k in seq_along(step_sizes)) {
    if (ascending & current_start > end_value) {
      stop("For ascending search, start_value must be <= end_value.")
    }
    if (!ascending & current_start < end_value) {
      stop("For descending search, start_value must be >= end_value.")
    }
    step_size <- step_sizes[k]
    parameter <- current_start
    prev_error <- Inf
    inc_streak <- 0

    step_history <- data.frame(
      parameter = numeric(),
      error     = numeric(),
      step_size = numeric()
    )

    while ((ascending && parameter <= end_value) ||
      (!ascending && parameter >= end_value)) {
      out <- eval_fun(parameter, fit = fit, ...)
      fit <- out$fit
      # print(fit$d[1])
      error <- out$error

      step_history <- rbind(
        step_history,
        data.frame(
          parameter = parameter,
          error     = error,
          step_size = step_size
        )
      )

      if (error < best_overall$error) {
        best_overall$parameter <- parameter
        best_overall$error <- error
        best_overall$fit <- out$fit
      }

      if (error >= prev_error) {
        inc_streak <- inc_streak + 1
      } else {
        inc_streak <- 0
      }

      if (inc_streak >= inc_streak_to_stop) {
        break
      }

      prev_error <- error
      parameter <- parameter + direction * step_size
    }

    results <- rbind(results, step_history)
    ord <- order(step_history$error)
    best_idx <- ord[1]
    best_param <- step_history$parameter[best_idx]

    if (ascending) {
      current_start <- max(best_param - step_size, start_value)
      end_value <- min(best_param + step_size, end_value)
    } else {
      current_start <- min(best_param + step_size, start_value)
      end_value <- max(best_param - step_size, end_value)
    }
  }
  fit <- NULL # reset just in case
  list(
    best_parameter = best_overall$parameter,
    best_error     = best_overall$error,
    best_fit       = best_overall$fit,
    history        = results
  )
}
