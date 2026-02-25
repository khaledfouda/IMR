

#' @export
imr_hparameters <- function(beta_range = c(0, NA),
                             beta_default = 0,
                             gamma_range = c(0, NA),
                             gamma_default = 0,
                             rank_range = c(2, 30, 1, 2),
                             grid_length = 10,
                             laplace_range = c(0,10),
                             laplace_steps = c(1, 0.1)
) {


  structure(
    list(
      beta = list(
        min = beta_range[1],
        max = if (length(beta_range)==1){beta_range[1] }else
          if(is.na(beta_range[2])) "auto" else beta_range[2],
        length = if(length(beta_range)==1) 1 else grid_length,
        default = if(length(beta_range)==1) beta_range[1] else beta_default
      ),
      gamma = list(
        min = gamma_range[1],
        max = if (length(gamma_range)==1){gamma_range[1] }else
          if(is.na(gamma_range[2])) "auto" else gamma_range[2],
        length = if(length(gamma_range)==1) 1 else grid_length,
        default = if(length(gamma_range)==1) gamma_range[1] else gamma_default
      ),
      rank = list(
        min = rank_range[1],
        max = rank_range[2],
        step_sizes = if(length(rank_range)>2) rank_range[3] else 1,
        n_streaks = if(length(rank_range)>3) rank_range[4] else 2
      ),
      laplace = list(
        min = laplace_range[1],
        max = laplace_range[2],
        steps = laplace_steps
      )
    ),
    class = "imr_hparameters"
  )
}

#' @export
print.imr_hparameters <- function(x, ...) {
  cat("\n== IMR Hyperparameter Configuration ==\n")

  fmt_range <- function(param) {
    if (param$length == 1) {
      return(paste0("Fixed at ", param$min))
    }

    max_val <- if (identical(param$max, "auto")) "Auto" else param$max

    sprintf("Range: %s -> %s (Grid: %d points) | Default: %s",
            param$min, max_val, param$length, param$default)
  }

  cat(sprintf("%-18s %s\n", "Beta (Row Reg):", fmt_range(x$beta)))
  cat(sprintf("%-18s %s\n", "Gamma (Col Reg):", fmt_range(x$gamma)))
  cat(sprintf("%-18s Range: %d -> %d (Step: %s, Streaks: %s)\n",
              "Rank:",
              x$rank$min,
              x$rank$max,
              x$rank$step_sizes,
              x$rank$n_streaks))
  steps_str <- paste(x$laplace$steps, collapse = ", ")
  cat(sprintf("%-18s Range: %s -> %s (Steps: %s)\n",
              "Laplace:",
              x$laplace$min,
              x$laplace$max,
              steps_str))

  cat("========================================\n\n")

  invisible(x)
}

#-----------------------------------------------------
#' @export
get_lambda_m_max <-
  function(data,
           intercept_row = FALSE,
           intercept_col = FALSE,
           lambda_beta = 0,
           lambda_gamma = 0,
           shared_beta = FALSE,
           shared_gamma = FALSE,
           convergence = IMR::imr_convergence(trace = FALSE, ls_initial = FALSE)) {

    need_fit <- any(!is.null(data$Xq), !is.null(data$Zq), intercept_row, intercept_col)


    if (!need_fit) {
      return(IMR::svd_opt(data$Y,1)$d[1])
    }

    fit <- IMR::imr_fit(data, rank = 0,
                        lambda_beta = lambda_beta,
                        lambda_gamma = lambda_gamma,
                        intercept_row = intercept_row,
                        intercept_col = intercept_col,
                        shared_beta = shared_beta,
                        shared_gamma = shared_gamma,
                        convergence = convergence)

    # return largest singular value
    return(IMR::svd_opt(fit$residuals,1)$d[1])
  }
#------------------------------------------------
#' Find the minimum Lasso lambda that forces all covariates to zero
#' @export
get_lambda_lasso_max <- function(
    data,                     # Must be an 'imr_data' S3 object
    target = c("beta", "gamma"),
    rank = 2,
    lambda_m = 0,
    intercept_row = TRUE,
    intercept_col = TRUE,
    shared_effects = FALSE,   # Maps shared_beta/shared_gamma flags
    bisection_iter = 15,
    convergence = imr_convergence(),
    verbose = 0
) {

  stopifnot(inherits(data, "imr_data"))
  target <- match.arg(target, c("beta", "gamma"))
  is_beta <- target == "beta"
  # remove the other set of covariates from the data.
  if(is_beta) data$Zq = NULL else data$Xq = NULL

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
  resids <- IMR::reconstruct_partial(baseline_fit, data, data$obs_mask,trace = FALSE)
  resids@x <- data$Y@x - resids@x

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
  for(i in seq_len(bisection_iter)){
    test_model <- fit_test(upper)
    # Check if all coefficients are effectively zero
    coefs <- if (is_beta) test_model$coefficients$beta else test_model$coefficients$gamma
    zero_ratio <- mean(abs(coefs) < 1e-6)

    if (zero_ratio == 1) {
      # we achieved an upper-bound. Let's return.
      break
    } else {
      # It is not fully sparse. We must try a HIGHER lambda.
      upper <- lambda_kkt * (i+1)
    }
    if(verbose >= 2)
      message(sprintf("lambda = %.4f, zero ratio = %.2f",upper, zero_ratio))
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
    if(verbose >= 2)
      message(sprintf("lambda = %.4f, zero ratio = %.2f",mid, zero_ratio))
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
  fit = .warm_start
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
      out <- eval_fun(parameter, fit=fit,  ...)
      fit <- out$fit
      #print(fit$d[1])
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

