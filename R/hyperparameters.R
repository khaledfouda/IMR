

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
  target <- match.arg(target)
  is_beta <- target == "beta"

  if (is_beta && !data$meta$has_X) stop("Target is 'beta' but no X matrix found in data.")
  if (!is_beta && !data$meta$has_Z) stop("Target is 'gamma' but no Z matrix found in data.")

  # --- 1. Baseline Fit (Holding Target at 0) ---
  # Fit the model with rank and intercepts, but keep our target penalized to 0
  baseline_fit <- imr_fit(
    data = data,
    rank = rank,
    lambda_m = lambda_m,
    lambda_beta = 0,
    lambda_gamma = 0,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    shared_beta = if(is_beta) shared_effects else FALSE,
    shared_gamma = if(!is_beta) shared_effects else FALSE,
    convergence = convergence
  )

  # --- 2. Calculate Residuals for KKT Conditions ---
  # Reconstruct predictions (M + intercepts)
  preds <- reconstruct(baseline_fit, data, trace = FALSE)$estimates

  # Calculate true dense residuals for the observed elements
  # (Assuming data$Y is a sparse matrix, we extract the non-zero locations)
  # For a full KKT check, we need the matrix multiplication against X or Z
  resid_mat <- matrix(0, nrow = data$meta$dimensions[1], ncol = data$meta$dimensions[2])
  resid_mat[data$obs_mask@i + 1, data$obs_mask@p] <- data$Y@x # Example mapping, adjust based on your exact sparse class
  resid_mat <- resid_mat - preds

  # --- 3. Theoretical KKT Max (The Anchor) ---
  # The theoretical minimum lambda that forces all coefs to 0 is max(abs(Covariates^T * Residuals))
  if (is_beta) {
    # CRITICAL FIX: Added abs()
    lambda_kkt <- max(abs(crossprod(data$Xq, resid_mat)))
  } else {
    lambda_kkt <- max(abs(resid_mat %*% data$Zq))
  }

  # --- 4. Bisection Search for Practical Supremum ---
  # Due to shared effects or scaling, practical max might deviate slightly from theory.
  # Bisection search is exponentially faster than grid search.

  lower <- 0
  upper <- lambda_kkt * 1.5 # 50% buffer above theory just in case

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
      shared_beta = if(is_beta) shared_effects else FALSE,
      shared_gamma = if(!is_beta) shared_effects else FALSE,
      warm_start = baseline_fit, # Use baseline to instantly converge
      convergence = convergence
    )
  }

  for (i in seq_len(bisection_iter)) {
    mid <- (lower + upper) / 2
    test_model <- fit_test(mid)

    # Check if all coefficients are effectively zero
    coefs <- if (is_beta) test_model$coefficients$beta else test_model$coefficients$gamma
    zero_ratio <- mean(abs(coefs) < 1e-6)

    if (zero_ratio == 1) {
      # It is fully sparse. We can try a LOWER lambda.
      upper <- mid
    } else {
      # It is NOT fully sparse. We must try a HIGHER lambda.
      lower <- mid
    }
  }

  # The upper bound is guaranteed to be the minimum lambda that achieves 100% sparsity
  lambda_sup <- upper

  if (verbose > 0) {
    message(sprintf(
      "Target: %s | KKT Theory: %.3f | Empiric Sup: %.3f (%.1f%% of KKT)",
      ifelse(is_beta, "Beta", "Gamma"), lambda_kkt, lambda_sup, 100 * lambda_sup / lambda_kkt
    ))
  }

  return(lambda_sup)
}
#' this return that max value for either row covariates or column covariates.
#' you must either privde X or Z. do not provide both
#' @export
get_lambda_lasso_max <- function(
  y_train,
  X = NULL,
  Z = NULL,
  y_valid = NULL,
  # W_valid = NULL,
  y = NULL,
  # row_cov = TRUE,
  intercept_row = TRUE,
  intercept_col = TRUE,
  shared_information = FALSE,
  hpar = get_imr_default_hparams(),
  r = 5,
  interior_loop_length = 20,
  maxit = 100,
  thresh = 1e-3,
  init_maxit = 100,
  init_thresh = 1e-4,
  verbose = 0,
  tol = 1
) {
  if (!xor(is.null(X), is.null(Z))) {
    stop("Either X or Z must be provided, but not both or neither.")
  }
  lambda_beta <- lambda_gamma <- 0
  row_cov <- is.null(Z)
  nr <- nrow(y_train)
  nc <- ncol(y_train)
  #-------
  # step 0: get initial values for all fits.
  initial_fit <- IMR::imr.fit_no_low_rank(y_train, X, Z, 0, 0,
                                          intercept_row = intercept_row,
                                          intercept_col = intercept_col,
                                          shared_information = shared_information,
                                          maxit = init_maxit,
                                          thresh = init_thresh)
  init <- IMR::svd_opt(initial_fit$resid, r, nr, nc, FALSE, FALSE)
  initial_fit$u <- init$u
  initial_fit$d <- init$d
  initial_fit$v <- init$v
  initial_fit$resid =NULL
  # step 1: get an initial fit and find suitable lambda_m and rank before starting:
  if (is.null(y_valid)) {
    mfit <- list()
    # if no validation set provided then fit with a generic r and lambda_m
    mfit$fit <- IMR::imr.fit(
      Y = y_train,
      X = X,
      Z = Z,
      r = r,
      lambda_m = 0,
      lambda_beta = 0,
      lambda_gamma = 0,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_information = shared_information,
      warm_start = initial_fit,
      maxit = maxit,
      thresh = thresh,
      trace = F
    )
    lambda_m <- 0
    r <- 2
  } else {
    mfit <- IMR::imr.cv_M(
      y_train = y_train,
      y_valid = y_valid,
      X = X,
      Z = Z,
      Y_full = y,
      lambda_beta = 0,
      lambda_gamma = 0,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      trace = F,
      maxit = maxit,
      thresh = thresh,
      hpar = hpar
    )
    lambda_m <- mfit$lambda_m
    r <- max(2, mfit$rank_M) # do not want the rank to be below 2
  }
  ##  Compute max value using kkt ------------
  residuals <- y_train -
    mfit$fit$u %*% (mfit$fit$d * t(mfit$fit$v))
  if (intercept_row) {
    residuals <- residuals - mfit$fit$beta0 %*% matrix(1, 1, nc)
  }
  if (intercept_col) {
    residuals <- residuals - matrix(1, nr, 1) %*% mfit$fit$gamma0
  }
  if (row_cov) {
    lambda_max <- max(crossprod(X, residuals))
  } else {
    lambda_max <- max(residuals %*% Z)
  }
  ## main_fit function -------------------
  onefit <- function(lambda_beta, lambda_gamma) {
    IMR::imr.fit(
      Y             = y_train,
      X             = X,
      Z             = Z,
      r             = r,
      lambda_m      = lambda_m,
      lambda_beta   = lambda_beta,
      lambda_gamma  = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_information = shared_information,
      maxit         = maxit,
      warm_start    = initial_fit,
      trace         = F
    )
  }

  ## 3) Line search for supremum value ---------------------------
  lambda_sup <- lambda_max
  mid_pt <- lambda_max / 2
  # we fit at mid point
  mfit <- onefit(
    if (row_cov) mid_pt else lambda_beta,
    if (row_cov) lambda_gamma else mid_pt
  )
  zero_ratio <- if (row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

  #-- now either search the right side or the left side
  count <- 0L
  if (zero_ratio < 1) {
    # search in [mid_pt, lambda_max]
    old_lambda <- lambda_max
    for (lam in seq(lambda_max, mid_pt, length.out = interior_loop_length)) {
      mfit <- onefit(
        if (row_cov) lam else lambda_beta,
        if (row_cov) lambda_gamma else lam
      )

      zero_ratio <- if (row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

      if (zero_ratio < 1) {
        lambda_sup <- old_lambda
        count <- count + 1L
      } else {
        count <- 0L
      }
      if (count >= tol) break
      old_lambda <- lam
    }
  } else {
    # search in [0, mid_pt]
    old_lambda <- mid_pt
    for (lam in seq(mid_pt, 0, length.out = interior_loop_length)) {
      mfit <- onefit(
        if (row_cov) lam else lambda_beta,
        if (row_cov) lambda_gamma else lam
      )

      zero_ratio <- if (row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

      if (zero_ratio < 1) {
        lambda_sup <- old_lambda
        count <- count + 1L
      } else {
        count <- 0L
      }
      if (count >= tol) break
      old_lambda <- lam
    }
  }

  if (verbose > 0) {
    message(sprintf(
      "%s: lambda: max = %.3f; sup = %.3f (%.1f%% of max)",
      ifelse(row_cov, "Beta", "Gamma"),
      lambda_max,
      lambda_sup,
      100 * lambda_sup / lambda_max
    ))
  }

  lambda_sup
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

#---
#' #' @export
#' adaptive_tuner_parallel <- function(
#'   eval_fun,
#'   step_sizes = c(1, 0.1, 0.01),
#'   start_value = 0,
#'   end_value = 20,
#'   parallel_seed = TRUE,
#'   parallel_packages = NULL,
#'   parallel_progress = TRUE,
#'   ...
#' ){
#'   results <- data.frame()
#'   best_overall <- list(parameter = NA, error = Inf, fit = NULL)
#'   current_start <- start_value
#'
#'   for(k in seq_along(step_sizes)){
#'     if(current_start > end_value)
#'       stop("start_value must be <= end_value (descending search not yet implemented).")
#'     step_size <- step_sizes[k]
#'     param_seq <- seq(current_start, end_value, step_size)
#'     if(length(param_seq) == 0)
#'       break
#'
#'     res_list <- parallel_grid_g(
#'       grid = list(lambda = param_seq),
#'       f = function(lambda, ...) eval_fun(lambda, ...),
#'       combine = "list",
#'       .seed = parallel_seed,
#'       .packages = parallel_packages,
#'       .progress = parallel_progress,
#'       ...
#'     )
#'
#'     errors <- vapply(
#'       res_list,
#'       function(x) {
#'         err <- x[[2L]]
#'         if (!is.finite(err)) Inf else err
#'       },
#'       numeric(1L)
#'     )
#'     fits <- lapply(res_list, `[[`, 1L)
#'
#'     step_history <- data.frame(
#'       parameter        = param_seq,
#'       error            = errors,
#'       step_size        = step_size,
#'       resolution_level = k
#'     )
#'
#'     # accumulate history
#'     results <- rbind(results, step_history)
#'
#'     # best for this resolution
#'     best_idx_k <- which.min(errors)
#'     best_param <- param_seq[best_idx_k]
#'
#'     # update global best
#'     if (errors[best_idx_k] < best_overall$error) {
#'       best_overall$parameter <- best_param
#'       best_overall$error     <- errors[best_idx_k]
#'       best_overall$fit       <- fits[[best_idx_k]]
#'     }
#'
#'     # narrow the window around the best param (same idea as your sequential fn)
#'     current_start <- max(best_param - step_size, start_value)
#'     end_value     <- best_param + step_size
#'   }
#'
#'   list(
#'     best_parameter = best_overall$parameter,
#'     best_error     = best_overall$error,
#'     best_fit       = best_overall$fit,
#'     history        = results
#'   )
#' }
