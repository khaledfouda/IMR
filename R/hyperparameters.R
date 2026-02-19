#----------------------------------------------------------------
#' @export
generate_similarity <- function(x,
                             d = NULL,
                             matern_params = list(smoothness = 1.5, range = 1),
                             rbf_params = list(ell = 1),
                             jitter = 0,
                             invert = FALSE) {

    S <- NULL
    source_type <- "User Matrix"
    params_used <- list()

    if (is.matrix(x)) {
      if (nrow(x) != ncol(x)) stop("Input matrix 'x' must be square.")
      S <- x

    } else if (is.character(x)) {
      if (is.null(d)) stop("Distance matrix 'd' is required for kernel generation.")

      if (tolower(x) == "matern") {

        source_type <- "Matern Kernel"
        params_used <- matern_params

        S <- fields::Matern(d,
                            smoothness = matern_params$smoothness,
                            range = matern_params$range)

      } else if (tolower(x) == "rbf") {
        source_type <- "RBF Kernel"
        params_used <- rbf_params
        ell <- rbf_params$ell
        S <- exp(-(d^2) / (2 * ell^2))

      } else {
        stop("Unknown method. 'x' must be a matrix, 'matern', or 'RBF'.")
      }
      if (!invert) {
        warning(paste("Generated a raw Covariance matrix without inversion.",
                "fit_imr() expects the inverse."))
      }

    } else {
      stop("Input 'x' must be a matrix or a character string ('matern', 'RBF').")
    }

    if (invert) {
      S <- tryCatch({
        chol2inv(S)
      }, error = function(e) {
        stop("Matrix inversion failed (matrix might be singular).")
      })
      source_type <- paste(source_type, "(Inverted)")
    }
    if(is.numeric(jitter) && jitter > 0){
      S <- S + diag(jitter, nrow(S), ncol(S))
      #source_type <- paste(source_type, "(With Jitter)")
    }else jitter = 0

  decomp <- eigen(S, symmetric = TRUE)
  evals <- abs(decomp$values)
  cond_num <- max(evals) / min(evals[evals > 0])

    structure(
      list(
        U = decomp$vectors,
        d = decomp$values,
        meta = list(
          source = source_type,
          dim = dim(S),
          params = params_used,
          inverted = invert,
          jitter = jitter,
          condition_number = cond_num
        )
      ),
      class = "imr_similarity"
    )
}

#' @export
print.imr_similarity <- function(x, ...) {
  cat("\n== IMR Similarity Decomposition ==\n")
  cat(sprintf("Source:           %s\n", x$meta$source))
  cat(sprintf("Dimensions:       %d x %d\n", x$meta$dim[1], x$meta$dim[2]))
  cat(sprintf("Jitter value:     %s\n", x$meta$jitter))

  if (length(x$meta$params) > 0) {
    p_names <- names(x$meta$params)
    p_vals <- unlist(x$meta$params)
    cat(sprintf("Parameters:       %s\n", paste(p_names, p_vals, sep="=", collapse=", ")))
  }
  cond_fmt <- if (x$meta$condition_number > 1e4) "%.2e" else "%.2f"
  cat(sprintf("Condition Number: %s\n", sprintf(cond_fmt, x$meta$condition_number)))

  cat("Top 5 Eigenvalues:", paste(format(head(x$d, 5), digits = 3), collapse = ", "), "...\n")
  cat("==================================\n")

  invisible(x)
}



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



#' @export
get_imr_default_hparams <- function(similarity_row = NULL,
                                    similarity_col = NULL,
                                    lambda_row = 0,
                                    lambda_col = 0,
                                    tol = 1e-4) {
  laplacian_row <- if (is.null(similarity_row)) {
    list(U = NULL, d = NULL)
  } else {
    IMR::decompose_symmetric_matrix(similarity_row, lambda_row)
  }
  laplacian_col <- if (is.null(similarity_col)) {
    list(U = NULL, d = NULL)
  } else {
    IMR::decompose_symmetric_matrix(similarity_col, lambda_col)
  }

  list(
    beta = list(
      min   = 0,
      max = NULL, # if NULL, computed internally (recommended)
      #step_sizes = c(5,1),
      length = 10,
      #n_streaks = 2,
      value = 0 # if equal to max then no tuning to be done
    ),
    gamma = list(
      min   = 0,
      max = NULL, # if NULL, computed internally (recommended)
      #step_sizes = c(5,1),
      length = 10,
      #n_streaks = 2,
      value = 0 # if equal to max then no tuning to be done
    ),
    laplacian_row = laplacian_row,
    laplacian_col = laplacian_col,
    rank = list(
      min = 2,
      max = 30,
      default = 2,
      step_sizes = c(1),
      n_streaks = 2
    ),
    laplace = list(
      min = 0,
      max = 30,
      step_sizes = c(5, 2, 1, 0.1),
      n_streaks = 2
    )
  )
}
#-----------------------------------------------------
#' @export
decompose_symmetric_matrix <- function(x, lambda = 1) {
  stopifnot(isSymmetric(x))
  # if (lambda == 0) {
  #  return(list(U = NULL, d = NULL))
  # }
  xsvd <- base::eigen(x, symmetric = TRUE)
  return(list(U = xsvd$vectors, d = xsvd$values * lambda))
}
#-----------------------------------------------------
#' @export
get_lambda_M_max <-
  function(Y,
           X = NULL,
           Z = NULL,
           intercept_row = FALSE,
           intercept_col = FALSE,
           lambda_beta = NULL,
           lambda_gamma = NULL,

           maxit = 30) {
    need_fit <- any(!is.null(X), !is.null(Z), intercept_row, intercept_col)


    if (!need_fit) {
      return(svd::propack.svd(as.matrix(naive_MC(Y)),
        neig = 1, opts = list(kmax = maxit)
      )$d[1])
    }
    mfit <- imr.fit_no_low_rank(
      Y,
      X = X,
      Z = Z,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      trace = FALSE
    )
    # return largest singular value
    svd::propack.svd(as.matrix(naive_MC(mfit$resid)),
      neig = 1, opts = list(kmax = maxit)
    )$d[1]
  }
#------------------------------------------------
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
  init <- IMR::opt_svd(initial_fit$resid, r, nr, nc, FALSE, FALSE)
  initial_fit$u <- init$u
  initial_fit$d <- init$d
  initial_fit$v <- init$v
  initial_fit$resid =NULL
  # step 1: get an initial fit and find suitable lambda_M and rank before starting:
  if (is.null(y_valid)) {
    mfit <- list()
    # if no validation set provided then fit with a generic r and lambda_M
    mfit$fit <- IMR::imr.fit(
      Y = y_train,
      X = X,
      Z = Z,
      r = r,
      lambda_M = 0,
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
    lambda_M <- 0
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
    lambda_M <- mfit$lambda_M
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
      lambda_M      = lambda_M,
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
