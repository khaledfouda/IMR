#----------------------------------------------------------------
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
    M = list(
      lambda_max = NULL, # = lambda_{M,min}
      lambda_factor = 1 / 4, # ignored if lambda_init is provided
      n.lambda = 20, # sequence from lambda_init to 0 (inclusive)
      rank.init = 2,
      rank.max = 30,
      rank.min = 2,
      rank.step = 2,
      early.stopping = 1
    ),
    beta = list(
      lambda_max = NULL, # if NULL, computed internally (recommended)
      n.lambda   = 20,
      init.tol   = 3
    ),
    gamma = list(
      lambda_max = NULL,
      n.lambda   = 20,
      init.tol   = 3
    ),
    laplacian_row = laplacian_row,
    laplacian_col = laplacian_col,
    rank = list(
      step_sizes = c(2, 1),
      rank.min = 2,
      rank.max = 30,
      n_streaks = 2
    ),
    laplace = list(
      step_sizes = c(5, 2, 1, 0.1, 0.01),
      start_value = 0,
      end_value = 30,
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
    hpar = get_imr_default_hparams(),
    interior_loop_length = 20,
    maxit = 100,
    verbose = 0,
    tol = 1) {
  if (!xor(is.null(X), is.null(Z))) {
    stop("Either X or Z must be provided, but not both or neither.")
  }
  lambda_beta <- lambda_gamma <- 0
  row_cov <- is.null(Z)
  nr <- nrow(y_train)
  nc <- ncol(y_train)
  # step 1: get an initial fit and find suitable lambda_M and rank before starting:
  if (is.null(y_valid)) {
    mfit <- list()
    # if no validation set provided then fit with a generic r and lambda_M
    mfit$fit <- IMR::imr.fit(
      Y = y_train,
      X = X,
      Z = Z,
      r = 5,
      lambda_M = 0,
      lambda_beta = 0,
      lambda_gamma = 0,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = maxit,
      trace = F
    )
    lambda_M <- 0
    r <- 5
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
      maxit         = maxit,
      trace         = FALSE
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
    ... # all the other parameters being passed to eval_fun
    ) {
  results <- data.frame()
  best_overall <- list(parameter = NA, error = Inf, fit = NULL)
  ascending = start_value < end_value
  direction = if(ascending) 1 else -1
  current_start <- start_value

  for (k in seq_along(step_sizes)) {
    if (ascending & current_start > end_value) {
      stop("For ascending search, start_value must be <= end_value.")
    }
    if ( !ascending & current_start < end_value) {
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

    while (( ascending && parameter <= end_value) ||
           (!ascending && parameter >= end_value)) {
      out <- eval_fun(parameter, ...)
      fit <- out[[1]]
      error <- out[[2]]

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
        best_overall$fit <- fit
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

    if(ascending){
      current_start <- max(best_param - step_size, start_value)
      end_value <- min(best_param + step_size, end_value)
    }else{
      current_start <- min(best_param + step_size, start_value)
      end_value <- max(best_param - step_size, end_value)
    }

    # second_idx <- if (length(ord) > 1) ord[2] else best_idx
    # current_start <- step_history$parameter[second_idx]
    # end_value <- step_history$parameter[best_idx] + step_size
  }

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
