#----------------------------------------------------------
#' @noRd
imr_tune_nuclear_fast <- function(data,
                                  grid,
                                  lambda_beta = 0,
                                  lambda_gamma = 0,
                                  final_fit = TRUE,
                                  convergence = imr_convergence(),
                                  error_function = get_metric("rmse"),
                                  warm_start = NULL,
                                  verbose = 1,
                                  log_grid = TRUE,
                                  seed = NULL) {
  #-----------------------------------------------------
  # input verification
  if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
  stopifnot(
    inherits(data, "imr_data"),
    inherits(grid, "imr_tune_grid"),
    inherits(convergence, "imr_convergence")
  )
  stopifnot(
    is_incomplete(data$y_valid),
    is_incomplete(data$y_train)
  )
  stopifnot(is.numeric(grid$nuclear$max))
  stopifnot(
    grid$rank$min > 0,
    grid$rank$max >= grid$rank$min,
    grid$rank$step >= 0,
    grid$nuclear$length >= 1,
    grid$nuclear$max >= grid$nuclear$min
  )
  #---------------------------------------------------
  # indices
  if (grid$nuclear$max <= 0) grid$nuclear$max <- 1e-4
  if (grid$nuclear$min <= 0) grid$nuclear$min <- 1e-6
  if (log_grid) {
    lambda_seq <- exp(seq(log(grid$nuclear$max),
      log(grid$nuclear$min),
      length.out = grid$nuclear$length
    ))
  } else {
    lambda_seq <- seq(grid$nuclear$max, grid$nuclear$min,
      length.out = grid$nuclear$length
    )
  }
  virow <- data$y_valid@i
  vpcol <- data$y_valid@p
  reference <- data$y_valid@x
  # initial values
  mfit <- warm_start
  rank_max <- grid$rank$min
  history <- data.frame(
    lambda_m = lambda_seq,
    verror = rep(NA_real_, length(lambda_seq)),
    rank_in = rep(NA_integer_, length(lambda_seq)),
    rank_out = rep(NA_integer_, length(lambda_seq))
  )
  best_fit_obj <- NULL
  best_params <- NULL
  best_verror <- Inf
  no_improve_count <- 0
  # main loop
  for (i in seq_along(lambda_seq)) {
    mfit <- imr_fit(data,
      rank = rank_max,
      lambda_m = lambda_seq[i],
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      convergence = convergence,
      training = TRUE,
      warm_start = mfit
    )
    # compute validation error
    vestim <- reconstruct_partial(mfit, data, virow, vpcol)
    verror <- error_function(vestim, reference)
    # compute rank
    rank_out <- sum(mfit$coefficients$d > 1e-4)
    # verbose
    if (verbose >= 2) {
      message(sprintf(
        "%2d lambda_m=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
        i,
        lambda_seq[i],
        rank_max,
        rank_out,
        verror,
        mfit$meta$n_iter
      ))
    }
    # update history
    history$verror[i] <- verror
    history$rank_in[i] <- rank_max
    history$rank_out[i] <- rank_out
    # track best model & early stopping
    completely_sparse <- rank_out == 0
    if (verror < best_verror) {
      best_fit_obj <- mfit
      best_verror <- verror
      best_params <- history[i, ]
      no_improve_count <- 0
    } else if (!completely_sparse) {
      no_improve_count <- no_improve_count + 1
    }

    if (no_improve_count >= grid$nuclear$streaks) {
      if (verbose >= 2) {
        message(
          sprintf(
            "Early Stopping: no improvement in last %d iterations.",
            no_improve_count
          )
        )
      }
      break
    }
    # update rank_max for next iteration
    rank_max <- min(rank_out + grid$rank$step, grid$rank$max)
    rank_max <- max(rank_max, grid$rank$min)
  }
  if (verbose > 0) {
    message(sprintf(
      "Best fit: lambda_m=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
      best_params$lambda_m,
      best_params$rank_in,
      best_params$rank_out,
      best_verror,
      best_fit_obj$meta$n_iter
    ))
  }
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    best_fit_obj <- imr_fit(data,
      rank = best_params$rank_in,
      lambda_m = best_params$lambda_m,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      convergence = convergence,
      training = FALSE,
      warm_start = best_fit_obj
    )
  }
  # Clean history of NAs from early stopping
  history <- history[!is.na(history$verror), ]

  list(fit = best_fit_obj, params = best_params, history = history)
}
# ==============================================================================================
#----------------------------------------------------------
#' @noRd
imr_tune_nuclear_slow <- function(data,
                                  grid,
                                  lambda_beta = 0,
                                  lambda_gamma = 0,
                                  final_fit = TRUE,
                                  convergence = imr_convergence(),
                                  error_function = get_metric("rmse"),
                                  warm_start = NULL,
                                  log_grid = TRUE,
                                  verbose = 1,
                                  seed = NULL) {
  #-----------------------------------------------------
  # input verification
  if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
  stopifnot(
    inherits(data, "imr_data"),
    inherits(grid, "imr_tune_grid"),
    inherits(convergence, "imr_convergence")
  )
  stopifnot(
    is_incomplete(data$y_valid),
    is_incomplete(data$y_train)
  )
  stopifnot(is.numeric(grid$nuclear$max))
  stopifnot(
    grid$rank$min > 0,
    grid$rank$max >= grid$rank$min,
    grid$rank$step >= 0,
    grid$nuclear$length >= 1,
    grid$nuclear$max >= grid$nuclear$min
  )
  #---------------------------------------------------
  # training grids
  if (grid$nuclear$max <= 0) grid$nuclear$max <- 1e-4
  if (grid$nuclear$min <= 0) grid$nuclear$min <- 1e-6
  if (log_grid) {
    lambda_seq <- exp(seq(log(grid$nuclear$max),
                          log(grid$nuclear$min),
                          length.out = grid$nuclear$length
    ))
  } else {
    lambda_seq <- seq(grid$nuclear$max, grid$nuclear$min,
                      length.out = grid$nuclear$length
    )
  }
  rank_seq <- seq(grid$rank$min, grid$rank$max, grid$rank$step)
  #-----------------------------------------------------------
  # indices

  virow <- data$y_valid@i
  vpcol <- data$y_valid@p
  reference <- data$y_valid@x
  # initial values
  mfit <- warm_start
  history <- expand.grid(
    rank = rank_seq,
    lambda_m = lambda_seq
  )
  history <- history[, c("lambda_m", "rank")]
  history$verror <- rep(NA_real_, nrow(history))
  history$rank_in <- rep(NA_integer_, nrow(history))
  history$rank_out <- rep(NA_integer_, nrow(history))
  n_ranks <- length(rank_seq)
  best_fit_obj_1 <- NULL
  best_params_1 <- NULL
  best_verror_1 <- Inf
  no_improve_count_1 <- 0
  # main loop
  for (i in seq_along(lambda_seq)) {
    current_lambda <- lambda_seq[i]
    best_fit_obj_2 <- NULL
    best_params_2 <- NULL
    best_verror_2 <- Inf
    no_improve_count_2 <- 0

    for (j in seq_along(rank_seq)) {
      current_rank <- rank_seq[j]
      row_idx <- (i - 1) * n_ranks + j

      mfit <- imr_fit(data,
        rank = current_rank,
        lambda_m = current_lambda,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        convergence = convergence,
        training = TRUE,
        warm_start = mfit
      )
      # compute validation error
      vestim <- reconstruct_partial(mfit, data, virow, vpcol)
      verror <- error_function(vestim, reference)
      # compute rank
      rank_out <- sum(mfit$coefficients$d > 1e-4)
      # verbose
      if (verbose >= 2) {
        message(sprintf(
          "%2d lambda_m=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
          i,
          lambda_seq[i],
          current_rank,
          rank_out,
          verror,
          mfit$meta$n_iter
        ))
      }
      # update history
      history$verror[row_idx] <- verror
      history$rank_in[row_idx] <- current_rank
      history$rank_out[row_idx] <- rank_out
      # track best model & early stopping
      completely_sparse <- rank_out == 0
      if (verror < best_verror_2) {
        best_fit_obj_2 <- mfit
        best_verror_2 <- verror
        best_params_2 <- history[row_idx, ]
        no_improve_count_2 <- 0
      } else if (!completely_sparse) {
        no_improve_count_2 <- no_improve_count_2 + 1
      }

      if (no_improve_count_2 >= grid$rank$streaks) {
        if (verbose >= 2) {
          message(
            sprintf(
              "Early Stopping in inner loop: no improvement in last %d iterations.",
              no_improve_count_2
            )
          )
        }
        break
      }
    }
    #----
    if (best_verror_2 < best_verror_1) {
      best_fit_obj_1 <- best_fit_obj_2
      best_verror_1 <- best_verror_2
      best_params_1 <- best_params_2
      no_improve_count_1 <- 0
    } else if (best_params_2$rank_out != 0) {
      no_improve_count_1 <- no_improve_count_1 + 1
    }

    if (no_improve_count_1 >= grid$nuclear$streaks) {
      if (verbose >= 2) {
        message(
          sprintf(
            "Early Stopping in outer loop: no improvement in last %d iterations.",
            no_improve_count_1
          )
        )
      }
      break
    }
  }

  if (verbose > 0) {
    message(sprintf(
      "Best fit: lambda_m=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
      best_params_1$lambda_m,
      best_params_1$rank_in,
      best_params_1$rank_out,
      best_verror_1,
      best_fit_obj_1$meta$n_iter
    ))
  }
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    best_fit_obj_1 <- imr_fit(data,
      rank = best_params_1$rank_in,
      lambda_m = best_params_1$lambda_m,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      convergence = convergence,
      training = FALSE,
      warm_start = best_fit_obj_1
    )
  }
  # Clean history of NAs from early stopping
  history <- history[!is.na(history$verror), ]

  list(fit = best_fit_obj_1, params = best_params_1, history = history)
}
# ==============================================================================================
#' @noRd
imr_tune_lasso <- function(data,
                           grid,
                           target = c("beta", "gamma"),
                           fixed_other_lasso = 0,
                           final_fit = TRUE,
                           use_warm_in_final = TRUE,
                           convergence = imr_convergence(),
                           error_function = get_metric("rmse"),
                           warm_start = NULL,
                           verbose = 1,
                           n_cores = 4,
                           fast_nuclear = TRUE,
                           nuclear_log_scale = TRUE,
                           seed = NULL) {
  #-----------------------------------------------------
  # input verification and seed
  if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
  stopifnot(
    inherits(data, "imr_data"),
    inherits(grid, "imr_tune_grid"),
    inherits(convergence, "imr_convergence")
  )
  stopifnot(
    is_incomplete(data$y_valid),
    is_incomplete(data$y_train)
  )
  target <- stringr::str_to_lower(target)
  stopifnot(target %in% c("beta", "gamma"))
  is_beta <- target == "beta"
  lambda_obj <- if (is_beta) grid$beta else grid$gamma
  stopifnot(
    is.numeric(lambda_obj$max),
    lambda_obj$min >= 0,
    lambda_obj$max >= lambda_obj$min,
    lambda_obj$length >= 1
  )
  #---------------------------------------------------
  # Fallback for Windows [mclapply not supported on windows]
  if (.Platform$OS.type == "windows" && n_cores > 1) {
    warning("mclapply forking is not supported on Windows. Falling back to sequential execution.")
    n_cores <- 1L
  }
  if (verbose > 0) {
    message(sprintf(
      "Parallel Nested Search: %d lambda %s values using %d cores...",
      lambda_obj$length, target, n_cores
    ))
  }
  #--------------------------------------------------
  # Run the loop
  nuclear_function <- if (fast_nuclear) imr_tune_nuclear_fast else imr_tune_nuclear_slow
  # if(!final_fit) use_warm_in_final = FALSE # not needed.
  lambda_seq <- seq(lambda_obj$min, lambda_obj$max, length.out = lambda_obj$length)
  results_list <- parallel::mclapply(seq_along(lambda_seq), function(i) {
    lambda <- lambda_seq[i]
    if (verbose >= 2) {
      cat(sprintf("Worker started: %s = %.4f\n", target, lambda))
    }
    # run full tune_nuclear
    nuclear_res <- nuclear_function(
      data = data,
      grid = grid,
      lambda_beta = if (is_beta) lambda else fixed_other_lasso,
      lambda_gamma = if (is_beta) fixed_other_lasso else lambda,
      final_fit = FALSE,
      convergence = convergence,
      error_function = error_function,
      warm_start = warm_start,
      log_grid = nuclear_log_scale,
      verbose = 0,
      seed = seed
    )
    # return best results
    best_inner <- nuclear_res$params
    best_inner[[paste0("lambda_", target)]] <- lambda
    best_inner[[paste0("lambda_", if (is_beta) "gamma" else "beta")]] <- fixed_other_lasso
    # we also return the full history for diagnostics/plots
    history <- nuclear_res$history
    history[[paste0("lambda_", target)]] <- lambda
    history[[paste0("lambda_", if (is_beta) "gamma" else "beta")]] <- fixed_other_lasso

    out <- list(
      best_inner = best_inner,
      history = history
    )
    if (use_warm_in_final) {
      out$fit <- nuclear_res$fit
    }

    return(out)
  },
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  mc.set.seed = seed,
  mc.silent = if (n_cores > 1) TRUE else FALSE,
  mc.cleanup = TRUE
  )
  #-----------------------------------------------
  # Aggregate results
  best_inner <- do.call(rbind, lapply(results_list, function(x) x$best_inner))
  history <- do.call(rbind, lapply(results_list, function(x) x$history))
  #------------------------
  # find best values
  best_idx <- which.min(best_inner$verror)
  best_params <- best_inner[best_idx, ]
  if (verbose > 0) {
    message(sprintf(
      "Best parameters: %s: %.5f | nuclear: %.4f | Error: %.5f",
      target, best_params[[paste0("lambda_", target)]],
      best_params$lambda_m, best_params$verror
    ))
  }
  #---------------------------------------------
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    if (use_warm_in_final) {
      warm_start <- results_list[[best_idx]]$fit
    }
    best_fit_obj <- imr_fit(data,
      rank = best_params$rank_in,
      lambda_m = best_params$lambda_m,
      lambda_beta = best_params$lambda_beta,
      lambda_gamma = best_params$lambda_gamma,
      convergence = convergence,
      warm_start = warm_start
    )
  } else {
    best_fit_obj <- NULL
  }

  list(fit = best_fit_obj, params = best_params, history = history)
}
# =============================================================================================
#' @title Hyperparameter Tuning for IMR Models
#'
#' @description
#' Executes hyperparameter optimization for Incomplete Matrix Regression (IMR)
#' models. The procedure evaluates predictive performance on a validation set
#' (`y_valid`) while estimating the model on a training set (`y_train`).
#'
#' @param data An object of class `"imr_data"` containing the training and
#'   validation partitions.
#' @param grid An object of class `"imr_tune_grid"`. **Note:** Any `"auto"`
#'   specifications for maximum values within this grid must be resolved via
#'   \code{\link{imr_set_grid_limits}} prior to invoking `imr_tune`.
#' @param final_fit Logical. If `TRUE`, the function performs a final model
#'   estimation on the complete dataset (`Y`) using the identified optimal
#'   hyperparameters. Defaults to `TRUE`.
#' @param use_warm_in_final Internal. Pending deprecation.
#' @param fast_nuclear Logical. Specifies the algorithmic approach for tuning the
#'   low-rank component (nuclear norm penalty). Defaults to `TRUE`.
#' @param convergence An `"imr_convergence"` object specifying numerical tolerances
#'   for internal model estimation.
#' @param error_function A function used to evaluate prediction error on the
#'   validation set. Must accept two arguments, `(predicted, actual)`.
#'   Defaults to `get_metric("rmse")`.
#' @param warm_start Internal. Pending deprecation.
#' @param verbose Integer. Controls the level of diagnostic progress output.
#'   Defaults to `1`.
#' @param n_cores Integer. Number of CPU cores allocated for parallel execution
#'   across covariate grids. Defaults to `4`.
#' @param seed Integer. Seed for random number generation to ensure reproducibility.
#'   Defaults to `NULL`.
#' @param nuclear_log_scale Logical. If `TRUE`, the \eqn{\lambda_M} (`lambda_m`)
#'   grid is constructed on a logarithmic scale, increasing the density of evaluation
#'   points near the lower bound. Recommended when the maximum \eqn{\lambda_M}
#'   is large. Defaults to `TRUE`.
#' @param tune_maxit Integer. Maximum number of alternating optimization iterations
#'   permitted when simultaneously tuning all three parameters (Scenario 3).
#'   Defaults to `10`.
#' @param tune_tol Numeric. Relative tolerance threshold for early stopping during
#'   alternating optimization. Optimization terminates if the absolute relative
#'   change in validation error falls below this limit. Defaults to `1e-4`.
#'
#' @details
#' **Tuning Scenarios:**
#' The function selects an optimization trajectory based on the complexity of
#' the `grid` and model specification:
#' \enumerate{
#'   \item If \eqn{\lambda_\beta} and \eqn{\lambda_\Gamma} are fixed or
#'     absent from the model, optimization proceeds sequentially over the
#'     low-rank component parameters.
#'   \item  If exactly one covariate penalty
#'     (\eqn{\lambda_\beta} or \eqn{\lambda_\Gamma}) requires tuning, the function parallelizes
#'     execution across the corresponding covariate grid. For each value,
#'     Scenario 1 is executed to optimize the nuclear parameters.
#'   \item When \eqn{\lambda_\beta},
#'     \eqn{\lambda_\Gamma}, and the nuclear component all require tuning, an
#'     alternating loop is utilized. The algorithm fixes \eqn{\lambda_\Gamma} and
#'     optimizes \eqn{\lambda_\beta} + nuclear (Scenario 2), then fixes \eqn{\lambda_\beta}
#'     and optimizes \eqn{\lambda_\Gamma} + nuclear. This iterates until `tune_tol`
#'     or `tune_maxit` is reached.
#' }
#'
#' **Nuclear Tuning Modes (`fast_nuclear`):**
#' \itemize{
#'   \item \strong{Fast Mode (`TRUE`):} Commences at the maximum \eqn{\lambda_M}
#'     and minimum rank. The algorithm simultaneously decrements \eqn{\lambda_M}
#'     while increasing the rank, terminating when validation performance
#'     stagnates according to the grid's patience parameter.
#'   \item \strong{Slow Mode (`FALSE`):} Constructs a nested grid.
#'     For each \eqn{\lambda_M} value, the function iteratively evaluates ranks
#'     until predictive performance degrades. This mode provides a comprehensive
#'     mapping of the parameter space.
#' }
#'
#' @return A list comprising:
#' \itemize{
#'   \item \code{all_params}: A data frame recording the optimal parameter
#'     configurations identified at each major iteration of the tuning process.
#'   \item \code{history}: A comprehensive data frame of all evaluated parameter
#'     combinations and their respective validation errors.
#'   \item \code{fit}: The final estimated `"imr_fit"` object evaluated at the
#'     global optimal parameters (if `final_fit = TRUE`).
#'   \item \code{params}: A data frame containing the global optimal
#'     hyperparameter combination.
#'   \item \code{time_secs}: Total execution time in seconds.
#' }
#' @seealso [imr_tune_grid()], [imr_set_grid_limits()], [imr_data()]
#' @examples
#' # create sample data
#' Y <- matrix(
#'   c(2, NA, 3, 4,
#'     4, .5, NA, 4,
#'     NA, NA, 5, 3), 3, byrow= TRUE
#' )
#' # create a data object
#' data <- imr_data(Y =  Y, val_prop = 0.2)
#' # create a grid of hyperparameters
#' grid <- imr_tune_grid(nuclear = c(0, NA, 5, 2),
#'                      rank = c(2, 5, 1, 2))
#'
#' # get the KKT max value for the nuclear parameter
#' grid <- imr_set_grid_limits(data, grid, bisection_iter=0)
#'
#' # tune the parameters lambda_m and r on the model Y = M
#' cv_out <- imr_tune(data, grid, fast_nuclear = TRUE, n_cores = 1)
#'
#' @export
imr_tune <- function(data,
                     grid,
                     final_fit = TRUE,
                     use_warm_in_final = TRUE,
                     fast_nuclear = TRUE,
                     convergence = imr_convergence(),
                     error_function = get_metric("rmse"),
                     warm_start = NULL,
                     verbose = 1,
                     n_cores = 4,
                     seed = NULL,
                     nuclear_log_scale = TRUE,
                     tune_maxit = 10,
                     tune_tol = 1e-4) {
  if (!is.null(seed) && is.numeric(seed)) set.seed(seed)
  stopifnot(inherits(data, "imr_data"), inherits(grid, "imr_tune_grid"))
  t_start_global <- Sys.time()
  # --- Determine which parameters to tune
  tune_beta <- data$model$row_covariates && data$meta$has_X && grid$beta$length > 1
  tune_gamma <- data$model$col_covariates && data$meta$has_Z && grid$gamma$length > 1

  # if(data$model$row_covariates && grid$beta$length == 1)
  #   default_lambda_beta <- grid$beta$min
  # if(data$model$col_covariates && grid$gamma$length == 1)
  #   default_lambda_gamma <- grid$gamma$min

  #----------------------------------------------------------
  # Scenario 1: Tune nuclear Only
  #------------------------------------------------------------
  if (!tune_beta && !tune_gamma) {
    if (verbose > 0) message("Tuning  nuclear (M) only...")
    nuclear_function <- if (fast_nuclear) imr_tune_nuclear_fast else imr_tune_nuclear_slow
    out_obj <- nuclear_function(
      data = data, grid = grid,
      lambda_beta = grid$beta$min, lambda_gamma = grid$gamma$min,
      final_fit = final_fit, convergence = convergence,
      error_function = error_function, warm_start = warm_start,
      log_grid = nuclear_log_scale,
      verbose = verbose, seed = seed
    )
    t_total <- round(difftime(Sys.time(), t_start_global), 2)
    if (verbose > 0) {
      message(sprintf("\nTotal Tuning Time: %s", format(t_total)))
    }
    out_obj$time_secs <- as.numeric(t_total,"secs")
    return(out_obj)
  }
  #---------------------------------------------------------------------
  # Scenario 2: Tuning nuclear + one of lambda_beta or lambda_gamma
  #-----------------------------------------------------------------
  if (tune_beta != tune_gamma) {
    target <- if (tune_beta) "beta" else "gamma"
    fixed_other <- if (tune_beta) grid$gamma$min else grid$beta$min

    if (verbose > 0) message(sprintf("Tuning %s + nuclear...", (target)))
    out_obj <- imr_tune_lasso(
      data = data, grid = grid, target = target, fixed_other_lasso = fixed_other,
      final_fit = final_fit, use_warm_in_final = use_warm_in_final,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose, fast_nuclear = fast_nuclear,
      n_cores = n_cores, seed = seed
    )
    t_total <- difftime(Sys.time(), t_start_global)
    if (verbose > 0) {
      message(sprintf("\nTotal Tuning Time: %s", format(t_total)))
    }
    out_obj$time_secs <- as.numeric(t_total,"secs")
    return(out_obj)
  }
  #-----------------------------------------------------------------------
  # Scenario 3: Tuning all 3 parameters. An iterative method will be used.
  #------------------------------------------------------------------------
  if (verbose > 0) message("Scenario 3: Alternating optimization for all 3 parameters...")

  # Initialize
  cur_gamma <- grid$gamma$min
  # cur_beta <- diff_beta <- diff_gamma <- Inf
  old_verror <- 9999
  all_history <- all_params <- data.frame()
  # one_more_fit <- FALSE # this is an extra for for the purpose of final_fit=TRUE
  last_output <- NULL
  keep_iterating <- TRUE
  iter <- 1
  while (keep_iterating) {
    if (verbose > 0) message(sprintf("\n--- Tune Iteration %d ---", iter))

    # old_beta <- cur_beta
    # old_gamma <- cur_gamma

    # --- Step A: Tune Beta (given current Gamma) ---
    t_start_iter <- Sys.time()
    res_beta <- imr_tune_lasso(
      data = data, grid = grid, target = "beta", fixed_other_lasso = cur_gamma,
      final_fit = final_fit, # if (one_more_fit) final_fit else FALSE,
      use_warm_in_final = use_warm_in_final, fast_nuclear = fast_nuclear,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose - 1,
      n_cores = n_cores, seed = seed
    )
    cur_beta <- res_beta$params$lambda_beta
    iter_time_secs <- as.numeric(difftime(Sys.time(), t_start_iter, units = "secs"))
    #---------------------------------------------------------------------------
    # track parameter squared difference for convergence.
    # diff_beta <- (cur_beta - old_beta)^2
    # diff <- diff_beta + diff_gamma
    diff <- abs(res_beta$params$verror - old_verror) / old_verror
    old_verror <- res_beta$params$verror
    # track history for debugging/plotting
    # res_beta$history$step <- 1
    res_beta$history$iter <- iter
    all_history <- rbind(all_history, res_beta$history)
    # track best performance rows as well
    # res_beta$params$step <- 1
    res_beta$params$iter <- iter
    res_beta$params$diff <- diff
    all_params <- rbind(all_params, res_beta$params)
    # convergence check
    if (verbose > 0 && iter > 1) {
      message(sprintf(
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | nuclear: %.4f | Diff: %.6f | Time: %.2fs",
        res_beta$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_m, diff, iter_time_secs
      ))
    }

    if (diff <= tune_tol) {
      # if (one_more_fit || !final_fit) {
      best_fit <- res_beta$fit
      last_output <- res_beta
      break
      # }
      # one_more_fit <- TRUE
    }
    #---------------------------------------------------------------------------
    # --- Step B: Tune Gamma (given current Beta) ---
    t_start_iter <- Sys.time()
    res_gamma <- imr_tune_lasso(
      data = data, grid = grid, target = "gamma", fixed_other_lasso = cur_beta,
      final_fit = final_fit, # if (one_more_fit) final_fit else FALSE,
      use_warm_in_final = use_warm_in_final, fast_nuclear = fast_nuclear,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose - 1,
      n_cores = n_cores, seed = seed
    )
    cur_gamma <- res_gamma$params$lambda_gamma
    iter_time_secs <- as.numeric(difftime(Sys.time(), t_start_iter, units = "secs"))
    #--------------------------------------------------
    # track parameter squared difference for convergence.
    # diff_gamma <- (cur_gamma - old_gamma)^2
    # diff <- diff_beta + diff_gamma
    diff <- abs(res_gamma$params$verror - old_verror) / old_verror
    old_verror <- res_gamma$params$verror
    # track history for debugging/plotting
    # res_gamma$history$step <- 1
    res_gamma$history$iter <- iter + 0.5
    all_history <- rbind(all_history, res_gamma$history)
    # track best performance rows as well
    # res_gamma$params$step <- 1
    res_gamma$params$iter <- iter + 0.5
    res_gamma$params$diff <- diff
    all_params <- rbind(all_params, res_gamma$params)
    # convergence check
    if (verbose > 0) {
      message(sprintf(
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | nuclear: %.4f | Diff: %.6f | Time: %.2fs",
        res_gamma$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_m, diff, iter_time_secs
      ))
    }

    if (diff <= tune_tol) {
      # if (one_more_fit || !final_fit) {
      best_fit <- res_gamma$fit
      last_output <- res_gamma
      break
      # }
      # one_more_fit <- TRUE
    }
    #------------------------------
    # if you reach the final iteration, do one more.
    if (iter >= tune_maxit) {
      # if (final_fit) {
      #  one_more_fit <- TRUE
      # } else {
      keep_iterating <- FALSE
      # }
    }
    iter <- iter + 1
  }
  #---------------------------------------------------------------------------
  t_total <- difftime(Sys.time(), t_start_global)
  if (verbose > 0) {
    if (iter >= tune_maxit) {
      message(">> Alternating tuning reached maxed iterations. <<")
    } else {
      message(">> Alternating tuning converged. <<")
    }
    message(sprintf("\nTotal Tuning Time: %s", format(t_total)))
  }
  #---------------------------------------------------------------------------
  return(list(
    all_params = all_params,
    history = all_history,
    fit = last_output$fit, # if (final_fit) best_fit else NULL,
    params = last_output$params,
    time_secs = as.numeric(t_total,"secs")
  ))
}
