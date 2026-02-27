#----------------------------------------------------------
imr_tune_laplace <- function(data,
                             grid,
                             lambda_beta = 0,
                             lambda_gamma = 0,
                             intercept_row = FALSE,
                             intercept_col = FALSE,
                             shared_beta = FALSE,
                             shared_gamma = FALSE,
                             final_fit = TRUE,
                             convergence = IMR::imr_convergence(),
                             error_function = IMR::error_metrics$rmse,
                             warm_start = NULL,
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
    is.Incomplete(data$y_valid),
    is.Incomplete(data$y_train)
  )
  stopifnot(is.numeric(grid$laplace$max))
  stopifnot(
    grid$rank$min > 0,
    grid$rank$max >= grid$rank$min,
    grid$rank$step >= 0,
    grid$laplace$length >= 1,
    grid$laplace$max >= grid$laplace$min
  )
  #---------------------------------------------------
  # indices
  if (grid$laplace$max <= 0) grid$laplace$max <- 1e-4
  if (grid$laplace$min <= 0) grid$laplace$min <- 1e-6
  lambda_seq <- exp(seq(log(grid$laplace$max),
    log(grid$laplace$min),
    length.out = grid$laplace$length
  ))
  virow <- data$y_valid@i
  vpcol <- data$y_valid@p
  reference <- data$y_valid@x
  # initial values
  mfit <- warm_start
  rank_max <- grid$rank$min
  history <- data.frame(
    lambda_laplace = lambda_seq,
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
    mfit <- IMR::imr_fit(data,
      rank = rank_max,
      lambda_m = lambda_seq[i],
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_beta = shared_beta,
      shared_gamma = shared_gamma,
      convergence = convergence,
      training = TRUE,
      warm_start = mfit
    )
    # compute validation error
    vestim <- IMR::reconstruct_partial(mfit, data, data$valid_mask)@x
    verror <- error_function(vestim, reference)
    # compute rank
    rank_out <- sum(mfit$coefficients$d > 1e-4)
    # verbose
    if (verbose >= 2) {
      message(sprintf(
        "%2d lambda_laplace=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
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
    if (verror <= best_verror) {
      best_fit_obj <- mfit
      best_verror <- verror
      best_params <- history[i, ]
      no_improve_count <- 0
    } else {
      no_improve_count <- no_improve_count + 1
    }

    if (no_improve_count >= grid$laplace$streaks) {
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
      "Best fit: lambda_laplace=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
      best_params$lambda_seq[i],
      best_params$rank_in,
      best_params$rank_out,
      best_verror,
      best_fit_obj$meta$n_iter
    ))
  }
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    best_fit_obj <- IMR::imr_fit(data,
      rank = best_params$rank_in,
      lambda_m = best_params$lambda_laplace,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_beta = shared_beta,
      shared_gamma = shared_gamma,
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
imr_tune_lasso <- function(data,
                           grid,
                           target = c("beta", "gamma"),
                           fixed_other_lasso = 0,
                           intercept_row = FALSE,
                           intercept_col = FALSE,
                           shared_beta = FALSE,
                           shared_gamma = FALSE,
                           final_fit = TRUE,
                           use_warm_in_final = TRUE,
                           convergence = IMR::imr_convergence(),
                           error_function = IMR::error_metrics$rmse,
                           warm_start = NULL,
                           verbose = 1,
                           n_cores = 4,
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
    is.Incomplete(data$y_valid),
    is.Incomplete(data$y_train)
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
  # if(!final_fit) use_warm_in_final = FALSE # not needed.
  lambda_seq <- seq(lambda_obj$min, lambda_obj$max, length.out = lambda_obj$length)
  results_list <- parallel::mclapply(seq_along(lambda_seq), function(i) {
    lambda <- lambda_seq[i]
    if (verbose >= 2) {
      cat(sprintf("Worker started: %s = %.4f\n", target, lambda))
    }
    # run full tune_laplace
    laplace_res <- imr_tune_laplace(
      data = data,
      grid = grid,
      lambda_beta = if (is_beta) lambda else fixed_other_lasso,
      lambda_gamma = if (is_beta) fixed_other_lasso else lambda,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_beta = shared_beta,
      shared_gamma = shared_gamma,
      final_fit = FALSE,
      convergence = convergence,
      error_function = error_function,
      warm_start = warm_start,
      verbose = 0,
      seed = seed
    )
    # return best results
    best_inner <- laplace_res$params
    best_inner[[paste0("lambda_", target)]] <- lambda
    best_inner[[paste0("lambda_", if (is_beta) "gamma" else "beta")]] <- fixed_other_lasso
    # we also return the full history for diagnostics/plots
    history <- laplace_res$history
    history[[paste0("lambda_", target)]] <- lambda
    history[[paste0("lambda_", if (is_beta) "gamma" else "beta")]] <- fixed_other_lasso

    out <- list(
      best_inner = best_inner,
      history = history
    )
    if (use_warm_in_final) {
      out$fit <- laplace_res$fit
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
      "Best parameters: %s: %.5f | Laplace: %.4f | Error: %.5f",
      target, best_params[[paste0("lambda_", target)]],
      best_params$lambda_laplace, best_params$verror
    ))
  }
  #---------------------------------------------
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    if (use_warm_in_final) {
      warm_start <- results_list[[best_idx]]$fit
    }
    best_fit_obj <- IMR::imr_fit(data,
      rank = best_params$rank_in,
      lambda_m = best_params$lambda_laplace,
      lambda_beta = best_params$lambda_beta,
      lambda_gamma = best_params$lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_beta = shared_beta,
      shared_gamma = shared_gamma,
      convergence = convergence,
      warm_start = warm_start
    )
  } else {
    best_fit_obj <- NULL
  }

  list(fit = best_fit_obj, params = best_params, history = history)
}
# =============================================================================================
#' @export
imr_tune <- function(data,
                     grid,
                     default_lambda_beta = 0,
                     default_lambda_gamma = 0,
                     intercept_row = FALSE,
                     intercept_col = FALSE,
                     shared_beta = FALSE,
                     shared_gamma = FALSE,
                     final_fit = TRUE,
                     use_warm_in_final = TRUE,
                     convergence = IMR::imr_convergence(),
                     error_function = IMR::error_metrics$rmse,
                     warm_start = NULL,
                     verbose = 1,
                     n_cores = 4,
                     seed = NULL,
                     tune_maxit = 10,
                     tune_tol = 1e-4) {
  if (!is.null(seed) && is.numeric(seed)) set.seed(seed)
  stopifnot(inherits(data, "imr_data"), inherits(grid, "imr_tune_grid"))
  t_start_global <- Sys.time()
  # --- Determine which parameters to tune
  tune_beta <- data$meta$has_X && grid$beta$length > 1
  tune_gamma <- data$meta$has_Z && grid$gamma$length > 1

  #----------------------------------------------------------
  # Scenario 1: Tune Laplace Only
  #------------------------------------------------------------
  if (!tune_beta && !tune_gamma) {
    if (verbose > 0) message("Tuning  Laplace (M) only...")
    out_obj <- imr_tune_laplace(
      data = data, grid = grid,
      lambda_beta = default_lambda_beta, lambda_gamma = default_lambda_gamma,
      intercept_row = intercept_row, intercept_col = intercept_col,
      shared_beta = shared_beta, shared_gamma = shared_gamma,
      final_fit = final_fit, convergence = convergence,
      error_function = error_function, warm_start = warm_start,
      verbose = verbose, seed = seed
    )
    t_total <- round(difftime(Sys.time(), t_start_global), 2)
    if (verbose > 0) {
      message(sprintf("\nTotal Tuning Time: %s", format(t_total)))
    }
    out_obj$time <- t_total
    return(out_obj)
  }
  #---------------------------------------------------------------------
  # Scenario 2: Tuning Laplace + one of lambda_beta or lambda_gamma
  #-----------------------------------------------------------------
  if (tune_beta != tune_gamma) {
    target <- if (tune_beta) "beta" else "gamma"
    fixed_other <- if (tune_beta) default_lambda_gamma else default_lambda_beta

    if (verbose > 0) message(sprintf("Tuning %s + Laplace...", (target)))
    out_obj <- imr_tune_lasso(
      data = data, grid = grid, target = target, fixed_other_lasso = fixed_other,
      intercept_row = intercept_row, intercept_col = intercept_col,
      shared_beta = shared_beta, shared_gamma = shared_gamma,
      final_fit = final_fit, use_warm_in_final = use_warm_in_final,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose,
      n_cores = n_cores, seed = seed
    )
    t_total <- difftime(Sys.time(), t_start_global)
    if (verbose > 0) {
      message(sprintf("\nTotal Tuning Time: %s", format(t_total)))
    }
    out_obj$time <- t_total
    return(out_obj)
  }
  #-----------------------------------------------------------------------
  # Scenario 3: Tuning all 3 parameters. An iterative method will be used.
  #------------------------------------------------------------------------
  if (verbose > 0) message("Scenario 3: Alternating optimization for all 3 parameters...")

  # Initialize
  cur_gamma <- grid$gamma$max
  cur_beta <- diff_beta <- diff_gamma <- Inf
  all_history <- all_params <- data.frame()
  one_more_fit <- FALSE # this is an extra for for the purpose of final_fit=TRUE
  keep_iterating <- TRUE
  iter <- 1
  while (keep_iterating) {
    if (verbose > 0) message(sprintf("\n--- Tune Iteration %d ---", iter))

    old_beta <- cur_beta
    old_gamma <- cur_gamma

    # --- Step A: Tune Beta (given current Gamma) ---
    t_start_iter <- Sys.time()
    res_beta <- imr_tune_lasso(
      data = data, grid = grid, target = "beta", fixed_other_lasso = cur_gamma,
      intercept_row = intercept_row, intercept_col = intercept_col,
      shared_beta = shared_beta, shared_gamma = shared_gamma,
      final_fit = if (one_more_fit) final_fit else FALSE,
      use_warm_in_final = use_warm_in_final,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose - 1,
      n_cores = n_cores, seed = seed
    )
    cur_beta <- res_beta$params$lambda_beta
    iter_time_secs <- as.numeric(difftime(Sys.time(), t_start_iter, units = "secs"))
    #---------------------------------------------------------------------------
    # track parameter squared difference for convergence.
    diff_beta <- (cur_beta - old_beta)^2
    diff <- diff_beta + diff_gamma
    # track history for debugging/plotting
    res_beta$history$step <- 1
    res_beta$history$iter <- iter
    all_history <- rbind(all_history, res_beta$history)
    # track best performance rows as well
    res_beta$params$step <- 1
    res_beta$params$iter <- iter
    res_beta$params$diff <- diff
    all_params <- rbind(all_params, res_beta$params)
    # convergence check
    if (verbose > 0 && iter > 1) {
      message(sprintf(
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | Laplace: %.4f | Diff: %.6f | Time: %.2fs",
        res_beta$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_laplace, diff, iter_time_secs
      ))
    }

    if (diff < tune_tol) {
      if (one_more_fit || !final_fit) {
        best_fit <- res_beta$fit
        break
      }
      one_more_fit <- TRUE
    }
    #---------------------------------------------------------------------------
    # --- Step B: Tune Gamma (given current Beta) ---
    t_start_iter <- Sys.time()
    res_gamma <- imr_tune_lasso(
      data = data, grid = grid, target = "gamma", fixed_other_lasso = cur_beta,
      intercept_row = intercept_row, intercept_col = intercept_col,
      shared_beta = shared_beta, shared_gamma = shared_gamma,
      final_fit = if (one_more_fit) final_fit else FALSE,
      use_warm_in_final = use_warm_in_final,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose - 1,
      n_cores = n_cores, seed = seed
    )
    cur_gamma <- res_gamma$params$lambda_gamma
    iter_time_secs <- as.numeric(difftime(Sys.time(), t_start_iter, units = "secs"))
    #--------------------------------------------------
    # track parameter squared difference for convergence.
    diff_gamma <- (cur_gamma - old_gamma)^2
    diff <- diff_beta + diff_gamma
    # track history for debugging/plotting
    res_gamma$history$step <- 1
    res_gamma$history$iter <- iter
    all_history <- rbind(all_history, res_gamma$history)
    # track best performance rows as well
    res_gamma$params$step <- 1
    res_gamma$params$iter <- iter
    res_gamma$params$diff <- diff
    all_params <- rbind(all_params, res_gamma$params)
    # convergence check
    if (verbose > 0) {
      message(sprintf(
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | Laplace: %.4f | Diff: %.6f | Time: %.2fs",
        res_gamma$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_laplace, diff, iter_time_secs
      ))
    }

    if (diff < tune_tol) {
      if (one_more_fit || !final_fit) {
        best_fit <- res_gamma$fit
        break
      }
      one_more_fit <- TRUE
    }
    #------------------------------
    # if you reach the final iteration, do one more.
    if (iter >= tune_maxit) {
      if (final_fit) {
        one_more_fit <- TRUE
      } else {
        keep_iterating <- FALSE
      }
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
    params = all_params,
    history = all_history,
    fit = if (final_fit) best_fit else NULL,
    time = t_total
  ))
}


#
#
#
# imr.cv_laplace <- function(
#   data,
#   intercept_row = FALSE,
#   intercept_col = FALSE,
#   hpar = IMR::get_imr_default_hparams(),
#   error_function = IMR:::error_metric$rmse,
#   thresh = 1e-4,
#   maxit = 300,
#   trace = 1,
#   warm_start = NULL,
#   ls_initial = TRUE,
#   shared_information = FALSE,
#   num_cores = 6,
#   final_fit = TRUE,
#   final_thresh = 1e-6,
#   final_maxit = 1000,
#   seed = NULL
# ) {
#   stopifnot(is.Incomplete(data$Y))
#   stopifnot(is.Incomplete(data$y_train))
#   stopifnot(is.Incomplete(data$y_valid))
#   if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
#   if (is.numeric(num_cores) && num_cores > 0) IMR::initialize_parallel_workers(num_cores)
#   #---------------------------------------------------
#
#   # fixed: all. variable: none. number of fits: 1.
#   # fits a single fit and returns [fit, error]; with all lambdas fixed.
#   # this is a single fit where all parameters are fixed but it returns validation error
#   rank_fit_function <- function(r, fdata, hpar, shared_information,
#                                 lambda_laplace,
#                                 intercept_row, intercept_col,
#                                 trace, thresh, maxit,
#                                 ls_initial, fit = NULL,
#                                 error_function = IMR:::error_metric$rmse) {
#     fit <- IMR::imr.fit(
#       Y = fdata$y_train,
#       X = fdata$Xq,
#       Z = fdata$Zq,
#       r = r,
#       lambda_m = lambda_laplace,
#       lambda_beta = hpar$beta$value,
#       lambda_gamma = hpar$gamma$value,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       shared_information = shared_information,
#       Ur = hpar$laplacian_row$U,
#       dr = hpar$laplacian_row$d,
#       Uc = hpar$laplacian_col$U,
#       dc = hpar$laplacian_col$d,
#       warm_start = fit,
#       trace = F,
#       thresh = thresh,
#       maxit = maxit,
#       ls_initial = ls_initial
#     )
#
#
#     vestim <- IMR:::reconstruct_partial(fit, fdata, fdata$y_valid, shared_information)
#     verror <- error_function(fdata$y_valid@x, vestim@x)
#     # verbose
#
#     if (trace >= 3) message(sprintf("rank=%d | err=%.5f ", r, verror))
#     fit$r <- r
#     return(list(fit = fit, error = verror))
#   }
#
#   #---
#   # this function takes a single lambda_laplace and finds the optimal rank. everything else fixed.
#   laplace_cv_lambda_function <- function(lambda_laplace, data, hpar,
#                                          intercept_row, shared_information,
#                                          intercept_col, trace, thresh,
#                                          maxit, ls_initial, fit = NULL, warm_start = NULL) {
#     if (!is.null(data$similarity_rows) && lambda_laplace > 0) {
#       hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row, lambda_laplace)
#     }
#     if (!is.null(data$similarity_cols) && lambda_laplace > 0) {
#       hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col, lambda_laplace)
#     }
#
#     results <- IMR::adaptive_tuner(rank_fit_function,
#       step_sizes = hpar$rank$step_sizes,
#       start_value = hpar$rank$min,
#       end_value = hpar$rank$max,
#       inc_streak_to_stop = hpar$rank$n_streaks,
#       fdata = data,
#       hpar = hpar,
#       lambda_laplace = lambda_laplace,
#       shared_information = shared_information,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       trace = trace,
#       thresh = thresh,
#       .warm_start = warm_start,
#       maxit = maxit,
#       ls_initial = ls_initial
#     )
#     if (trace >= 2) {
#       message(sprintf(
#         " lambda_laplace = %.3f | best rank = %.0f | error = %.5f",
#         lambda_laplace,
#         results$best_parameter,
#         results$best_error
#       ))
#     }
#     results$best_fit$lambda_laplace <- lambda_laplace
#     return(list(fit = results$best_fit, error = results$best_error, history = results$history))
#   }
#
#   #----------------
#   # the following goes over a grid of lambda_laplace and then calls what's above.
#   # this one is supposed to run in parallel!! implement that.
#   # output: results$best_parameter: best lambda_laplace, and error. it also returns results$fit
#   results <- IMR:::parallel_grid_1d_adaptive(
#     param_min = hpar$laplace$min,
#     param_max = hpar$laplace$max,
#     step_sizes = hpar$laplace$step_sizes,
#     f = laplace_cv_lambda_function,
#     .progress = TRUE,
#     .trace = trace >= 4,
#     # .packages = c("IMR"),
#     .seed = if (is.null(seed)) FALSE else seed,
#     data = data,
#     hpar = hpar,
#     shared_information = shared_information,
#     intercept_row = intercept_row,
#     intercept_col = intercept_col,
#     trace = trace,
#     thresh = thresh,
#     maxit = maxit,
#     warm_start = warm_start,
#     ls_initial = ls_initial
#   )
#
#
#   if (trace >= 1) {
#     message(sprintf(
#       "Best lambda_laplace = %.3f | best rank = %.0f | error = %.5f",
#       results$best_fit$lambda_laplace,
#       results$best_fit$r,
#       results$best_error
#     ))
#   }
#
#
#   # (optional) retrain on the full data
#   if (!is.null(data$Y) & final_fit) {
#     if (!is.null(data$similarity_rows) && results$best_fit$lambda_laplace > 0) {
#       hpar$laplacian_row <- IMR::decompose_symmetric_matrix(
#         data$similarity_row,
#         results$best_fit$lambda_laplace
#       )
#     }
#     if (!is.null(data$similarity_cols) && results$best_fit$lambda_laplace > 0) {
#       hpar$laplacian_col <- IMR::decompose_symmetric_matrix(
#         data$similarity_col,
#         results$best_fit$lambda_laplace
#       )
#     }
#     fit <- IMR::imr.fit(
#       Y = data$Y,
#       X = data$Xq,
#       Z = data$Zq,
#       r = results$best_fit$r,
#       lambda_m = results$best_fit$lambda_laplace,
#       lambda_beta = hpar$beta$value,
#       lambda_gamma = hpar$gamma$value,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       shared_information = shared_information,
#       Ur = hpar$laplacian_row$U,
#       dr = hpar$laplacian_row$d,
#       Uc = hpar$laplacian_col$U,
#       dc = hpar$laplacian_col$d,
#       warm_start = results$best_fit,
#       trace = trace >= 3,
#       thresh = final_thresh,
#       maxit = final_maxit,
#       ls_initial = ls_initial
#     )
#     results$best_fit[names(fit)] <- fit
#   }
#
#   # record final hyper-parameters and return
#   results$lambda_beta <- hpar$beta$value
#   results$lambda_gamma <- hpar$gamma$value
#   results$fit <- results$best_fit
#   results$fit$params <-
#     list(
#       lambda_beta = results$lambda_beta,
#       lambda_gamma = results$lambda_gamma,
#       lambda_laplace = results$fit$lambda_laplace,
#       rank = results$fit$r
#     )
#   return(results)
# }

#-----------------------------------------------------------
#' @export
imr.cv <- function(
  data,
  intercept_row = FALSE,
  intercept_col = FALSE,
  hpar = IMR::get_imr_default_hparams(),
  error_function = IMR:::error_metric$rmse,
  thresh = 1e-6,
  maxit = 500,
  trace = 0,
  ls_initial = FALSE,
  shared_information = FALSE,
  num_cores = num_cores,
  warm_start = NULL,
  seed = NULL,
  # fast.cv = FALSE,
  layer_1_parallel = TRUE,
  separate_tuning = FALSE
) {
  # if (fast.cv) {
  #   return(
  #     IMR:::nlrr.cv(
  #       data = data,
  #       intercept_row = intercept_row,
  #       intercept_col = intercept_col,
  #       lambda_beta = lambda_beta,
  #       lambda_gamma = lambda_gamma,
  #       hpar = hpar,
  #       error_function = error_function,
  #       thresh = thresh,
  #       maxit = maxit,
  #       verbose = verbose,
  #       seed = seed,
  #       ls_initial = ls_initial
  #     )
  #   )
  # }

  #-------------------
  stopifnot(is.Incomplete(data$Y))
  stopifnot(is.Incomplete(data$y_train))
  stopifnot(is.Incomplete(data$y_valid))
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  #-------------------------------
  # set flags
  beta_flag <- !(is.null(data$Xq))
  gamma_flag <- !(is.null(data$Zq))
  tune_beta <- is.null(hpar$beta$max) || hpar$beta$max != hpar$beta$value
  tune_gamma <- is.null(hpar$gamma$max) || hpar$gamma$max != hpar$gamma$value

  # if neither beta or gamma are provided
  if (!(beta_flag | gamma_flag)) {
    return(IMR:::imr.cv_laplace(
      data = data,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      hpar = hpar,
      error_function = error_function,
      thresh = thresh,
      trace = trace,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed,
      num_cores = num_cores
    ))
  }

  # obtain upperbounds to the lambda hyperparameters
  if (beta_flag & is.null(hpar$beta$max)) {
    hpar$beta$max <- IMR::get_lambda_lasso_max(
      y_train = data$y_train,
      X = data$Xq,
      # y_valid = data$y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = trace
    )
    if (is.function(hpar$beta$step_sizes)) {
      hpar$beta$step_sizes <- hpar$beta$step_sizes(hpar$beta$min, hpar$beta$max)
    }
  }
  if (gamma_flag & is.null(hpar$gamma$max)) {
    hpar$gamma$max <- IMR::get_lambda_lasso_max(
      y_train = data$y_train,
      Z = data$Zq,
      # y_valid = data$y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = trace
    )
    if (is.function(hpar$gamma$step_sizes)) {
      hpar$gamma$step_sizes <- hpar$gamma$step_sizes(hpar$gamma$min, hpar$gamma$max)
    }
  }


  #---------------------------
  # parallel setup
  # the following function takes
  single_fit <- function(parameter, type = "rows", data, intercept_row, intercept_col,
                         shared_information,
                         hpar, error_function, thresh, trace, maxit, ls_initial,
                         seed, num_cores, fit = NULL) {
    if (type == "rows") {
      hpar$beta$value <- parameter
    } else {
      hpar$gamma$value <- parameter
    }
    #---------------------------------------
    results <- IMR:::imr.cv_laplace(
      data = data,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      hpar = hpar,
      shared_information = shared_information,
      error_function = error_function,
      thresh = thresh,
      trace = trace - 2,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed,
      num_cores = num_cores,
      # warm_start = fit,
      final_fit = FALSE
    )
    if (trace >= 2) {
      message(sprintf(
        paste0(
          "lambda_beta = %.5f | lambda_gamma = %.5f | ",
          "best lambda_laplace = %.3f | best rank = %.0f | error = %.6f"
        ),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$r,
        results$best_error
      ))
    }
    results$best_fit$lambda_beta <- hpar$beta$value
    results$best_fit$lambda_gamma <- hpar$gamma$value
    return(list(fit = results$best_fit, error = results$best_error))
  }
  #----------------------------------------------------------------------------
  # if (separate_tuning & gamma_flag & beta_flag) {
  message("Fitting lambda_beta and lambda_gamma separately...")
  # if separate then tune beta first followed by gamma. Meanwhile,
  # keep lambda_gamma at 10% of its maximum.
  if (beta_flag) {
    if (layer_1_parallel) {
      results <- IMR:::parallel_grid_1d_adaptive(
        param_min = hpar$beta$min,
        param_max = hpar$beta$max,
        step_sizes = hpar$beta$step_sizes,
        f = single_fit,
        .progress = TRUE,
        .trace = trace >= 3,
        # .packages = c("IMR"),
        .seed = if (is.null(seed)) FALSE else seed,
        type = "rows",
        data = data,
        intercept_row = intercept_row,
        intercept_col = intercept_col,
        hpar = hpar,
        shared_information = shared_information,
        error_function = error_function,
        thresh = thresh,
        trace = trace,
        maxit = maxit,
        ls_initial = ls_initial,
        num_cores = 0,
        seed = seed
        # .warm_start = warm_start
      )
    } else {
      results <- IMR::adaptive_tuner(
        single_fit,
        step_sizes = hpar$beta$step_sizes,
        end_value = hpar$beta$max,
        start_value = 0,
        inc_streak_to_stop = hpar$beta$n_streaks,
        type = "rows",
        data = data,
        intercept_row = intercept_row,
        intercept_col = intercept_col,
        hpar = hpar,
        shared_information = shared_information,
        error_function = error_function,
        thresh = thresh,
        trace = trace,
        maxit = maxit,
        ls_initial = ls_initial,
        seed = seed,
        .warm_start = warm_start,
        num_cores = num_cores # ,
        # fit = NULL
      )
    }
    hpar$beta$value <- results$best_fit$lambda_beta
    if (trace >= 1) {
      message(sprintf(
        paste0(
          "best lambda_beta = %.5f | lambda_gamma = %.5f | ",
          "best lambda_laplace = %.4f | best rank = %.0f | error = %.5f"
        ),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$r,
        results$best_error
      ))
    }
  }
  if (gamma_flag) {
    # we now do the same to gamma
    results <- IMR::adaptive_tuner(
      single_fit,
      step_sizes = hpar$gamma$step_sizes,
      end_value = hpar$gamma$max,
      start_value = 0,
      inc_streak_to_stop = hpar$gamma$n_streaks,
      type = "cols",
      data = data,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      hpar = hpar,
      shared_information = shared_information,
      error_function = error_function,
      thresh = thresh,
      trace = trace,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed,
      num_cores = num_cores,
      .warm_start = results$best_fit
    )

    if (trace >= 1) {
      message(sprintf(
        paste0(
          "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
          "best lambda_laplace = %.2f |  best rank = %.0f | error = %.5f"
        ),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$r,
        results$best_error
      ))
    }
  }
  if (!is.null(data$Y)) {
    if (!is.null(data$similarity_rows) && results$best_fit$lambda_laplace > 0) {
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(
        data$similarity_row,
        results$best_fit$lambda_laplace
      )
    }
    if (!is.null(data$similarity_cols) && results$best_fit$lambda_laplace > 0) {
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(
        data$similarity_col,
        results$best_fit$lambda_laplace
      )
    }

    fit <- IMR::imr.fit(
      Y = data$Y,
      X = data$Xq,
      Z = data$Zq,
      r = results$best_fit$r,
      lambda_m = results$best_fit$lambda_laplace,
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_information = shared_information,
      Ur = hpar$laplacian_row$U,
      dr = hpar$laplacian_row$d,
      Uc = hpar$laplacian_col$U,
      dc = hpar$laplacian_col$d,
      warm_start = results$best_fit,
      trace = trace >= 3,
      thresh = thresh,
      maxit = maxit,
      ls_initial = ls_initial
    )
    results$best_fit[names(fit)] <- fit
  }
  return(results)
}
#
#     results <- parallel_grid(grid, IMR::imr.cv_M,
#       "list",
#       .packages = "IMR",
#       .progress = TRUE,
#       .seed = seed,
#       y_train = data$y_train,
#       y_valid = data$y_valid,
#       X = data$Xq,
#       # Z = data$Zq,
#       Y_full = NULL,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       hpar = hpar,
#       error_function = error_function,
#       thresh = thresh,
#       maxit = maxit,
#       trace = verbose >= 1,
#       ls_initial = ls_initial,
#       seed = seed
#     )
#
#     # Select the best fit
#     errors <- vapply(results, `[[`, numeric(1), "error")
#     best_idx <- which.min(errors)
#     best_fit <- results[[best_idx]]
#     best_fit$fit$gamma <- matrix(0, nrow(data$Y), ncol(data$Zq))
#     #----
#     # we now tune lambda gamma
#     grid <- list(
#       lambda_gamma = lambda_gamma_grid,
#       lambda_beta = best_fit$lambda_beta
#     )
#     results <- parallel_grid(grid, IMR::imr.cv_M,
#       "list",
#       .packages = "IMR",
#       .progress = TRUE,
#       .seed = seed,
#       y_train = data$y_train,
#       y_valid = data$y_valid,
#       X = data$Xq,
#       Z = data$Zq,
#       Y_full = data$Y,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       hpar = hpar,
#       error_function = error_function,
#       thresh = thresh,
#       maxit = maxit,
#       trace = verbose >= 1,
#       ls_initial = ls_initial,
#       old_fit = best_fit$fit,
#       seed = seed
#     )
#
#
#     # Select the best fit
#     errors <- vapply(results, `[[`, numeric(1), "error")
#     best_idx <- which.min(errors)
#     best_fit <- results[[best_idx]]
#   } else {
#     message("Fitting lambda_beta and lambda_gamma simultaneously ...")
#     grid <- list(
#       lambda_beta = lambda_beta_grid,
#       lambda_gamma = lambda_gamma_grid
#     )
#
#     results <- parallel_grid(grid, IMR::imr.cv_M,
#       "list",
#       .packages = "IMR",
#       .progress = TRUE,
#       .seed = seed,
#       y_train = data$y_train,
#       y_valid = data$y_valid,
#       X = data$Xq,
#       Z = data$Zq,
#       Y_full = data$Y,
#       intercept_row = intercept_row,
#       intercept_col = intercept_col,
#       hpar = hpar,
#       error_function = error_function,
#       thresh = thresh,
#       maxit = maxit,
#       trace = verbose >= 1,
#       ls_initial = ls_initial,
#       seed = seed
#     )
#
#
#     # Select the best fit
#     errors <- vapply(results, `[[`, numeric(1), "error")
#     best_idx <- which.min(errors)
#     best_fit <- results[[best_idx]]
# }
# #--------------------
# # message >>
# if (verbose >= 1) {
#   for (res in results) {
#     message(sprintf(
#       "<< lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
#       res$lambda_beta,
#       sum(res$fit$beta == 0) / length(res$fit$beta),
#       res$lambda_gamma,
#       sum(res$fit$gamma == 0) / length(res$fit$gamma),
#       res$error,
#       res$loop_size,
#       res$rank_M,
#       res$lambda_m
#     ))
#   }
#   message(sprintf(
#     "<< Best fit >> lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
#     best_fit$lambda_beta,
#     sum(best_fit$fit$beta == 0) / length(best_fit$fit$beta),
#     best_fit$lambda_gamma,
#     sum(best_fit$fit$gamma == 0) / length(best_fit$fit$gamma),
#     best_fit$error,
#     best_fit$loop_size,
#     best_fit$rank_M,
#     best_fit$lambda_m
#   ))
#   best_fit$init_hparams <- hpar
# }
# rm(results)
# return(best_fit)
