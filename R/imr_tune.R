#----------------------------------------------------------
imr_tune_laplace_fast <- function(data,
                                  grid,
                                  lambda_beta = 0,
                                  lambda_gamma = 0,
                                  final_fit = TRUE,
                                  convergence = IMR::imr_convergence(),
                                  error_function = IMR::error_metrics$rmse,
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
  if (log_grid) {
    lambda_seq <- exp(seq(log(grid$laplace$max),
      log(grid$laplace$min),
      length.out = grid$laplace$length
    ))
  } else {
    lambda_seq <- seq(grid$laplace$max, grid$laplace$min,
      length.out = grid$laplace$length
    )
  }
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
    completely_sparse <- rank_out == 0
    if (verror < best_verror) {
      best_fit_obj <- mfit
      best_verror <- verror
      best_params <- history[i, ]
      no_improve_count <- 0
    } else if (!completely_sparse) {
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
      best_params$lambda_laplace,
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
imr_tune_laplace_slow <- function(data,
                                  grid,
                                  lambda_beta = 0,
                                  lambda_gamma = 0,
                                  final_fit = TRUE,
                                  convergence = IMR::imr_convergence(),
                                  error_function = IMR::error_metrics$rmse,
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
  # training grids
  if (grid$laplace$max <= 0) grid$laplace$max <- 1e-4
  if (grid$laplace$min <= 0) grid$laplace$min <- 1e-6
  if (log_grid) {
    lambda_seq <- exp(seq(log(grid$laplace$max),
                          log(grid$laplace$min),
                          length.out = grid$laplace$length
    ))
  } else {
    lambda_seq <- seq(grid$laplace$max, grid$laplace$min,
                      length.out = grid$laplace$length
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
    lambda_laplace = lambda_seq
  )
  history <- history[, c("lambda_laplace", "rank")]
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

      mfit <- IMR::imr_fit(data,
        rank = current_rank,
        lambda_m = current_lambda,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
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

    if (no_improve_count_1 >= grid$laplace$streaks) {
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
      "Best fit: lambda_laplace=%.6f | rank_in=%d => rank_out=%d | verr=%.5f | fit_iter=%d",
      best_params_1$lambda_laplace,
      best_params_1$rank_in,
      best_params_1$rank_out,
      best_verror_1,
      best_fit_obj_1$meta$n_iter
    ))
  }
  # (optional) final fit on full dataset
  if (final_fit) {
    if (verbose > 0) message("Fitting final model on full dataset...")
    best_fit_obj_1 <- IMR::imr_fit(data,
      rank = best_params_1$rank_in,
      lambda_m = best_params_1$lambda_laplace,
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

imr_tune_lasso <- function(data,
                           grid,
                           target = c("beta", "gamma"),
                           fixed_other_lasso = 0,
                           final_fit = TRUE,
                           use_warm_in_final = TRUE,
                           convergence = IMR::imr_convergence(),
                           error_function = IMR::error_metrics$rmse,
                           warm_start = NULL,
                           verbose = 1,
                           n_cores = 4,
                           fast_laplace = TRUE,
                           laplace_log_scale = TRUE,
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
  laplace_function <- if (fast_laplace) imr_tune_laplace_fast else imr_tune_laplace_slow
  # if(!final_fit) use_warm_in_final = FALSE # not needed.
  lambda_seq <- seq(lambda_obj$min, lambda_obj$max, length.out = lambda_obj$length)
  results_list <- parallel::mclapply(seq_along(lambda_seq), function(i) {
    lambda <- lambda_seq[i]
    if (verbose >= 2) {
      cat(sprintf("Worker started: %s = %.4f\n", target, lambda))
    }
    # run full tune_laplace
    laplace_res <- laplace_function(
      data = data,
      grid = grid,
      lambda_beta = if (is_beta) lambda else fixed_other_lasso,
      lambda_gamma = if (is_beta) fixed_other_lasso else lambda,
      final_fit = FALSE,
      convergence = convergence,
      error_function = error_function,
      warm_start = warm_start,
      log_grid = laplace_log_scale,
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
                     final_fit = TRUE,
                     use_warm_in_final = TRUE,
                     fast_laplace = TRUE,
                     convergence = IMR::imr_convergence(),
                     error_function = IMR::error_metrics$rmse,
                     warm_start = NULL,
                     verbose = 1,
                     n_cores = 4,
                     seed = NULL,
                     laplace_log_scale = TRUE,
                     tune_maxit = 10,
                     tune_tol = 1e-4) {
  if (!is.null(seed) && is.numeric(seed)) set.seed(seed)
  stopifnot(inherits(data, "imr_data"), inherits(grid, "imr_tune_grid"))
  t_start_global <- Sys.time()
  # --- Determine which parameters to tune
  tune_beta <- data$model$row_covariates && data$meta$has_X && grid$beta$length > 1
  tune_gamma <- data$model$col_covariates && data$meta$has_Z && grid$gamma$length > 1

  #----------------------------------------------------------
  # Scenario 1: Tune Laplace Only
  #------------------------------------------------------------
  if (!tune_beta && !tune_gamma) {
    if (verbose > 0) message("Tuning  Laplace (M) only...")
    laplace_function <- if (fast_laplace) imr_tune_laplace_fast else imr_tune_laplace_slow
    out_obj <- laplace_function(
      data = data, grid = grid,
      lambda_beta = default_lambda_beta, lambda_gamma = default_lambda_gamma,
      final_fit = final_fit, convergence = convergence,
      error_function = error_function, warm_start = warm_start,
      log_grid = laplace_log_scale,
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
  # Scenario 2: Tuning Laplace + one of lambda_beta or lambda_gamma
  #-----------------------------------------------------------------
  if (tune_beta != tune_gamma) {
    target <- if (tune_beta) "beta" else "gamma"
    fixed_other <- if (tune_beta) default_lambda_gamma else default_lambda_beta

    if (verbose > 0) message(sprintf("Tuning %s + Laplace...", (target)))
    out_obj <- imr_tune_lasso(
      data = data, grid = grid, target = target, fixed_other_lasso = fixed_other,
      final_fit = final_fit, use_warm_in_final = use_warm_in_final,
      convergence = convergence, error_function = error_function,
      warm_start = warm_start, verbose = verbose, fast_laplace = fast_laplace,
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
      use_warm_in_final = use_warm_in_final, fast_laplace = fast_laplace,
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
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | Laplace: %.4f | Diff: %.6f | Time: %.2fs",
        res_beta$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_laplace, diff, iter_time_secs
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
      use_warm_in_final = use_warm_in_final, fast_laplace = fast_laplace,
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
        "verror: %.4f | Beta: %.4f | Gamma: %.4f | Laplace: %.4f | Diff: %.6f | Time: %.2fs",
        res_gamma$params$verror, cur_beta, cur_gamma, res_gamma$params$lambda_laplace, diff, iter_time_secs
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
