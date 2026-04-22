
#----------------------------------------------------------
#' @export
imr.cv_laplace <- function(
  data,
  row_intercept = FALSE,
  col_intercept = FALSE,
  hpar = IMR::get_imr_default_hparams(),
  error_function = IMR:::error_metric$rmse,
  thresh = 1e-6,
  maxit = 300,
  trace = 1,
  warm_start = NULL,
  ls_initial = FALSE,
  shared_information = FALSE,
  num_cores = 4,
  final_fit = TRUE,
  seed = NULL
) {
  stopifnot(is_incomplete(data$Y))
  stopifnot(is_incomplete(data$y_train))
  stopifnot(is_incomplete(data$y_valid))
  if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
  if( is.numeric(num_cores) && num_cores > 0) IMR::initialize_parallel_workers(num_cores)
  #---------------------------------------------------

  # fixed: all. variable: none. number of fits: 1.
  # fits a single "rank" and returns [fit, error]; with all lambdas fixed.
  # this is a single fit where all parameters are fixed but it returns validation error
  rank_fit_function <- function(r, data, hpar, shared_information,
                                lambda_laplace,
                                row_intercept, col_intercept,
                                trace, thresh, maxit,
                                ls_initial, fit=NULL) {
    if(shared_information){

    }else{
      fit <- IMR::imr.fit(
        Y = data$y_train,
        X = data$Xq,
        Z = data$Zq,
        r = r,
        lambda_m = lambda_laplace,
        lambda_beta = hpar$beta$value,
        lambda_gamma = hpar$gamma$value,
        row_intercept = row_intercept,
        col_intercept = col_intercept,
        Ur = hpar$laplacian_row$U,
        dr = hpar$laplacian_row$d,
        Uc = hpar$laplacian_col$U,
        dc = hpar$laplacian_col$d,
        warm_start = fit,
        trace = F,
        thresh = thresh,
        maxit = maxit,
        ls_initial = ls_initial
      )
    }
    vestim <- IMR:::reconstruct_partial(fit, data, data$y_valid)
    verror <- error_function(data$y_valid@x, vestim@x)
    # verbose

    if (trace >= 4) message(sprintf("rank=%d | err=%.5f ", r, verror))
    fit$r = r
    return(list(fit = fit, error = verror))
  }

  #---


  #--- we fit the function above twice, once for laplace rows and once for columns
  #---
  if (is.null(hpar$laplacian_row$U) || is.null(hpar$laplacian_col$U) )
    stop("Laplace matrices must be initialized")

  #  variable: alpha, lambda_laplace, rank. Fixed: lambda_beta/gamma. Supposed to be parallel!!
  # we run on a grid of (alpha, lambda_laplace)
  #' this function takes a single alpha and finds the optimal rank. everything else fixed.
  laplace_cv_alpha_function <- function(alpha, data, hpar, lambda_laplace,
                                        row_intercept,
                                        col_intercept, trace, thresh,
                                        maxit, ls_initial, fit = NULL, warm_start = NULL) {
    if(!is.null(data$similarity_rows) && lambda_laplace > 0)
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row, lambda_laplace * alpha)
    if(!is.null(data$similarity_cols) && lambda_laplace > 0)
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col, lambda_laplace * (1 - alpha))

    results <- IMR::adaptive_tuner(rank_fit_function,
      step_sizes = hpar$rank$step_sizes,
      start_value = hpar$rank$min,
      end_value = hpar$rank$max,
      inc_streak_to_stop = hpar$rank$n_streaks,
      data = data,
      hpar = hpar,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      trace = trace,
      thresh = thresh,
      .warm_start = warm_start,
      maxit = maxit,
      ls_initial = ls_initial
    )
    if (trace >= 3) {
      message(sprintf(
        " lambda_laplace = %.3f | alpha = %.3f | best rank = %.0f | error = %.5f",
        lambda_laplace,
        alpha,
        results$best_parameter,
        results$best_error
      ))
    }
    results$best_fit$alpha = alpha
    results$best_fit$lambda_laplace = lambda_laplace
    return(list(fit = results$best_fit, error = results$best_error))
  }

  #-- layer two >>
  #' This function takes a single lambda_laplace and finds the optimal (alpha, rank) combination
  #' for a fixed lambda_beta, lambda_gamma
  laplace_cv_lambda_function <- function(lambda_laplace, data, hpar,
                                         row_intercept,
                                         col_intercept, trace, thresh,
                                         maxit, ls_initial, fit = NULL,
                                         warm_start = NULL) {
    results <- IMR::adaptive_tuner(laplace_cv_alpha_function,
      step_sizes = hpar$laplace$alpha_step_sizes,
      start_value = if(lambda_laplace > 0) hpar$laplace$alpha_min else 0,
      end_value = if(lambda_laplace > 0) hpar$laplace$alpha_max else 0,
      inc_streak_to_stop = hpar$laplace$n_streaks,
      lambda_laplace = lambda_laplace,
      data = data,
      hpar = hpar,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      trace = trace,
      thresh = thresh,
      warm_start = warm_start,
      maxit = maxit,
      ls_initial = ls_initial
    )
    if (trace >= 2) {
      message(sprintf(
        "lambda_laplace = %.3f | best alpha = %.3f | best rank = %.0f | error = %.5f",
        lambda_laplace,
        results$best_parameter,
        results$best_fit$r,
        results$best_error
      ))
    }
    results$best_fit$lambda_laplace = lambda_laplace
    return(list(fit = results$best_fit, error = results$best_error))
  }
  #----------------
  #' the following goes over a grid of lambda_laplace and then calls what's above.
  #' this one is supposed to run in parallel!! implement that.
  #' output: results$best_parameter: best lambda_laplace, and error. it also returns results$fit
  results <- IMR:::parallel_grid_1d_adaptive(
    param_min = hpar$laplace$lambda_min,
    param_max = hpar$laplace$lambda_max,
    step_sizes = hpar$laplace$lambda_step_sizes,
    f = laplace_cv_lambda_function,
    .progress = TRUE,
    .trace = trace >= 5,
    .packages = c("IMR"),
    .seed = TRUE,
    data = data,
    hpar = hpar,
    row_intercept = row_intercept,
    col_intercept = col_intercept,
    trace = trace,
    thresh = thresh,
    maxit = maxit,
    warm_start = warm_start,
    ls_initial = ls_initial
  )


  if (trace >= 1) {
    message(sprintf(
      "Best lambda_laplace = %.3f | best alpha = %.3f | best rank = %.0f | error = %.5f",
      results$best_fit$lambda_laplace,
      results$best_fit$alpha,
      results$best_fit$r,
      results$best_error
    ))
  }





  # (optional) retrain on the full data
  if (!is.null(data$Y) & final_fit) {
    hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row,
                                                          results$best_fit$lambda_laplace *
                                                            results$best_fit$alpha)
    hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col,
                                                          results$best_fit$lambda_laplace *
                                                            (1 - results$best_fit$alpha))
    fit <- IMR::imr.fit(
      Y = data$Y,
      X = data$Xq,
      Z = data$Zq,
      r = results$best_fit$r,
      lambda_m = 0,
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
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
    results$best_fit[names(fit)] = fit
  }

  # record final hyper-parameters and return
  results$lambda_beta <- hpar$beta$value
  results$lambda_gamma <- hpar$gamma$value
  return(results)
}

#-----------------------------------------------------------
#' @export
imr.cv <- function(
  data,
  row_intercept = FALSE,
  col_intercept = FALSE,
  hpar = IMR::get_imr_default_hparams(),
  error_function = IMR:::error_metric$rmse,
  thresh = 1e-6,
  maxit = 500,
  trace = 0,
  ls_initial = FALSE,
  num_cores = num_cores,
  warm_start = NULL,
  seed = NULL,
  #fast.cv = FALSE,
  separate_tuning = FALSE
) {
  # if (fast.cv) {
  #   return(
  #     IMR:::nlrr.cv(
  #       data = data,
  #       row_intercept = row_intercept,
  #       col_intercept = col_intercept,
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
  stopifnot(is_incomplete(data$Y))
  stopifnot(is_incomplete(data$y_train))
  stopifnot(is_incomplete(data$y_valid))
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
      row_intercept = row_intercept,
      col_intercept = col_intercept,
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
      y_valid = data$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 100,
      verbose = trace
    )
  }
  if (gamma_flag & is.null(hpar$gamma$max)) {
    hpar$gamma$max <- get_lambda_lasso_max(
      y_train = data$y_train,
      Z = data$Zq,
      y_valid = data$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 100,
      verbose = trace
    )
  }


  #---------------------------
  # parallel setup
  # the following function takes
  single_fit <- function(parameter, type="rows", data, row_intercept, col_intercept,
                         hpar, error_function, thresh, trace, maxit, ls_initial,
                         seed, num_cores, fit = NULL){

    if(type == "rows"){
      hpar$beta$value <- parameter
    }else
      hpar$gamma$value <- parameter
    #---------------------------------------
    results <- IMR:::imr.cv_laplace(
      data = data,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      hpar = hpar,
      error_function = error_function,
      thresh = thresh,
      trace = trace - 2,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed,
      num_cores = num_cores,
      warm_start = fit,
      final_fit = FALSE
    )
    if (trace >= 2) {
      message(sprintf(
        paste0( "lambda_beta = %.2f | lambda_gamma = %.2f | ",
          "best lambda_laplace = %.2f | best alpha = %.2f | best rank = %.0f | error = %.3f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$alpha,
        results$best_fit$r,
        results$best_error
      ))
    }
    results$best_fit$lambda_beta = hpar$beta$value
    results$best_fit$lambda_gamma = hpar$gamma$value
    return(list(fit = results$best_fit, error = results$best_error))
  }
  #----------------------------------------------------------------------------
  #if (separate_tuning & gamma_flag & beta_flag) {
    message("Fitting lambda_beta and lambda_gamma separately...")
    # if separate then tune beta first followed by gamma. Meanwhile,
    # keep lambda_gamma at 10% of its maximum.
  if(beta_flag){

    results <- IMR::adaptive_tuner(
      single_fit,
      step_sizes = hpar$beta$step_sizes,
      start_value = hpar$beta$max,
      end_value = 0,
      inc_streak_to_stop = hpar$beta$n_streaks,
      type = "rows",
      data = data,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      hpar = hpar,
      error_function = error_function,
      thresh = thresh,
      trace = trace,
      maxit = maxit,
      ls_initial = ls_initial,
      seed = seed,
      .warm_start = warm_start,
      num_cores = num_cores#,
      #fit = NULL
    )
    hpar$beta$value <- results$best_fit$lambda_beta
    if (trace >= 1) {
      message(sprintf(
        paste0( "best lambda_beta = %.2f | lambda_gamma = %.2f | ",
                "best lambda_laplace = %.2f | best alpha = %.2f | best rank = %.0f | error = %.3f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$alpha,
        results$best_fit$r,
        results$best_error
      ))
    }

  }
  if(gamma_flag){
    # we now do the same to gamma
    results <- IMR::adaptive_tuner(
      single_fit,
      step_sizes = hpar$gamma$step_sizes,
      start_value = hpar$gamma$max,
      end_value = 0,
      inc_streak_to_stop = hpar$gamma$n_streaks,
      type = "cols",
      data = data,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      hpar = hpar,
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
        paste0( "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
                "best lambda_laplace = %.2f | best alpha = %.2f | best rank = %.0f | error = %.3f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$alpha,
        results$best_fit$r,
        results$best_error
      ))
    }

}
    if (!is.null(data$Y)) {
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row,
                                                            results$best_fit$lambda_laplace *
                                                              results$best_fit$alpha)
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col,
                                                            results$best_fit$lambda_laplace *
                                                              (1 - results$best_fit$alpha))
      fit <- IMR::imr.fit(
        Y = data$Y,
        X = data$Xq,
        Z = data$Zq,
        r = results$best_fit$r,
        lambda_m = 0,
        lambda_beta = hpar$beta$value,
        lambda_gamma = hpar$gamma$value,
        row_intercept = row_intercept,
        col_intercept = col_intercept,
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
      results$best_fit[names(fit)] = fit
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
#       row_intercept = row_intercept,
#       col_intercept = col_intercept,
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
#       row_intercept = row_intercept,
#       col_intercept = col_intercept,
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
#       row_intercept = row_intercept,
#       col_intercept = col_intercept,
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


