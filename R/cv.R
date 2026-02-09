#----------------------------------------------------------
#' @export
imr.cv_laplace <- function(
  data,
  intercept_row = FALSE,
  intercept_col = FALSE,
  hpar = IMR::get_imr_default_hparams(),
  error_function = IMR:::error_metric$rmse,
  thresh = 1e-4,
  maxit = 300,
  trace = 1,
  warm_start = NULL,
  ls_initial = TRUE,
  shared_information = FALSE,
  num_cores = 6,
  final_fit = TRUE,
  final_thresh = 1e-6,
  final_maxit = 1000,
  seed = NULL
) {
  stopifnot(is.Incomplete(data$Y))
  stopifnot(is.Incomplete(data$y_train))
  stopifnot(is.Incomplete(data$y_valid))
  if ((!is.null(seed)) && is.numeric(seed)) set.seed(seed)
  if( is.numeric(num_cores) && num_cores > 0) IMR::initialize_parallel_workers(num_cores)
  #---------------------------------------------------

  # fixed: all. variable: none. number of fits: 1.
  # fits a single fit and returns [fit, error]; with all lambdas fixed.
  # this is a single fit where all parameters are fixed but it returns validation error
  rank_fit_function <- function(r, fdata, hpar, shared_information,
                                lambda_laplace,
                                intercept_row, intercept_col,
                                trace, thresh, maxit,
                                ls_initial, fit=NULL,
                                error_function = IMR:::error_metric$rmse) {

    fit <- IMR::imr.fit(
      Y = fdata$y_train,
      X = fdata$Xq,
      Z = fdata$Zq,
      r = r,
      lambda_M = lambda_laplace,
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      shared_information = shared_information,
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


    vestim <- IMR:::reconstruct_partial(fit, fdata, fdata$y_valid, shared_information)
    verror <- error_function(fdata$y_valid@x, vestim@x)
    # verbose

    if (trace >= 3) message(sprintf("rank=%d | err=%.5f ", r, verror))
    fit$r = r
    return(list(fit = fit, error = verror))
  }

  #---
  # this function takes a single lambda_laplace and finds the optimal rank. everything else fixed.
  laplace_cv_lambda_function <- function(lambda_laplace, data, hpar,
                                        intercept_row, shared_information,
                                        intercept_col, trace, thresh,
                                        maxit, ls_initial, fit = NULL, warm_start = NULL) {
    if(!is.null(data$similarity_rows) && lambda_laplace > 0)
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row, lambda_laplace)
    if(!is.null(data$similarity_cols) && lambda_laplace > 0)
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col, lambda_laplace)

    results <- IMR::adaptive_tuner(rank_fit_function,
      step_sizes = hpar$rank$step_sizes,
      start_value = hpar$rank$min,
      end_value = hpar$rank$max,
      inc_streak_to_stop = hpar$rank$n_streaks,
      fdata = data,
      hpar = hpar,
      lambda_laplace = lambda_laplace,
      shared_information = shared_information,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      trace = trace,
      thresh = thresh,
      .warm_start = warm_start,
      maxit = maxit,
      ls_initial = ls_initial
    )
    if (trace >= 2) {
      message(sprintf(
        " lambda_laplace = %.3f | best rank = %.0f | error = %.5f",
        lambda_laplace,
        results$best_parameter,
        results$best_error
      ))
    }
    results$best_fit$lambda_laplace = lambda_laplace
    return(list(fit = results$best_fit, error = results$best_error, history = results$history))
  }

  #----------------
  # the following goes over a grid of lambda_laplace and then calls what's above.
  # this one is supposed to run in parallel!! implement that.
  # output: results$best_parameter: best lambda_laplace, and error. it also returns results$fit
  results <- IMR:::parallel_grid_1d_adaptive(
    param_min = hpar$laplace$min,
    param_max = hpar$laplace$max,
    step_sizes = hpar$laplace$step_sizes,
    f = laplace_cv_lambda_function,
    .progress = TRUE,
    .trace = trace >= 4,
    #.packages = c("IMR"),
    .seed = if(is.null(seed)) FALSE else seed,
    data = data,
    hpar = hpar,
    shared_information = shared_information,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    trace = trace,
    thresh = thresh,
    maxit = maxit,
    warm_start = warm_start,
    ls_initial = ls_initial
  )


  if (trace >= 1) {
    message(sprintf(
      "Best lambda_laplace = %.3f | best rank = %.0f | error = %.5f",
      results$best_fit$lambda_laplace,
      results$best_fit$r,
      results$best_error
    ))
  }





  # (optional) retrain on the full data
  if (!is.null(data$Y) & final_fit) {
    if(!is.null(data$similarity_rows) && results$best_fit$lambda_laplace > 0)
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row,
                                                          results$best_fit$lambda_laplace)
    if(!is.null(data$similarity_cols) && results$best_fit$lambda_laplace > 0)
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col,
                                                          results$best_fit$lambda_laplace)
    fit <- IMR::imr.fit(
      Y = data$Y,
      X = data$Xq,
      Z = data$Zq,
      r = results$best_fit$r,
      lambda_M = results$best_fit$lambda_laplace,
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
      thresh = final_thresh,
      maxit = final_maxit,
      ls_initial = ls_initial
    )
    results$best_fit[names(fit)] = fit
  }

  # record final hyper-parameters and return
  results$lambda_beta <- hpar$beta$value
  results$lambda_gamma <- hpar$gamma$value
  results$fit <- results$best_fit
  results$fit$params <-
    list(
      lambda_beta = results$lambda_beta,
      lambda_gamma = results$lambda_gamma,
      lambda_laplace = results$fit$lambda_laplace,
      rank = results$fit$r
    )
  return(results)
}

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
  #fast.cv = FALSE,
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
      #y_valid = data$y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = trace
    )
    if (is.function(hpar$beta$step_sizes))
      hpar$beta$step_sizes <- hpar$beta$step_sizes(hpar$beta$min, hpar$beta$max)
  }
  if (gamma_flag & is.null(hpar$gamma$max)) {
    hpar$gamma$max <- IMR::get_lambda_lasso_max(
      y_train = data$y_train,
      Z = data$Zq,
      #y_valid = data$y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = trace
    )
    if(is.function(hpar$gamma$step_sizes))
      hpar$gamma$step_sizes <- hpar$gamma$step_sizes(hpar$gamma$min, hpar$gamma$max)
  }


  #---------------------------
  # parallel setup
  # the following function takes
  single_fit <- function(parameter, type="rows", data, intercept_row, intercept_col,
                         shared_information,
                         hpar, error_function, thresh, trace, maxit, ls_initial,
                         seed, num_cores, fit = NULL){

    if(type == "rows"){
      hpar$beta$value <- parameter
    }else
      hpar$gamma$value <- parameter
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
      #warm_start = fit,
      final_fit = FALSE
    )
    if (trace >= 2) {
      message(sprintf(
        paste0( "lambda_beta = %.5f | lambda_gamma = %.5f | ",
          "best lambda_laplace = %.3f | best rank = %.0f | error = %.6f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
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
    if(layer_1_parallel){

    results <- IMR:::parallel_grid_1d_adaptive(
      param_min = hpar$beta$min,
      param_max = hpar$beta$max,
      step_sizes = hpar$beta$step_sizes,
      f = single_fit,
      .progress = TRUE,
      .trace = trace >= 3,
      #.packages = c("IMR"),
      .seed = if(is.null(seed)) FALSE else seed,
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
      #.warm_start = warm_start
    )

    }else{
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
      num_cores = num_cores#,
      #fit = NULL
    )
    }
    hpar$beta$value <- results$best_fit$lambda_beta
    if (trace >= 1) {
      message(sprintf(
        paste0( "best lambda_beta = %.5f | lambda_gamma = %.5f | ",
                "best lambda_laplace = %.4f | best rank = %.0f | error = %.5f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
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
        paste0( "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
                "best lambda_laplace = %.2f |  best rank = %.0f | error = %.5f"),
        hpar$beta$value,
        hpar$gamma$value,
        results$best_fit$lambda_laplace,
        results$best_fit$r,
        results$best_error
      ))
    }

}
    if (!is.null(data$Y)) {
      if(!is.null(data$similarity_rows) && results$best_fit$lambda_laplace > 0)
        hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row,
                                                              results$best_fit$lambda_laplace)
      if(!is.null(data$similarity_cols) && results$best_fit$lambda_laplace > 0)
        hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col,
                                                              results$best_fit$lambda_laplace)

      fit <- IMR::imr.fit(
        Y = data$Y,
        X = data$Xq,
        Z = data$Zq,
        r = results$best_fit$r,
        lambda_M = results$best_fit$lambda_laplace,
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
  #       res$lambda_M
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
  #     best_fit$lambda_M
  #   ))
  #   best_fit$init_hparams <- hpar
  # }
  # rm(results)
  # return(best_fit)


