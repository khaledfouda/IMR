#' @importFrom foreach foreach %dopar% %:%
#' @export
imr.cv_21 <- function(
  data,
  row_intercept = FALSE,
  col_intercept = FALSE,
  hpar = IMR::get_imr_default_hparams(),
  error_function = IMR:::error_metric$rmse,
  thresh = 1e-3,
  maxit = 300,
  trace = 0,
  ls_initial = TRUE,
  shared_information = FALSE,
  num_cores = 4,
  warm_start = NULL,
  final_thresh = 1e-6,
  final_maxit = 1000,
  init_thresh = 1e-4,
  init_maxit = 500,
  seed = NULL,

  separate = FALSE
) {
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
      shared_information = shared_information,
      warm_start = warm_start,
      final_fit = TRUE,
      final_thresh = final_thresh,
      final_maxit = final_maxit,
      num_cores = num_cores
    ))
  }

  # obtain upperbounds to the lambda hyperparameters
  if (beta_flag & is.null(hpar$beta$max)) {
    hpar$beta$max <- IMR::get_lambda_lasso_max(
      y_train = data$y_train,
      X = data$Xq,
      # y_valid = data$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 50,
      thresh = 1e-3,
      init_maxit = 100,
      init_thresh = 1e-4,
      r = 5,
      verbose = trace
    )
  }
  if (gamma_flag & is.null(hpar$gamma$max)) {
    hpar$gamma$max <- IMR::get_lambda_lasso_max(
      y_train = data$y_train,
      Z = data$Zq,
      # y_valid = data$y_valid,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      maxit = 50,
      verbose = trace,
      thresh = 1e-3,
      init_maxit = 100,
      init_thresh = 1e-4,
      r = 5,
      tol = thresh
    )
  }
  if(!gamma_flag) hpar$gamma$max <- hpar$gamma$length <- 1
  if(!beta_flag) hpar$beta$max   <- hpar$beta$length <- 1
  # =================================================================================
  #---------------------------------------------------
  # fixed: all. variable: none. number of fits: 1.
  # fits a single fit and returns [fit, error]; with all lambdas fixed.
  # this is a single fit where all parameters are fixed but it returns validation error
  rank_fit_function <- function(r, fdata, hpar, shared_information,
                                lambda_laplace,
                                row_intercept, col_intercept,
                                trace, thresh, maxit,
                                ls_initial, fit = NULL,
                                error_function = IMR:::error_metric$rmse) {
    fit <- IMR::imr.fit(
      Y = fdata$y_train,
      X = fdata$Xq,
      Z = fdata$Zq,
      r = r,
      lambda_m = lambda_laplace,
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
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
    fit$r <- r
    return(list(fit = fit, error = verror))
  }
  # ==================================================================================
  # parallel setup
  # the following function takes
  single_fit <- function(lambda_laplace, lambda_beta, lambda_gamma,
                         data, row_intercept, col_intercept,
                         shared_information,
                         hpar, error_function, thresh, trace, maxit, ls_initial,
                         seed, fit = NULL, warm_start = NULL, return_fit = FALSE) {
    hpar$beta$value <- lambda_beta
    hpar$gamma$value <- lambda_gamma
    if (!is.null(data$similarity_rows) && lambda_laplace > 0) {
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row, lambda_laplace)
    }
    if (!is.null(data$similarity_cols) && lambda_laplace > 0) {
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col, lambda_laplace)
    }

    #---------------------------------------
    results <- IMR::adaptive_tuner(rank_fit_function,
      step_sizes = hpar$rank$step_sizes,
      start_value = hpar$rank$min,
      end_value = hpar$rank$max,
      inc_streak_to_stop = hpar$rank$n_streaks,
      fdata = data,
      hpar = hpar,
      lambda_laplace = lambda_laplace,
      shared_information = shared_information,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      trace = trace,
      thresh = thresh,
      .warm_start = warm_start,
      maxit = maxit,
      ls_initial = ls_initial
    )
    if (trace >= 2) {
      message(sprintf(
        paste0(
          "lambda_beta = %.5f | lambda_gamma = %.5f | ",
          "lambda_laplace = %.3f | best rank = %.0f | error = %.6f"
        ),
        hpar$beta$value,
        hpar$gamma$value,
        lambda_laplace,
        results$best_parameter,
        results$best_error
      ))
    }

    data.frame(
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      lambda_laplace = lambda_laplace,
      rank = results$best_parameter,
      error = results$best_error
    ) -> res
    if(return_fit)
      return(list(res = res, fit = results$best_fit))
    return(res)
  }

  #================================================================================
  initial_fit <- IMR::imr.fit_no_low_rank(data$y_train, data$Xq, data$Zq,
                                          hpar$beta$value,
                                          hpar$gamma$value,
                                          row_intercept = row_intercept,
                                          col_intercept = col_intercept,
                                          shared_information = shared_information,
                                          maxit = init_maxit,
                                          thresh = init_thresh)
  init <- IMR::svd_opt(initial_fit$resid, hpar$rank$default, rthin =  FALSE, cthin =  FALSE)
  initial_fit$u <- init$u
  initial_fit$d <- init$d
  initial_fit$v <- init$v
  initial_fit$resid = NULL
  #================================================================================
  # what we have >> beta$min -> beta$max with length
  #                 gamma$min -> gamma$max with length
  #                 laplace$min -> laplace$max  with step_sizes[1]

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  #parallel::clusterExport(cl, varlist = c("imr.fit"))

  beta_seq = seq(hpar$beta$min, hpar$beta$max, length.out = hpar$beta$length)
  gamma_seq = seq(hpar$gamma$min, hpar$gamma$max, length.out = hpar$gamma$length)

  if(separate){
    gamma_seq_old = gamma_seq
    gamma_seq = c(hpar$gamma$value)
  }
  results <-
    foreach::foreach(
      lambda_beta = beta_seq,
      .combine = rbind
    ) %:%
    foreach::foreach(
      lambda_gamma = gamma_seq,
      .combine = rbind
    ) %dopar% {
      lamseq = seq(hpar$laplace$max, hpar$laplace$min, -hpar$laplace$step_sizes[1])
      warm_start = initial_fit
      results = data.frame()
      for(lambda_laplace in lamseq){

      res <- single_fit(
        lambda_laplace = lambda_laplace,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        data = data,
        row_intercept = row_intercept,
        col_intercept = col_intercept,
        hpar = hpar,
        shared_information = shared_information,
        error_function = error_function,
        thresh = thresh,
        trace = trace,
        maxit = maxit,
        ls_initial = ls_initial,
        seed = seed,
        return_fit = TRUE,
        warm_start = warm_start
      )
      warm_start = res$fit
      results <- rbind(results, res$res)
      }
      results
    }


  # ============================================================================
  dplyr::arrange(results, error, desc(lambda_laplace), desc(lambda_beta), desc(lambda_gamma)) ->
  results
  hpar$beta$value <- results[1, ]$lambda_beta
  hpar$gamma$value <- results[1, ]$lambda_gamma
  if (trace >= 1) {
    message(sprintf(
      paste0( "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
              "best lambda_laplace = %.2f |  best rank = %.0f | error = %.5f"),
      hpar$beta$value,
      hpar$gamma$value,
      results[1,]$lambda_laplace,
      results[1,]$rank,
      results[1,]$error
    ))
  }
  # ===========================================================================
  if(separate){
    gamma_seq = gamma_seq_old
    lambda_beta = hpar$beta$value
    results <-
      foreach::foreach(
        lambda_gamma = gamma_seq,
        .combine = rbind
      ) %dopar% {
        lamseq = seq(hpar$laplace$min, hpar$laplace$max, hpar$laplace$step_sizes[1])
        warm_start = initial_fit
        results = data.frame()
        for(lambda_laplace in lamseq){

          res <- single_fit(
            lambda_laplace = lambda_laplace,
            lambda_beta = lambda_beta,
            lambda_gamma = lambda_gamma,
            data = data,
            row_intercept = row_intercept,
            col_intercept = col_intercept,
            hpar = hpar,
            shared_information = shared_information,
            error_function = error_function,
            thresh = thresh,
            trace = trace,
            maxit = maxit,
            ls_initial = ls_initial,
            seed = seed,
            return_fit = TRUE,
            warm_start = warm_start
          )
          warm_start = res$fit
          results <- rbind(results, res$res)
        }
        results
      }
    # ============================================================================
    dplyr::arrange(results, error, desc(lambda_laplace), desc(lambda_beta), desc(lambda_gamma)) ->
      results
    hpar$beta$value <- results[1, ]$lambda_beta
    hpar$gamma$value <- results[1, ]$lambda_gamma
    # ===========================================================================
    if (trace >= 1) {
      message(sprintf(
        paste0( "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
                "best lambda_laplace = %.2f |  best rank = %.0f | error = %.5f"),
        hpar$beta$value,
        hpar$gamma$value,
        results[1,]$lambda_laplace,
        results[1,]$rank,
        results[1,]$error
      ))
    }
  }


  # ========================== if there's a second loop for lambda_m go for it
  nsteps <- length(hpar$laplace$step_sizes)
  if (nsteps > 1) {
    steps <- hpar$laplace$step_sizes # [1:nsteps_laplace]

    for (i in 2:nsteps) {

      startp <- max(hpar$laplace$min, results[1, ]$lambda_laplace - steps[i-1])
      endp <- min(hpar$laplace$max, results[1, ]$lambda_laplace + steps[i-1])
      grid <- seq(startp, endp, steps[i])
      # remove those that already exist
      grid <- setdiff(grid, unique(results$lambda_laplace))

      results <- rbind(
        results,
        foreach::foreach(lambda_laplace = grid, .combine = rbind) %dopar% {
          single_fit(
            lambda_laplace = lambda_laplace,
            lambda_beta = hpar$beta$value,
            lambda_gamma = hpar$gamma$value,
            data = data,
            row_intercept = row_intercept,
            col_intercept = col_intercept,
            hpar = hpar,
            shared_information = shared_information,
            error_function = error_function,
            thresh = thresh,
            trace = trace,
            maxit = maxit,
            ls_initial = ls_initial,
            seed = seed,
            warm_start = initial_fit
          )
        }
      )
      dplyr::arrange(results, error, desc(lambda_laplace), desc(lambda_beta), desc(lambda_gamma)) ->
        results
    }
  }
  parallel::stopCluster(cl)
  #========= tuning is over, we now do one final fit and return >>
  if (trace >= 1) {
    message(sprintf(
      paste0( "best lambda_beta = %.2f | best lambda_gamma = %.2f | ",
              "best lambda_laplace = %.2f |  best rank = %.0f | error = %.5f"),
      hpar$beta$value,
      hpar$gamma$value,
      results[1,]$lambda_laplace,
      results[1,]$rank,
      results[1,]$error
    ))
  }


  if (!is.null(data$Y)) {
    if(!is.null(data$similarity_rows) && results[1,]$lambda_laplace > 0)
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(data$similarity_row,
                                                            results[1,]$lambda_laplace)
    if(!is.null(data$similarity_cols) && results[1,]$lambda_laplace > 0)
      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(data$similarity_col,
                                                            results[1,]$lambda_laplace)

    fit <- IMR::imr.fit(
      Y = data$Y,
      X = data$Xq,
      Z = data$Zq,
      r = results[1,]$rank,
      lambda_m = results[1,]$lambda_laplace,
      lambda_beta = hpar$beta$value,
      lambda_gamma = hpar$gamma$value,
      row_intercept = row_intercept,
      col_intercept = col_intercept,
      shared_information = shared_information,
      Ur = hpar$laplacian_row$U,
      dr = hpar$laplacian_row$d,
      Uc = hpar$laplacian_col$U,
      dc = hpar$laplacian_col$d,
      warm_start = initial_fit,
      trace = trace >= 3,
      thresh = final_thresh,
      maxit = final_maxit,
      ls_initial = ls_initial
    )


    fit$params <- results[1,]
    return(list(history = results, fit=fit, hpar=hpar))

}
}
