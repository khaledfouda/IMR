#-----------------------------------------------------------
#' @export
nlrr.cv <- function(
    Y = NULL,
    y_train = NULL,
    y_valid = NULL,
    X = NULL,
    Z = NULL,
    intercept_row = FALSE,
    intercept_col = FALSE,
    lambda_beta = NULL,
    lambda_gamma = NULL,
    val_prop = 0.2,
    hpar = get_imr_default_hparams(),
    error_function = error_metric$rmse,
    thresh = 1e-6,
    maxit = 300,
    verbose = 0,
    max_cores = 8,
    seed = NULL
)
{
  #-------------------
  # check input, perform train/test split if needed, and set seed
  if(is.null(Y) & (is.null(y_train)|is.null(y_valid)))
    stop("You must either provide Y or (y_train,y_valid), or both")
  if(is.null(y_train)| is.null(y_valid)){
    stopifnot(is.Incomplete(Y))
    message("Performing train/valid split")
    obs_mask <- as.matrix(Y != 0)
    valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
    y_train <- as(Y * (1-valid_mask), "Incomplete")
    y_valid <- as(Y * (valid_mask), "Incomplete")
    rm(obs_mask)
    #rm(valid_mask)
  }else{
    stopifnot(is.Incomplete(y_train))
    stopifnot(is.Incomplete(y_valid))
  }
  if((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  #-------------------------------
  # set flags
  beta_flag <- !(is.null(X))
  gamma_flag <- !(is.null(Z))
  # if neither beta or gamma are provided then send to cv_M
  if(! (beta_flag | gamma_flag) )
    stop("Covariates must be provided.")

  # obtain upperbounds to the lambda hyperparameters
  if(beta_flag & is.null(hpar$beta$lambda_max) & is.null(lambda_beta)){
    hpar$beta$lambda_max <- get_lambda_lasso_max(
      y_train = y_train,
      X = X,
      y_valid = y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = verbose
    )
  }
  if(gamma_flag & is.null(hpar$gamma$lambda_max) & is.null(lambda_gamma)){
    hpar$gamma$lambda_max <- get_lambda_lasso_max(
      y_train = y_train,
      Z = Z,
      y_valid = y_valid,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit = 100,
      verbose = verbose
    )
  }

  if(beta_flag & is.null(lambda_beta)){
    lambda_beta_grid <- seq(
      from = hpar$beta$lambda_max,
      to  = 0,
      length.out = hpar$beta$n.lambda
    )
  }else
    lambda_beta_grid <- c(if(is.null(lambda_beta)) 0 else lambda_beta)

  if(gamma_flag & is.null(lambda_gamma)){
    lambda_gamma_grid <- seq(
      from = hpar$gamma$lambda_max,
      to = 0,
      length.out = hpar$gamma$n.lambda
    )
  }else
    lambda_gamma_grid <- c(if(is.null(lambda_gamma)) 0 else lambda_gamma)

  #---------------------------
  # parallel setup
  inner_trace = verbose > 2
  grid <- list(lambda_beta  = lambda_beta_grid,
               lambda_gamma = lambda_gamma_grid)

  nllr.fit.step <- function(lambda_beta, lambda_gamma, ...){
    if(! is.null(seed)) set.seed(seed)
    fit <- IMR::imr.fit_no_low_rank(
    Y = y_train,
    X = X,
    Z = Z,
    lambda_beta = lambda_beta,
    lambda_gamma = lambda_gamma,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    thresh = thresh,
    maxit = maxit,
    trace = FALSE)
    estim = 0
    if(!is.null(fit$beta))
      estim <- partial_crossprod(X, fit$beta, y_valid@i, y_valid@p)
    if(!is.null(fit$gamma))
      estim <- estim + partial_crossprod(fit$gamma, Z, y_valid@i, y_valid@p, TRUE)
    if(!is.null(fit$beta0))
      estim <- estim + partial_crossprod(matrix(fit$beta0,ncol=1),
                                         matrix(1,1,ncol(y_train)),
                                         y_valid@i, y_valid@p)
      # add_to_rows_inplace_cpp(y_valid@x, y_valid@i, old_fit$beta0)
    if(!is.null(fit$gamma0))
      estim <- estim + partial_crossprod(matrix(1,nrow(y_train),1),
                                         matrix(fit$gamma0, nrow=1),
                                         y_valid@i, y_valid@p)
      # add_to_cols_inplace_cpp(y_valid@x, y_valid@p, old_fit$gamma0)

    fit$error <- error_function(estim, y_valid@x)
    fit$lambda_beta <- lambda_beta
    fit$lambda_gamma <- lambda_gamma
    fit
  }

  results <- parallel_grid(grid, nllr.fit.step,
                           "list",
                           .packages = "IMR",
                           .progress = TRUE,
                           .seed    = seed,
                           y_train = y_train,
                           y_valid = y_valid,
                           X = X,
                           Z = Z,
                           intercept_row = intercept_row,
                           intercept_col = intercept_col,
                           thresh = thresh,
                           maxit = maxit,
                           trace = inner_trace,
                           seed = seed,
                           error_function = error_function
  )


  # Select the best fit
  errors <- vapply(results, `[[`, numeric(1), "error")
  best_idx <- which.min(errors)
  best_fit <- results[[best_idx]]
  #--------------------
  # message >>
  if (verbose >= 1) {
    for (res in results) {
      message(sprintf(
        "<< lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f  >>",
        res$lambda_beta,
        sum(res$beta == 0) / length(res$beta),
        res$lambda_gamma,
        sum(res$gamma == 0) / length(res$gamma),
        res$error
      ))
    }
    message(sprintf(
      "<< Best fit >> lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f >>",
      best_fit$lambda_beta,
      sum(best_fit$beta == 0) / length(best_fit$beta),
      best_fit$lambda_gamma,
      sum(best_fit$gamma == 0) / length(best_fit$gamma),
      best_fit$error

    ))
    best_fit$init_hparams <- hpar
  }
  rm(results)

 nfit <- IMR::imr.fit_no_low_rank(
    Y = Y,
    X = X,
    Z = Z,
    lambda_beta  = best_fit$lambda_beta,
    lambda_gamma = best_fit$lambda_gamma,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    thresh = thresh,
    maxit = maxit,
    trace = FALSE)


  y_train <- as(nfit$resid * (1-valid_mask), "Incomplete")
  y_valid <- as(nfit$resid * (valid_mask), "Incomplete")
  mfit <- IMR::imr.cv_M(y_train, y_valid, Y_full = nfit$resid, hpar=hpar,
                        error_function = error_function)

  nfit$u <- mfit$fit$u
  nfit$v <- mfit$fit$v
  nfit$d <- mfit$fit$d

  best_fit$fit <- IMR::imr.fit(
    Y = Y,
    X = X,
    Z = Z,
    r = mfit$rank_M,
    lambda_M = mfit$lambda_M,
    lambda_beta = best_fit$lambda_beta,
    lambda_gamma = best_fit$lambda_gamma,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    thresh = thresh,
    maxit = maxit,
    trace = inner_trace,
    warm_start = nfit
    )
  mfit$rank_M -> best_fit$rank_M
  mfit$lambda_M -> best_fit$lambda_M

  return(best_fit)
}

