#' @export
imr.cv_M <- function(
  y_train,
  y_valid,
  X = NULL,
  Z = NULL,
  Y_full = NULL,
  lambda_beta = 0,
  lambda_gamma = 0,
  intercept_row = FALSE,
  intercept_col = FALSE,
  hpar = get_imr_default_hparams(),
  error_function = error_metric$rmse,
  thresh = 1e-6,
  maxit = 300,
  trace = TRUE,
  old_fit = NULL,
  seed = NULL
){
  # set seed and check input matrix type
  if( (!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  stopifnot(is.Incomplete(y_train))
  stopifnot(is.Incomplete(y_valid))
  if(!is.null(Y_full)) stopifnot(is.Incomplete(Y_full))

  # lambda lambda_m sequence
  if(is.null(hpar$M$lambda_max)){
    hpar$M$lambda_max <- get_lambda_M_max(y_train, X, Z, T, T,
                                           lambda_beta, lambda_gamma) *
      hpar$M$lambda_factor
  }
  lambda_seq <- seq(
    from = hpar$M$lambda_max,
    to  = 0,
    length.out = hpar$M$n.lambda
  )

  # extract indices from the mask and create flags
  virow <- y_valid@i
  vpcol <- y_valid@p
  reference <- y_valid@x[]
  beta_flag <- !(is.null(lambda_beta) | is.null(X))
  gamma_flag <- !(is.null(lambda_gamma) | is.null(Z))

  # initialize vars
  rank_max <- hpar$M$rank.init
  best_fit <- list(error = Inf)
  no_improve_count <- 0
  loop_size <- 0

  for(i in seq_along(lambda_seq)){
    loop_size <- loop_size + 1
    # fit

    old_fit <- IMR::imr.fit(
      Y = y_train,
      X = X,
      Z = Z,
      r = rank_max,
      lambda_M = lambda_seq[i],
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      lambda_a = hpar$laplacian$lambda_a,
      lambda_b = hpar$laplacian$lambda_b,
      L_a = hpar$laplacian$L_a,
      L_b = hpar$laplacian$L_b,
      warm_start = old_fit,
      trace = F,
      thresh = thresh,
      maxit = maxit
    )

    # compute validation error
    y_valid@x <- partial_crossprod(old_fit$u, old_fit$d * t(old_fit$v),virow,vpcol)
    if(beta_flag)
      y_valid@x <- y_valid@x + partial_crossprod(X, old_fit$beta, virow, vpcol)
    if(gamma_flag)
      y_valid@x <- y_valid@x + partial_crossprod(old_fit$gamma, Z, virow, vpcol, TRUE)
    if(intercept_row)
      add_to_rows_inplace_cpp(y_valid@x, y_valid@i, old_fit$beta0)
    if(intercept_col)
      add_to_cols_inplace_cpp(y_valid@x, y_valid@p, old_fit$gamma0)

    verror <- error_function(y_valid@x, reference)

    # compute rank
    current_rank <- sum(round(old_fit$d, 4) > 0)

    # verbose
    if (trace) {
      message(sprintf(
        "%2d lambda=%.4g | rank_max=%d => rank=%d | err=%.5f | iters=%d",
        i,
        lambda_seq[i],
        rank_max,
        current_rank,
        verror,
        old_fit$n_iter
      ))
    }

    # track best model & early stopping
    if (verror < best_fit$error) {
      best_fit <- list(
        error     = verror,
        rank_M    = current_rank,
        lambda_M  = lambda_seq[i],
        rank_max  = rank_max,
        fit       = old_fit
      )
      no_improve_count <- 0
    } else {
      no_improve_count <- no_improve_count + 1
    }

    if (no_improve_count >= hpar$M$early.stopping) {
      if (trace) {
        message(
          sprintf(
            "Early stopping: no improvement in last %d lambda’s.",
            no_improve_count
          )
        )
      }
      break
    }

    # update rank_max (r) for next iteration
    rank_max <- min(
      current_rank + hpar$M$rank.step,
      hpar$M$rank.max
    )
    rank_max <- max(
      rank_max,
      hpar$M$rank.min
    )
    # end of loop
  }
  # (optional) retrain on the full data
  if(!is.null(Y_full)){
    best_fit$fit <- imr.fit(
      Y = Y_full,
      X = X,
      Z = Z,
      r = best_fit$rank_max,
      lambda_M = best_fit$lambda_M,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      lambda_a = hpar$laplacian$lambda_a,
      lambda_b = hpar$laplacian$lambda_b,
      L_a = hpar$laplacian$L_a,
      L_b = hpar$laplacian$L_b,
      warm_start = old_fit,
      trace = FALSE,
      thresh = thresh,
      maxit = maxit
    )
    loop_size <- loop_size + 1
  }

  # record final hyper-parameters and return
  best_fit$lambda_beta <- lambda_beta
  best_fit$lambda_gamma <- lambda_gamma
  best_fit$loop_size  <- loop_size
  return(best_fit)
}

#-----------------------------------------------------------
#' @export
imr.cv <- function(
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
    rm(valid_mask)
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
    return(IMR::imr.cv_M(
      y_train = y_train,
      y_valid = y_valid,
      Y_full = Y,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      hpar = hpar,
      error_function = error_function,
      thresh = thresh,
      trace = verbose > 0,
      maxit = maxit,
      seed = seed
    ))

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

  results <- parallel_grid(grid, IMR::imr.cv_M,
                           "list",
                           .packages = "IMR",
                           .progress = TRUE,
                           .seed    = seed,
                           y_train = y_train,
                           y_valid = y_valid,
                           X = X,
                           Z = Z,
                           Y_full = Y,
                           intercept_row = intercept_row,
                           intercept_col = intercept_col,
                           hpar = hpar,
                           error_function = error_function,
                           thresh = thresh,
                           maxit = maxit,
                           trace = inner_trace,
                           seed = seed
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
        "<< lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
        res$lambda_beta,
        sum(res$fit$beta == 0) / length(res$fit$beta),
        res$lambda_gamma,
        sum(res$fit$gamma == 0) / length(res$fit$gamma),
        res$error,
        res$loop_size,
        res$rank_M,
        res$lambda_M
      ))
    }
    message(sprintf(
      "<< Best fit >> lambda_beta=%.4g | sparsity=%.2f | lambda_gamma=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
      best_fit$lambda_beta,
      sum(best_fit$fit$beta == 0) / length(best_fit$fit$beta),
      best_fit$lambda_gamma,
      sum(best_fit$fit$gamma == 0) / length(best_fit$fit$gamma),
      best_fit$error,
      best_fit$loop_size,
      best_fit$rank_M,
      best_fit$lambda_M
    ))
    best_fit$init_hparams <- hpar
  }
  rm(results)
  return(best_fit)
}





