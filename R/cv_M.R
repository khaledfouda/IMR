#' @export
imr.cv.M <- function(
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
    hpar$M$lambda_max <- lambdaM.max(y_train, X, Z, T, T,
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
    old_fit <- imr.fit(
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
      trace = FALSE,
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

#----------------------------------------------------------------
#' @export
get_imr_default_hparams <- function(){
  list(
    M = list(
      lambda_max    = NULL,   # = lambda_{M,min}
      lambda_factor  = 1 / 4,  # ignored if lambda_init is provided
      n.lambda       = 20,     # sequence from lambda_init to 0 (inclusive)
      rank.init      = 2,
      rank.max       = 30,
      rank.min       = 2,
      rank.step      = 2,
      early.stopping = 1
    ),
    beta = list(
      lambda_max = NULL,  # if NULL, computed internally (recommended)
      n.lambda   = 20,
      init.tol   = 3
    ),
    gamma = list(
      lambda_max = NULL,
      n.lambda   = 20,
      init.tol   = 3
    ),
    laplacian = list(
      lambda_a = 0,
      L_a      = NULL,
      lambda_b = 0,
      S_b      = NULL
    )
  )
}


#-----------------------------------------------------------
#' @export
lambdaM.max <-
  function(Y,
           X = NULL,
           Z = NULL,
           intercept_row = FALSE,
           intercept_col = FALSE,
           lambda_beta = NULL,
           lambda_gamma = NULL,
           maxit = 30){
    need_fit <- any(!is.null(X), !is.null(Z), intercept_row, intercept_col)


    if (need_fit) {
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
    }
    # return largest singular value
    svd::propack.svd(naive_MC(as.matrix(mfit$resid)),neig =  1, opts = list(kmax = maxit))$d[1]
  }
#------------------------------------------------

