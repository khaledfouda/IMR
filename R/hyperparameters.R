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
#' this return that max value for either row covariates or column covariates.
#' you must either privde X or Z. do not provide both
#' @export
get_lambda_lasso_max <- function(
    y_train,
    X       = NULL,
    Z       = NULL,
    y_valid = NULL,
    W_valid = NULL,
    y       = NULL,
    # row_cov = TRUE,
    intercept_row = TRUE,
    intercept_col = TRUE,
    hparams = get_imr_default_hparams(),
    interior_loop_length = 20,
    maxit   = 100,
    verbose = 0,
    tol     = 1
) {

  if( ((!is.null(X)) & (!is.null(Z))) & (is.null(X)&is.null(Z)))
    stop("Either X or Z must be provided. You cannot provide both or neither.")
  lambda_beta <- lambda_gamma <- 0
  # lambda_beta <- if(row_cov) 0 else NULL
  # lambda_gamma <- if(row_cov) NULL else 0
  # if(row_cov){
  #   stopifnot(! is.null(X))
  #   Z = NULL
  # }else{
  #   stopifnot(! is.null(Z))
  #   X = NULL
  # }
  nr <- nrow(y_train)
  nc <- ncol(y_train)
  #-...
  if (is.null(y_valid) || is.null(W_valid)) {
    mfit <- list()
    mfit$fit <- CAMC_fit(
      y           = y_train,
      X           = X,
      Z           = Z,
      J           = 10,
      lambda_M    = 1,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit       = maxit,
      trace       = F
    )
  } else {
    mfit <- CAMC_tv_M(
      y_train     = y_train,
      X           = X,
      Z           = Z,
      y_valid     = y_valid,
      W_valid     = W_valid,
      y           = y,
      lambda_beta = lambda_beta,
      lambda_gamma = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      trace       = F,
      maxit       = maxit,
      hpar        = hparams
    )
  }
  lambda_M = mfit$lambda_M
  J        = max(2, mfit$rank_M) # do not want the rank to be below 2
  ## 2) Compute λ_max  --------------------------------
  residuals  <- y_train -
    mfit$fit$u %*% (mfit$fit$d * t(mfit$fit$v))
  if(intercept_row)
    residuals <- residuals - mfit$fit$phi.a %*% matrix(1, 1, nc)
  if(intercept_col)
    residuals <- residuals - matrix(1, nr, 1) %*% mfit$fit$phi.b
  if(row_cov){
    lambda_max <- max(crossprod(X, residuals))
  }else
    lambda_max <- max(residuals %*% Z)
  ## main_fit function -------------------
  onefit <- function(lambda_beta, lambda_gamma){
    CAMC_fit(
      y             = y_train,
      X             = X,
      Z             = Z,
      J             = J,
      lambda_M      = lambda_M,
      lambda_beta   = lambda_beta,
      lambda_gamma  = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit         = maxit,
      trace         = FALSE
    )
  }

  ## 3) Line search for supremum λ -------------------------------------------
  lambda_sup <- lambda_max
  mid_pt <- lambda_max / 2
  mfit <- onefit(if(row_cov) mid_pt else lambda_beta,
                 if(row_cov) lambda_gamma else mid_pt)
  zero_ratio <- if(row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

  #-- now either search the right side or the left side
  count <- 0L
  if (zero_ratio < 1) {
    # search in [mid_pt, lambda_max]
    old_lambda <- lambda_max
    for (lam in seq(lambda_max, mid_pt, length.out = interior_loop_length)) {

      mfit <- onefit(if(row_cov) lam else lambda_beta,
                     if(row_cov) lambda_gamma else lam)

      zero_ratio <- if(row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

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

      mfit <- onefit(if(row_cov) lam else lambda_beta,
                     if(row_cov) lambda_gamma else lam)

      zero_ratio <- if(row_cov) mean(mfit$beta == 0) else mean(mfit$gamma == 0)

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
      "%s: λ_max = %.3f; λ_sup = %.3f (%.1f%% of λ_max)",
      ifelse(row_cov, "Beta", "Gamma"),
      lambda_max,
      lambda_sup,
      100 * lambda_sup / lambda_max
    ))
  }

  lambda_sup
}






