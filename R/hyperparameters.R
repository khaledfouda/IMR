#----------------------------------------------------------------
#' @export
get_imr_default_hpar <- function(){
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
    # W_valid = NULL,
    y       = NULL,
    # row_cov = TRUE,
    intercept_row = TRUE,
    intercept_col = TRUE,
    hpar = get_imr_default_hpar(),
    interior_loop_length = 20,
    maxit   = 100,
    verbose = 0,
    tol     = 1
) {

  if(!xor(is.null(X), is.null(Z)))
    stop("Either X or Z must be provided, but not both or neither.")
  lambda_beta <- lambda_gamma <- 0
  row_cov <- is.null(Z)
  nr <- nrow(y_train)
  nc <- ncol(y_train)
  # step 1: get an initial fit and find suitable lambda_M and rank before starting:
  if (is.null(y_valid)) {
    mfit <- list()
    # if no validation set provided then fit with a generic r and lambda_M
    mfit$fit <- IMR::imr.fit(
      Y           = y_train,
      X           = X,
      Z           = Z,
      r           = 5,
      lambda_M    = 0,
      lambda_beta = 0,
      lambda_gamma = 0,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit       = maxit,
      trace       = F
    )
    lambda_M <- 0
    r = 5
  } else {
    mfit <- IMR::imr.cv_M(
      y_train     = y_train,
      y_valid     = y_valid,
      X           = X,
      Z           = Z,
      Y_full      = y,
      lambda_beta = 0,
      lambda_gamma = 0,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      trace       = F,
      maxit       = maxit,
      hpar        = hpar
    )
    lambda_M = mfit$lambda_M
    J        = max(2, mfit$rank_M) # do not want the rank to be below 2
  }

  ##  Compute max value using kkt ------------
  residuals  <- y_train -
    mfit$fit$u %*% (mfit$fit$d * t(mfit$fit$v))
  if(intercept_row)
    residuals <- residuals - mfit$fit$beta0 %*% matrix(1, 1, nc)
  if(intercept_col)
    residuals <- residuals - matrix(1, nr, 1) %*% mfit$fit$gamma0
  if(row_cov){
    lambda_max <- max(crossprod(X, residuals))
  }else
    lambda_max <- max(residuals %*% Z)
  ## main_fit function -------------------
  onefit <- function(lambda_beta, lambda_gamma){
    IMR::imr.fit(
      Y             = y_train,
      X             = X,
      Z             = Z,
      r             = r,
      lambda_M      = lambda_M,
      lambda_beta   = lambda_beta,
      lambda_gamma  = lambda_gamma,
      intercept_row = intercept_row,
      intercept_col = intercept_col,
      maxit         = maxit,
      trace         = FALSE
    )
  }

  ## 3) Line search for supremum value ---------------------------
  lambda_sup <- lambda_max
  mid_pt <- lambda_max / 2
  # we fit at mid point
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
      "%s: lambda: max = %.3f; sup = %.3f (%.1f%% of max)",
      ifelse(row_cov, "Beta", "Gamma"),
      lambda_max,
      lambda_sup,
      100 * lambda_sup / lambda_max
    ))
  }

  lambda_sup
}






