#' Fit Incomplete Matrix Regression (IMR) Model
#'
#' \code{imr.fit} fits the model to the given data and hyper-parameters until
#' convergence is achieved.
#'
#' @param Y An incomplete matrix (class \code{Incomplete}; see \code{\link{as.incomplete}}).
#'   The target matrix to be completed (n by m).
#' @param X Optional matrix of row covariates (n by p). Default is \code{NULL}.
#' @param Z Optional matrix of column covariates (m by q). Default is \code{NULL}.
#' @param intercept_row Logical. Include row-level intercepts? Default is \code{FALSE}.
#' @param intercept_col Logical. Include column-level intercepts? Default is \code{FALSE}.
#' @param r Integer. The rank (number of latent factors/columns in A and B).
#'   Default is 2.
#' @param lambda_M Numeric scalar. Controls the nuclear penalty. Default is 0.
#' @param lambda_beta Numeric scalar. Controls the Lasso penalty on the row
#'   covariates. Default is 0.
#' @param lambda_gamma Numeric scalar. Controls the Lasso penalty on the column
#'   covariates. Default is 0.
#' @param Ur,dr Optional matrix (Ur) and vector (dr) containing the eigenvectors
#'   and eigenvalues of the row-level similarity matrix \eqn{S_r}.
#' @param Vc,dc Optional matrix (Vc) and vector (dc) containing the eigenvectors
#'   and eigenvalues of the column-level similarity matrix \eqn{S_c}.
#' @param maxit Integer. Maximum number of iterations. Default is 300.
#' @param thresh Numeric scalar. Convergence threshold based on the Frobenius
#'   difference between updates. Default is 1e-5.
#' @param trace Logical. If \code{TRUE}, prints the objective function value at
#'   each step. Default is \code{FALSE}.
#' @param warm_start Optional list. A previous result object from \code{imr.fit}
#'   to use as initial values. Default is \code{NULL}.
#' @param ls_initial Logical. Used only if \code{warm_start} is \code{NULL}.
#'   If \code{TRUE} (default), uses least-squares initialization. If \code{FALSE},
#'   uses random initialization.
#'
#' @return A list containing the learned parameters:
#' \describe{
#'   \item{u, d, v}{The SVD decomposition components of the matrix, such that
#'     \deqn{M = u \cdot \textrm{diag}(d) \cdot v^T}}
#'   \item{beta}{Matrix of row covariate coefficients. \code{NULL} if \code{X} is \code{NULL}.}
#'   \item{gamma}{Matrix of column covariate coefficients. \code{NULL} if \code{Z} is \code{NULL}.}
#'   \item{beta0}{Vector of row-level intercepts. \code{NULL} if \code{intercept_row} is \code{FALSE}.}
#'   \item{gamma0}{Vector of column-level intercepts. \code{NULL} if \code{intercept_col} is \code{FALSE}.}
#'   \item{n_iter}{Integer. The number of iterations performed.}
#' }
#'
#' @export
imr.fit <- function(
    Y,
    X = NULL,
    Z = NULL,
    intercept_row = FALSE,
    intercept_col = FALSE,
    r = 2,
    lambda_M = 0,
    lambda_beta = 0,
    lambda_gamma = 0,
    Ur = NULL,
    dr = NULL,
    Uc = NULL,
    dc = NULL,
    maxit = 300,
    thresh = 1e-5,
    trace = FALSE,
    warm_start = NULL,
    ls_initial = TRUE) {
  # Input checks & setup ----------------------------------------------------
  stopifnot(is.Incomplete(Y))

  dims <- dim(Y)
  nr <- dims[1]
  nc <- dims[2]
  nz <- Matrix::nnzero(Y, na.counted = TRUE)

  irow <- Y@i
  pcol <- Y@p


  # Laplacian flags (L_* expected as eigendecompositions) -------------------
  beta_flag <- !(is.null(X))
  gamma_flag <- !(is.null(Z))
  laplace_r_flag <- !(is.null(Ur) | is.null(dr))
  laplace_c_flag <- !(is.null(Uc) | is.null(dc))
  # initial everything to null ------------------------
  beta <- gamma <- beta0 <- gamma0 <- NULL

  # 3) Warm-start or initialize ------------------------------------------------
  warm_start <- verify_warm_start(warm_start, r)

  if (!is.null(warm_start)) {
    required <- c("u", "d", "v", "beta", "gamma", "beta0", "gamma0")
    if (!all(required %in% names(warm_start))) {
      stop("warm_start missing components: u, d, v, beta, gamma, beta0, gamma0")
    }

    if (beta_flag) {
      beta <- warm_start$beta
      xb_obs <- partial_crossprod(X, beta, irow, pcol)
    }
    if (gamma_flag) {
      gamma <- warm_start$gamma
      zg_obs <- partial_crossprod(gamma, Z, irow, pcol, TRUE)
    }
    if (intercept_row) {
      beta0 <- warm_start$beta0
    }

    if (intercept_col) {
      gamma0 <- warm_start$gamma0
    }

    U <- warm_start$u
    V <- warm_start$v
    Dsq <- warm_start$d
  } else {
    if (ls_initial) {
      mfit <- IMR::imr.fit_no_low_rank(Y, X, Z,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        intercept_row = intercept_row,
        intercept_col = intercept_col
      )
      if (beta_flag) {
        beta <- mfit$beta
        xb_obs <- partial_crossprod(X, beta, irow, pcol)
      }
      if (gamma_flag) {
        gamma <- mfit$gamma
        zg_obs <- partial_crossprod(gamma, Z, irow, pcol, TRUE)
      }
      if (intercept_row) {
        beta0 <- mfit$beta0
      }
      if (intercept_col) {
        gamma0 <- mfit$gamma0
      }

      init <- opt_svd(naive_MC(mfit$resid), r, nr, nc, FALSE, FALSE)
    } else {
      if (beta_flag) {
        beta <- matrix(0, ncol(X), nc)
        xb_obs <- rep(0, nz)
      }
      if (gamma_flag) {
        gamma <- matrix(0, nr, ncol(Z))
        zg_obs <- rep(0, nz)
      }
      if (intercept_row) {
        beta0 <- rep(0, nr)
      }
      if (intercept_col) {
        gamma0 <- rep(0, nc)
      }

      init <- opt_svd(naive_MC(Y), r, nr, nc, FALSE, FALSE)
    }

    U <- init$u
    Dsq <- init$d
    V <- init$v
    rm(init)
  }


  #  Update residuals (first iteration only)  -----------------------------
  M_obs <- partial_crossprod(U, t(t(V) *diag(Dsq)), irow, pcol, TRUE)
  Y@x <- Y@x - M_obs
  if (!is.null(warm_start)) {
    if (beta_flag) Y@x <- Y@x - xb_obs
    if (gamma_flag) Y@x <- Y@x - zg_obs
    if (intercept_row) add_to_rows_inplace_cpp(Y@x, Y@i, beta0, -1)
    if (intercept_col) add_to_cols_inplace_cpp(Y@x, Y@p, gamma0, -1)
  }

  #  Main loop ---------------------------------------------------------------
  ratio <- Inf
  iter <- 0
  while (ratio > thresh && iter < maxit) {
    iter <- iter + 1

    U_old <- U
    V_old <- V
    D_old <- Dsq
    # Intercepts (row/column) ---------------------------------------------
    # Row-level intercepts (beta0), then apply delta to residuals.
    if (intercept_row) {
      old_val <- beta0
      beta0 <- row_means_cpp(Y, nc) + beta0
      change <- old_val - beta0
      add_to_rows_inplace_cpp(Y@x, Y@i, change)
    }

    # Column-level intercepts (gamma0), then apply delta to residuals.
    if (intercept_col) {
      old_val <- gamma0
      gamma0 <- col_means_cpp(Y, nr) + gamma0
      change <- old_val - gamma0
      add_to_cols_inplace_cpp(Y@x, Y@p, change)
    }

    #  Update beta via soft-threshold --------------------------------------
    if (beta_flag) {
      beta <- soft_threshold_cpp(
        as.matrix((crossprod(X, Y)) + beta),
        lambda_beta
      )
      old_val <- xb_obs
      xb_obs <- partial_crossprod(X, beta, irow, pcol)
      Y@x <- Y@x + old_val - xb_obs
    }



    #  Update gamma via soft-threshold -------------------------------------
    if (gamma_flag) {
      gamma <- soft_threshold_cpp(
        as.matrix(Y %*% Z + gamma),
        lambda_gamma
      )

      old_val <- zg_obs
      zg_obs <- partial_crossprod(gamma, (Z), irow, pcol, TRUE)
      Y@x <- Y@x + old_val - zg_obs
    }

    #  Update (V, Dsq, U) from the "B" side --------------------------------
    # B_mat = BD
    if(laplace_c_flag){
      BD <- IMR:::update_B_sim_cpp(Y, U, V, Dsq, lambda_M, Uc, dc)
    }else
      BD <- IMR:::update_B_cpp(Y, U, V, Dsq, lambda_M)

    BD <- IMR:::svd_small_nc_cpp(BD)
    V <- BD$u
    Dsq <- tidyr::replace_na(BD$d, 0)
    U <- U %*% BD$v

    # update Y
    old_val <- M_obs
    M_obs <- partial_crossprod(U, t(t(V)*Dsq), irow, pcol, TRUE)
    Y@x <- Y@x + old_val - M_obs



    # 4.6 Update (U, Dsq, V) from the "A" side --------------------------------
    # A_mat <- AD
    if(laplace_r_flag){
      AD <- IMR:::update_A_sim_cpp(Y, U, V, Dsq, lambda_M, Ur, dr)
    }else
      AD <- IMR:::update_A_cpp(Y, U, V, Dsq, lambda_M)

    AD <- IMR:::svd_small_nc_cpp(AD)
    U <- AD$u
    Dsq <- tidyr::replace_na(AD$d, 0)
    V <- V %*% AD$v

    # update Y
    old_val <- M_obs
    M_obs <- partial_crossprod(U, t(t(V)*Dsq), irow, pcol, TRUE)
    Y@x <- Y@x + old_val - M_obs



    # 4.7 Convergence check ----------------------------------------------------
    ratio <- frob_ratio_cpp(U_old, D_old, V_old, U, Dsq, V)
    if (trace) {
      obj <- (0.5 * sum(Y@x^2) + lambda_M * sum(Dsq) +
        ifelse(beta_flag, lambda_beta * sum(abs(beta)), 0) +
        ifelse(gamma_flag, lambda_gamma * sum(abs(gamma)), 0)
      ) / nz
      cat(iter, " obj=", round(obj, 5), " ratio=", ratio, "\n")
    }
  }
  if (iter == maxit && trace) {
    warning("Did not converge in ", maxit, " iterations.")
  }

  # 5) Trim effective rank and return -----------------------------------------
  r_eff <- min(max(1, sum(Dsq > 0)), r)

  list(
    u = U[, seq_len(r_eff), drop = FALSE],
    d = Dsq[seq_len(r_eff)],
    v = V[, seq_len(r_eff), drop = FALSE],
    beta = beta,
    gamma = gamma,
    beta0 = beta0,
    gamma0 = gamma0,
    n_iter = iter
  )
}

#----------------------------------
#' Fit Incomplete Matrix Regression (IMR) Model without the low-rank \eqn{M}
#'
#' \code{imr.fit_no_low_rank} is similar to \code{imr.fit} except that it does not fit the low-rank matrix structure.
#'
#' @inheritParams imr.fit
#'
#' @return A list containing the learned parameters:
#' \describe{
#'   \item{resid}{An incomplete matrix of the last iteration's residuals (i.e., the model's training errors)}
#'   \item{beta}{Matrix of row covariate coefficients. \code{NULL} if \code{X} is \code{NULL}.}
#'   \item{gamma}{Matrix of column covariate coefficients. \code{NULL} if \code{Z} is \code{NULL}.}
#'   \item{beta0}{Vector of row-level intercepts. \code{NULL} if \code{intercept_row} is \code{FALSE}.}
#'   \item{gamma0}{Vector of column-level intercepts. \code{NULL} if \code{intercept_col} is \code{FALSE}.}
#'   \item{n_iter}{Integer. The number of iterations performed.}
#' }
#' @export
imr.fit_no_low_rank <- function(
    Y,
    X = NULL,
    Z = NULL,
    lambda_beta = NULL,
    lambda_gamma = NULL,
    intercept_row = FALSE,
    intercept_col = FALSE,
    maxit = 300,
    thresh = 1e-5,
    warm_start = NULL,
    trace = FALSE) {
  # Input checks & setup ----------------------------------------------------
  stopifnot(is.Incomplete(Y))

  dims <- dim(Y)
  nr <- dims[1]
  nc <- dims[2]
  nz <- Matrix::nnzero(Y, na.counted = TRUE)

  irow <- Y@i
  pcol <- Y@p


  # Laplacian flags (L_* expected as eigendecompositions) -------------------
  beta_flag <- !(is.null(lambda_beta) | is.null(X))
  gamma_flag <- !(is.null(lambda_gamma) | is.null(Z))
  # initial everything to null ------------------------
  beta <- gamma <- beta0 <- gamma0 <- NULL

  # 3) Warm-start or initialize ------------------------------------------------
  if(! is.null(warm_start)){
    if (beta_flag) {
      beta <- warm_start$beta
      xb_obs <- partial_crossprod(X, beta, irow, pcol)
    }
    if (gamma_flag) {
      gamma <- warm_start$gamma
      zg_obs <- partial_crossprod(gamma, Z, irow, pcol, TRUE)
    }
    if (intercept_row) {
      beta0 <- warm_start$beta0
    }

    if (intercept_col) {
      gamma0 <- warm_start$gamma0
    }
  }else{

  if (beta_flag) {
    beta <- matrix(0, ncol(X), nc)
    xb_obs <- rep(0, nz)
  }
  if (gamma_flag) {
    gamma <- matrix(0, nr, ncol(Z))
    zg_obs <- rep(0, nz)
  }
  if (intercept_row) {
    beta0 <- rep(0, nr)
  }
  if (intercept_col) {
    gamma0 <- rep(0, nc)
  }

  }
  if (!is.null(warm_start)) {
    if (beta_flag) Y@x <- Y@x - xb_obs
    if (gamma_flag) Y@x <- Y@x - zg_obs
    if (intercept_row) add_to_rows_inplace_cpp(Y@x, Y@i, beta0, -1)
    if (intercept_col) add_to_cols_inplace_cpp(Y@x, Y@p, gamma0, -1)
  }
  #  Main loop ---------------------------------------------------------------
  ratio <- Inf
  iter <- 0
  while (ratio > thresh && iter < maxit) {
    iter <- iter + 1
    old_err <- Y@x[]

    # Intercepts (row/column) ---------------------------------------------
    # Row-level intercepts (beta0), then apply delta to residuals.
    if (intercept_row) {
      old_val <- beta0
      beta0 <- row_means_cpp(Y, nc) + beta0
      change <- old_val - beta0
      add_to_rows_inplace_cpp(Y@x, Y@i, change)
    }

    # Column-level intercepts (gamma0), then apply delta to residuals.
    if (intercept_col) {
      old_val <- gamma0
      gamma0 <- col_means_cpp(Y, nr) + gamma0
      change <- old_val - gamma0
      add_to_cols_inplace_cpp(Y@x, Y@p, change)
    }




    #  Update gamma via soft-threshold -------------------------------------
    if (gamma_flag) {
      gamma <- soft_threshold_cpp(
        as.matrix(Y %*% Z + gamma),
        lambda_gamma
      )

      old_val <- zg_obs
      zg_obs <- partial_crossprod(gamma, (Z), irow, pcol, TRUE)
      Y@x <- Y@x + old_val - zg_obs
    }

    #  Update beta via soft-threshold --------------------------------------
    if (beta_flag) {
      beta <- soft_threshold_cpp(
        as.matrix((crossprod(X, Y)) + beta),
        lambda_beta
      )
      old_val <- xb_obs
      xb_obs <- partial_crossprod(X, beta, irow, pcol)
      Y@x <- Y@x + old_val - xb_obs
    }
    # 4.7 Convergence check ----------------------------------------------------
    ratio <- mean((Y@x - old_err)^2)

    if (trace) {
      obj <- (0.5 * sum(Y@x^2) +
        ifelse(beta_flag, lambda_beta * sum(abs(beta)), 0) +
        ifelse(gamma_flag, lambda_gamma * sum(abs(gamma)), 0)
      ) / nz
      cat(iter, " obj=", round(obj, 5), " ratio=", ratio, "\n")
    }
  }

  if (iter == maxit && trace) {
    warning("Did not converge in ", maxit, " iterations.")
  }

  #  return -----------------------------------------

  list(
    resid = Y,
    beta = beta,
    gamma = gamma,
    beta0 = beta0,
    gamma0 = gamma0,
    n_iter = iter
  )
}


#--------------------------------------

