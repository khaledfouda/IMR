#' @export
imr.fit <- function(
    Y,
    X = NULL,
    Z = NULL,
    r = 2,
    lambda_M = 0,
    lambda_beta = 0,
    lambda_gamma = 0,
    intercept_row = FALSE,
    intercept_col = FALSE,
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
  M_obs <- partial_crossprod(U, V %*% diag(Dsq), irow, pcol, TRUE)
  Y@x <- Y@x - M_obs
  if (!is.null(warm_start)) {
    if (beta_flag) Y@x <- Y@x - xb_obs
    if (gamma_flag) Y@x <- Y@x - zg_obs
    if (intercept_row) add_to_rows_inplace_cpp(Y@x, Y@i, beta0, -1)
    if (intercept_col) add_to_cols_inplace_cpp(Y@x, Y@p, gamma0, -1)
  }
  B_mat <- matrix(NA, ncol(Y), r)
  A_mat <- matrix(NA, nrow(Y), r)

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
      partial = crossprod(Y, U) + sweep(V, 2L, Dsq, `*`) # V %*% Dsq
      partial = crossprod(Uc, sweep(partial, 2L, Dsq, `*`))
      # partial = crossprod(Uc, partial %*% diag(Dsq))
      coef = 1 / outer(dc + lambda_M, Dsq, `+`)
      B_mat = Uc %*% (partial * coef)
      # for(j in seq_len(r)){
      #   B_mat[,j] <- Uc %*% diag(1/(dc+Dsq[j]+lambda_M)) %*% partial[,j]
      # }
    }else{
      B_mat <- update_B_cpp(Y, U, V, Dsq, lambda_M) # output: DB
    }
      BD_decomp <- svd_small_nc_cpp(B_mat)

    Dsq <- trim_eig(BD_decomp$d, eig.tol)
    numEig <- length(Dsq)
    V <- BD_decomp$u[, seq_len(numEig), drop = FALSE]
    U <- U[, seq_len(numEig), drop = FALSE] %*% BD_decomp$v[seq_len(numEig), seq_len(numEig), drop = FALSE]
    # V <- BD_decomp$u
    # Dsq <- BD_decomp$d
    # U <- U %*% BD_decomp$v

    old_val <- M_obs
    M_obs <- partial_crossprod(U, V %*% diag(Dsq, numEig, numEig), irow, pcol, TRUE)
    Y@x <- Y@x + old_val - M_obs



    # 4.6 Update (U, Dsq, V) from the "A" side --------------------------------
    # A_mat <- AD
    if(laplace_r_flag){
      partial = Y %*% V + sweep(U, 2L, Dsq, `*`)  #U %*% diag(Dsq)
      partial = crossprod(Ur, sweep(partial, 2L, Dsq, `*`))
      # partial = crossprod(Ur, partial %*% diag(Dsq))
      coef = 1 / outer(dr + lambda_M, Dsq, `+`)
      A_mat = Ur %*% (partial * coef)
      # for(j in seq_len(r)){
      #   A_mat[,j] <- Ur %*% diag(1/(dr+Dsq[j]+lambda_M)) %*% partial[,j]
      # }
    }else{
      # output here is
      A_mat <- update_A_cpp(Y, V, U, Dsq, lambda_M)
    }
    AD_decomp <- svd_small_nc_cpp(A_mat)

    # Dsq <- A_mat$d[A_mat$d > 0]
    Dsq <- trim_eig(AD_decomp$d, eig.tol)
    numEig <- length(Dsq)
    U <- AD_decomp$u[, seq_len(numEig), drop = FALSE]
    V <- V[, seq_len(numEig), drop = FALSE] %*% AD_decomp$v[seq_len(numEig), seq_len(numEig), drop = FALSE]

    old_val <- M_obs
    M_obs <- partial_crossprod(U, V %*% diag(Dsq, numEig, numEig), irow, pcol, TRUE)
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

