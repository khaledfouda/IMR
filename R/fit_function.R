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
    L_a = NULL,
    lambda_a = 0,
    L_b = NULL,
    lambda_b = 0,
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
                                       intercept_col = intercept_col)
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

      init <- opt_svd(naive_MC(as.matrix(mfit$resid)), r, nr, nc, FALSE, FALSE)
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

      init <- opt_svd(naive_MC(as.matrix(Y)), r, nr, nc, FALSE, FALSE)
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


    B_mat <- update_B_cpp(Y, U, V, Dsq, lambda_M)
    B_mat <- svd_small_nc_cpp(B_mat)

    V <- B_mat$u
    Dsq <- B_mat$d
    U <- U %*% B_mat$v

    old_val <- M_obs
    M_obs <- partial_crossprod(U, V %*% diag(Dsq[,1],r,r), irow, pcol, TRUE)
    Y@x <- Y@x + old_val - M_obs



    # 4.6 Update (U, Dsq, V) from the "A" side --------------------------------

    A_mat <- update_A_cpp(Y, V, U, Dsq, lambda_M)
    A_mat <- svd_small_nc_cpp(A_mat)
    U <- A_mat$u
    Dsq <- A_mat$d
    V <- V %*% A_mat$v

    old_val <- M_obs
    M_obs <- partial_crossprod(U, V %*% diag(Dsq[,1],r,r), irow, pcol, TRUE)
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
    u            = U[, seq_len(r_eff), drop = FALSE],
    d            = Dsq[seq_len(r_eff)],
    v            = V[, seq_len(r_eff), drop = FALSE],
    beta         = beta,
    gamma        = gamma,
    beta0        = beta0,
    gamma0        = gamma0,
    n_iter       = iter
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

    # 4.7 Convergence check ----------------------------------------------------
    ratio <- mean((Y@x-old_err)^2)

    if (trace) {
      obj <- (0.5 * sum(Y@x^2)  +
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
    resid        = Y,
    beta         = beta,
    gamma        = gamma,
    beta0        = beta0,
    gamma0        = gamma0,
    n_iter       = iter
  )
}


#--------------------------------------
#'@export
error_metric = list(

#--- error functions:
unexplained_variance = function(predicted, true, adjusted = FALSE, k = NA) {
   # SSE / SST
   if(! adjusted){
      return(sum((true - predicted) ^ 2) / sum((true - mean(true)) ^ 2))
   }else{
      n = length(true)
      stopifnot(is.numeric(k))
      return(((sum((true - predicted) ^ 2) /
           sum((true - mean(true)) ^ 2)) *
            (n - 1) / (n - k - 1)))
   }
},

mape = function(predicted, true) {
   mean(abs((true - predicted) / true), na.rm = TRUE) * 100
},
mae = function(predicted, true) {
   mean(abs(true - predicted), na.rm = TRUE)
},

rmse_normalized = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE)) / sd(true, na.rm = TRUE)
},

rmse = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE))
},

spearman_R2 = function(predicted, true) {
   cor(true, predicted, method = "spearman")
}


)
