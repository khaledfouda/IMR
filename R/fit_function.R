#' @export
imr_control <- function(maxit = 300,
                        thresh = 1e-5,
                        trace = FALSE,
                        ls_initial = TRUE) {
  structure(
    list(maxit = maxit,
         thresh = thresh,
         trace = trace,
         ls_initial = ls_initial),
    class = "imr_control"
  )
}

#' @export
fit_imr <- function(
    data,
    rank = 2,
    lambda_M = 0,
    lambda_beta = 0,
    lambda_gamma = 0,
    intercept_row = FALSE,
    intercept_col = FALSE,
    shared_beta = FALSE,
    shared_gamma = FALSE,
    control = imr_control(),
    warm_start = NULL){

  # validation
  if (!inherits(data, "imr_data")) stop("Data must be an 'imr_data' object.")
  if (!inherits(control, "imr_control")) control <- imr_control()


  if (inherits(warm_start, "imr_fit")) {
    warm_start <- warm_start$coefficients
  }



  result_list <- imr_solver(
    Y = data$Y,
    X = data$X,
    Z = data$Z,
    r = rank,
    lambda_M = lambda_M,
    lambda_beta = lambda_beta,
    lambda_gamma = lambda_gamma,
    Ur = data$similarity_row$U,
    dr = data$similarity_row$d,
    Uc = data$similarity_col$U,
    dc = data$similarity_col$d,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    shared_beta = shared_beta,
    shared_gamma = shared_gamma,
    control = control,
    warm_start = warm_start
  )

  # 5. Construct Final Object
  structure(
    list(
      coefficients = list(
        u = result_list$u,
        d = result_list$d,
        v = result_list$v,
        beta = result_list$beta,
        gamma = result_list$gamma,
        beta0 = result_list$beta0,
        gamma0 = result_list$gamma0
      ),
      residuals = result_list$residuals,
      meta = list(
        rank = rank,
        lambdas = c(M=lambda_M, beta=lambda_beta, gamma=lambda_gamma),
        intercepts = c(row = intercept_row, col = intercept_col),
        shared_effects = c(beta=shared_beta, gamma=shared_gamma),
        n_iter = result_list$n_iter,
        converged = result_list$n_iter < control$maxit,
        # statistic for print function
        tss = sum((data$Y@x - mean(data$Y@x))^2)
      ),
      control = control
    ),
    class = "imr_fit"
  )
}


#' Fit Incomplete Matrix Regression (IMR) Model
#'
#' Not exported
imr_solver <- function(
  Y, X, Z,
  intercept_row, intercept_col,
  shared_beta, shared_gamma,
  r, lambda_M, lambda_beta, lambda_gamma,
  Ur, dr, Uc, dc,
  control,
  warm_start) {
  # Input checks & setup ----------------------------------------------------
  stopifnot(is.Incomplete(Y))
  dims <- dim(Y)
  nr <- dims[1]
  nc <- dims[2]
  nz <- Matrix::nnzero(Y, na.counted = TRUE)
  irow <- Y@i
  pcol <- Y@p
  Y@x <- Y@x + 0 # force a copy so the original matrix doesn't get modified by C++

  # Unpack control ----------------------------------------------------------
  maxit  <- control$maxit
  thresh <- control$thresh
  trace  <- control$trace
  ls_initial <- control$ls_initial
  shared_beta <- control$shared_beta
  shared_gamma <- control$shared_gamma

  # Laplacian flags (L_* expected as eigendecompositions) -------------------
  beta_flag <- !(is.null(X))
  gamma_flag <- !(is.null(Z))
  laplace_r_flag <- !(is.null(Ur) || is.null(dr))
  laplace_c_flag <- !(is.null(Uc) || is.null(dc))
  low_rank_flag <- !(is.null(r) || is.null(lambda_M) || r <= 0)
  if(!low_rank_flag) ls_initial = FALSE
  #-------------------------------------------------
  if(laplace_r_flag) dr = dr * lambda_M
  if(laplace_c_flag) dc = dc * lambda_M
  #--------------------------------------------------
  # initial everything to null ------------------------
  beta <- gamma <- beta0 <- gamma0 <- U <- V <- Dsq <- NULL

  # 3) Warm-start or initialize ------------------------------------------------
  warm_start <- verify_warm_start(warm_start, r)

  if (!is.null(warm_start)) {
    required <- c("u", "d", "v", "beta", "gamma", "beta0", "gamma0")
    if (!all(required %in% names(warm_start))) {
      stop("warm_start missing components: u, d, v, beta, gamma, beta0, gamma0")
    }

    if (beta_flag) {
      beta <- warm_start$beta
      if(shared_beta){
        xbeta <- X %*% beta
      }else{
        xb_obs <- partial_crossprod(X, beta, irow, pcol)
      }
    }
    if (gamma_flag) {
      gamma <- warm_start$gamma
      if(shared_gamma){
        gammaz <- tcrossprod(gamma, Z)
      }else{
        zg_obs <- partial_crossprod(gamma, Z, irow, pcol, TRUE)
      }
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
      mift <- imr_solver(
        Y = Y, X = X, Z = Z,
        intercept_row = intercept_row, intercept_col = intercept_col,
        shared_beta = shared_beta, shared_gamma = shared_gamma,
        r = 0, lambda_M = NULL, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma,
        Ur = NULL, dr = NULL, Uc = NULL, dc = NULL,
        control = control,
        warm_start = NULL)

      if (beta_flag) {
        beta <- mfit$beta
        if(shared_beta){
          xbeta <- X %*% beta
        }else{
          xb_obs <- partial_crossprod(X, beta, irow, pcol)
        }
      }
      if (gamma_flag) {
        gamma <- mfit$gamma
        if(shared_gamma){
          gammaz <- tcrossprod(gamma, Z)
        }else{
          zg_obs <- partial_crossprod(gamma, Z, irow, pcol, TRUE)
        }
      }
      if (intercept_row) {
        beta0 <- mfit$beta0
      }
      if (intercept_col) {
        gamma0 <- mfit$gamma0
      }

      init <- IMR::opt_svd(mfit$residuals, r, nr, nc, FALSE, FALSE)
    } else {
      if (beta_flag) {
        if(shared_beta){
          beta <- rep(0, ncol(X))
          xbeta <- rep(0, nr)
        }else{
          beta <- matrix(0, ncol(X), nc)
          xb_obs <- rep(0, nz)
        }
      }
      if (gamma_flag) {
        if(shared_gamma){
          gamma <- rep(0, ncol(Z))
          gammaz <- rep(0, nc)
        }else{
          gamma <- matrix(0, nr, ncol(Z))
          zg_obs <- rep(0, nz)
        }
      }
      if (intercept_row) {
        beta0 <- rep(0, nr)
      }
      if (intercept_col) {
        gamma0 <- rep(0, nc)
      }
      if(low_rank_flag)
        init <- IMR::opt_svd(Y, r, nr, nc, FALSE, FALSE)
    }

    if(low_rank_flag){
      U <- init$u
      Dsq <- init$d
      V <- init$v
      rm(init)
    }
  }


  #  Update residuals (first iteration only)  -----------------------------
  if(low_rank_flag){
    M_obs <- partial_crossprod(U, t(t(V) * Dsq), irow, pcol, TRUE)
    Y@x <- Y@x - M_obs
  }
  if (!is.null(warm_start) || (is.null(warm_start) && ls_initial)) {
    if (beta_flag && !shared_beta) Y@x <- Y@x - xb_obs
    if (gamma_flag && !shared_gamma) Y@x <- Y@x - zg_obs
    if (beta_flag && shared_beta) add_to_rows_inplace_cpp(Y@x, irow, -xbeta)
    if (gamma_flag && shared_gamma) add_to_cols_inplace_cpp(Y@x, pcol, -gammaz)
    if (intercept_row) add_to_rows_inplace_cpp(Y@x, irow, -beta0)
    if (intercept_col) add_to_cols_inplace_cpp(Y@x, pcol, -gamma0)
  }

  #  Main loop ---------------------------------------------------------------
  ratio <- Inf
  iter <- 0
  while (ratio > thresh && iter < maxit) {
    iter <- iter + 1

    U_old <- U
    V_old <- V
    D_old <- Dsq
    beta_old <- beta
    gamma_old <- gamma

    if(low_rank_flag){
      #  Update (V, Dsq, U) from the "B" side --------------------------------
      # B_mat = BD
      if (laplace_c_flag) {
        BD <- IMR:::update_B_sim_cpp(Y, U, V, Dsq, lambda_M, Uc, dc)
      } else {
        BD <- IMR:::update_B_cpp(Y, U, V, Dsq, lambda_M)
      }

      BD <- IMR:::svd_small_nc_cpp(BD)
      V <- BD$u
      Dsq <- tidyr::replace_na(BD$d, 0)
      U <- U %*% BD$v

      # update Y
      old_val <- M_obs
      M_obs <- partial_crossprod(U, t(t(V) * Dsq), irow, pcol, TRUE)
      Y@x <- Y@x + old_val - M_obs


      # 4.6 Update (U, Dsq, V) from the "A" side --------------------------------
      # A_mat <- AD
      if (laplace_r_flag) {
        AD <- IMR:::update_A_sim_cpp(Y, U, V, Dsq, lambda_M, Ur, dr)
      } else {
        AD <- IMR:::update_A_cpp(Y, U, V, Dsq, lambda_M)
      }

      AD <- IMR:::svd_small_nc_cpp(AD)
      U <- AD$u
      Dsq <- tidyr::replace_na(AD$d, 0)
      V <- V %*% AD$v

      # update Y
      old_val <- M_obs
      M_obs <- partial_crossprod(U, t(t(V) * Dsq), irow, pcol, TRUE)
      Y@x <- Y@x + old_val - M_obs
    }

    #  Update beta via soft-threshold --------------------------------------
    if (beta_flag) {
      if(shared_beta){
        beta <- soft_threshold_cpp(
          crossprod(X, row_means_cpp(Y, nc) ) + beta,
          lambda_beta
        )
        old_val <- xbeta
        xbeta <- X %*% beta
        change <- old_val - xbeta
        add_to_rows_inplace_cpp(Y@x, irow, change)
      }else{
        beta <- soft_threshold_cpp(
          as.matrix((crossprod(X, Y)) + beta),
          lambda_beta
        )
        old_val <- xb_obs
        xb_obs <- partial_crossprod(X, beta, irow, pcol)
        Y@x <- Y@x + old_val - xb_obs
      }
    }


    #  Update gamma via soft-threshold -------------------------------------
    if (gamma_flag) {
      if(shared_gamma){
        gamma <- soft_threshold_cpp(
          col_means_cpp(Y, nr) %*% Z + gamma,
          lambda_gamma
        )
        old_val <- gammaz
        gammaz <- tcrossprod(gamma, Z)
        change <- old_val - gammaz
        add_to_cols_inplace_cpp(Y@x, pcol, change)
      }else{
        gamma <- soft_threshold_cpp(
          as.matrix(Y %*% Z + gamma),
          lambda_gamma
        )
        old_val <- zg_obs
        zg_obs <- partial_crossprod(gamma, (Z), irow, pcol, TRUE)
        Y@x <- Y@x + old_val - zg_obs
      }
    }


    # Intercepts (row/column) ---------------------------------------------
    # Row-level intercepts (beta0), then apply delta to residuals.
    if (intercept_row) {
      old_val <- beta0
      beta0 <- row_means_cpp(Y, nc) + beta0
      change <- old_val - beta0
      add_to_rows_inplace_cpp(Y@x, irow, change)
    }

    # Column-level intercepts (gamma0), then apply delta to residuals.
    if (intercept_col) {
      old_val <- gamma0
      gamma0 <- col_means_cpp(Y, nr) + gamma0
      change <- old_val - gamma0
      add_to_cols_inplace_cpp(Y@x, pcol, change)
    }

    # 4.7 Convergence check ----------------------------------------------------
    ratio <- 0
    if(low_rank_flag)
      ratio <- frob_ratio_cpp(U_old, D_old, V_old, U, Dsq, V)
    if(beta_flag){
      denom = sum(beta_old*beta_old)
      if(denom > 0)
        ratio <- ratio + sum((beta_old-beta)^2) / denom
    }
    if(gamma_flag){
      denom = sum(gamma_old*gamma_old)
      if(denom > 0)
        ratio <- ratio + sum((gamma_old-gamma)^2) / denom
    }
    if (trace) {
      obj <- (0.5 * sum(Y@x^2) +
        ifelse(low_rank_flag, lambda_M * sum(Dsq), 0) +
        ifelse(beta_flag, lambda_beta * sum(abs(beta)), 0) +
        ifelse(gamma_flag, lambda_gamma * sum(abs(gamma)), 0)
      ) / nz
      cat(iter, " obj=", round(obj, 5), " ratio=", ratio, "\n")
    }
  }
  if (iter == maxit) {
    warning("Did not converge in ", maxit, " iterations.")
  }

  # 5) Trim effective rank and return -----------------------------------------
  if(low_rank_flag){
    r_eff <- min(max(1, sum(Dsq > 0)), r)
    U = U[, seq_len(r_eff), drop = FALSE]
    Dsq = Dsq[seq_len(r_eff)]
    V = V[, seq_len(r_eff), drop = FALSE]
  }

  list(
    residuals = Y,
    u = u,
    d = Dsq,
    v = V,
    beta = beta,
    gamma = gamma,
    beta0 = beta0,
    gamma0 = gamma0,
    n_iter = iter
  )
}



#' @export
print.imr_fit <- function(x, ...) {

  cat("\n== IMR Fitted Model ==\n")

  has_low_rank <- !is.null(x$coefficients$d) && length(x$coefficients$d) > 0
  has_row_cov  <- !is.null(x$coefficients$beta)
  has_col_cov  <- !is.null(x$coefficients$gamma)
  has_row_int  <- !is.null(x$coefficients$beta0)
  has_col_int  <- !is.null(x$coefficients$gamma0)

  terms <- c()
  if (has_row_cov) {
    beta_type <- if (x$meta$shared_effects["beta"]) "vec" else "mat"
    terms <- c(terms, paste0("X\u00b7\u03b2(", beta_type, ")"))
  }
  if (has_col_cov) {
    gamma_type <- if (x$meta$shared_effects["gamma"]) "vec" else "mat"
    terms <- c(terms, paste0("\u03b3(", gamma_type, ")\u00b7Z'"))
  }
  if (has_low_rank) terms <- c(terms, "M(SVD)")
  if (has_row_int)  terms <- c(terms, "\u03b2\u2080")
  if (has_col_int)  terms <- c(terms, "\u03b3\u2080")
  if (length(terms) == 0) terms <- c("0")

  cat(sprintf("Equation:  Y ~ %s\n", paste(terms, collapse = " + ")))

  if (!is.null(x$residuals)) {
    rss <- sum(x$residuals@x^2)
    tss <- x$meta$tss
    n   <- x$meta$n_obs

    mse  <- rss / n
    rmse <- sqrt(mse)
    r2   <- 1 - (rss / tss)

    # Ensure R2 isn't negative (possible in non-OLS or heavy regularization)
    r2_str <- if (r2 < 0) "< 0 (Poor Fit)" else sprintf("%.4f", r2)

    cat("\n-- Fit Statistics --\n")
    cat(sprintf("RMSE:      %.5f\n", rmse))
    cat(sprintf("MSE:       %.5f\n", mse))
    cat(sprintf("R-squared: %s\n", r2_str))
  }

  # --- 4. Hyperparameters ---
  cat("\n-- Hyperparameters --\n")
  if (has_low_rank) {
    cat(sprintf("Rank (r):     %d\n", x$meta$rank))
    cat(sprintf("Lambda M:     %g\n", x$meta$lambdas["M"]))
  }
  if (has_row_cov) cat(sprintf("Lambda Beta:  %g\n", x$meta$lambdas["beta"]))
  if (has_col_cov) cat(sprintf("Lambda Gamma: %g\n", x$meta$lambdas["gamma"]))

  # --- 5. Convergence ---
  status <- if (x$meta$converged) "Converged" else "Did NOT converge"
  cat(sprintf("\nStatus: %s (in %d iterations)\n", status, x$meta$n_iter))
  cat("======================\n")

  invisible(x)
}

#' @export
print.imr_control <- function(x, ...) {
  cat("\n== IMR Convergence Parameters ==\n")

  cat(sprintf("Max Iterations: %d\n", x$maxit))
  cat(sprintf("Threshold:      %g\n", x$thresh))

  init_method <- if (x$ls_initial) "Least Squares" else "Randomized"
  cat(sprintf("Initialization: %s\n", init_method))

  # 3. Trace/Verbosity
  trace_status <- if (x$trace) "Enabled" else "Disabled"
  cat(sprintf("Trace Progress: %s\n", trace_status))

  cat("================================\n")
  invisible(x)
}
