#' @export
prepare_data <- function(Y, X = NULL, Z = NULL,
                         similarity_rows = NULL,
                         similarity_cols = NULL,
                         val_prop = 0.2, seed = 2025) {
  out <- list(model = list(Xq = NULL, Zq = NULL, Xr = NULL, Zr = NULL))
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  out$model$Y <- out$Y <- IMR::as.Incomplete(Y)
  message("Performing train/valid split")
  obs_mask <- as.matrix(Y != 0)
  out$model$valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
  out$model$train_mask <- IMR::as.Incomplete(obs_mask * (1 - out$model$valid_mask))
  out$model$valid_mask <- IMR::as.Incomplete(out$model$valid_mask)
  out$model$y_train <- IMR::as.Incomplete(Y * out$model$train_mask)
  out$model$y_valid <- IMR::as.Incomplete(Y * out$model$valid_mask)
  rm(obs_mask)

  if (!is.null(similarity_rows)) {
    out$model$similarity_rows <- similarity_rows
  } else {
    out$model$similarity_rows <- diag(1, nrow(Y), nrow(Y))
  }
  if (!is.null(similarity_cols)) {
    out$model$similarity_cols <- similarity_cols
  } else {
    out$model$similarity_cols <- diag(1, ncol(Y), ncol(Y))
  }

  if (!is.null(X)) {
    stopifnot(is.matrix(X))
    stopifnot(nrow(X) == nrow(Y))
    xqr <- qr(as.matrix(X))
    out$X <- X
    out$model$Xq <- qr.Q(xqr)
    out$Xr <- out$model$Xr <- qr.R(xqr)
  }
  if (!is.null(Z)) {
    stopifnot(is.matrix(Z))
    stopifnot(nrow(Z) == ncol(Y))
    Zqr <- qr(as.matrix(Z))
    out$Z <- Z
    out$model$Zq <- qr.Q(Zqr)
    out$Zr <- out$model$Zr <- qr.R(Zqr)
  }

  return(out)
}

#' @export
reconstruct <- function(fit, dat, trace = TRUE) {
  out <- list(beta = NA, gamma = NA, M = NA, xbeta = NA, gammaz = NA, estimates = 0)

  check_mat <- function(mat, is_matrix = TRUE) {
    if (any(is.na(mat))) {
      return(FALSE)
    }
    if (is_matrix & is.matrix(mat)) {
      return(TRUE)
    }
    if (is_matrix) {
      return(FALSE)
    }
    if (is.vector(mat)) {
      return(TRUE)
    }
    return(FALSE)
  }
  #-----
  if (check_mat(fit$u) & check_mat(fit$d, FALSE) & check_mat(fit$v)) {
    if (trace) message("Constructing M ...")
    out$M <- fit$u %*% (fit$d * t(fit$v))
    out$estimates <- out$M
  }
  if (check_mat(fit$beta) & check_mat(dat$X) & check_mat(dat$Xr)) {
    if (trace) message("Constructing XBeta ...")
    out$beta <- solve(dat$Xr) %*% fit$beta
    out$xbeta <- dat$X %*% out$beta
    out$estimates <- out$estimates + out$xbeta
  }
  if (check_mat(fit$gamma) & check_mat(dat$Z) & check_mat(dat$Zr)) {
    if (trace) message("Constructing GammaZ ...")
    out$gamma <- fit$gamma %*% solve(t(dat$Zr))
    out$gammaz <- out$gamma %*% t(dat$Z)
    out$estimates <- out$estimates + out$gammaz
  }
  if (check_mat(fit$beta0, FALSE)) {
    if (trace) message("Constructing row intercepts ...")
    out$estimates <- out$estimates + fit$beta0 %*% matrix(1, 1, ncol(out$estimates))
  }
  if (check_mat(fit$gamma0, FALSE)) {
    if (trace) message("Constructing column intercepts ...")
    out$estimates <- out$estimates + matrix(1, nrow(out$estimates), 1) %*% t(fit$gamma0)
  }
  if (trace) message("done.")
  return(out)
}
#-----------------------------
#' @export
reconstruct_partial <- function(fit, dat, target, trace = FALSE) {
  stopifnot(is.Incomplete(target))
  if (trace) message("Constructing M ...")
  target@x <- IMR:::partial_crossprod(fit$u, fit$d * t(fit$v), target@i, target@p)

  check_mat <- function(mat, is_matrix = TRUE) {
    if (any(is.na(mat))) {
      return(FALSE)
    }
    if (is_matrix & is.matrix(mat)) {
      return(TRUE)
    }
    if (is_matrix) {
      return(FALSE)
    }
    if (is.vector(mat)) {
      return(TRUE)
    }
    return(FALSE)
  }

  if (check_mat(fit$beta) & check_mat(dat$X) & check_mat(dat$Xr)) {
    if (trace) message("Constructing XBeta ...")
    target@x <- target@x + partial_crossprod(X, fit$beta, target@i, target@p)
  }
  if (check_mat(fit$gamma) & check_mat(dat$Z) & check_mat(dat$Zr)) {
    if (trace) message("Constructing GammaZ ...")
    target@x <- target@x + partial_crossprod(fit$gamma, Z, target@i, target@p, TRUE)
  }
  if (check_mat(fit$beta0, FALSE)) {
    if (trace) message("Constructing row intercepts ...")
    add_to_rows_inplace_cpp(target@x, target@i, fit$beta0)
  }
  if (check_mat(fit$gamma0, FALSE)) {
    if (trace) message("Constructing column intercepts ...")
    add_to_cols_inplace_cpp(target@x, target@p, fit$gamma0)
  }
  return(target)
}
