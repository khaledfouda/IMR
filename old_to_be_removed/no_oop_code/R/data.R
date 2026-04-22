#' @export
prepare_data <- function(Y, X = NULL, Z = NULL,
                         similarity_rows = NULL,
                         similarity_cols = NULL,
                         val_prop = 0.2, seed = 2025) {
  out <- list()
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  out$Y <- IMR::as_incomplete(Y)
  message("Performing train/valid split")
  obs_mask <- as.matrix(Y != 0)
  out$valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
  train_mask <- IMR::as_incomplete(obs_mask * (1 - out$valid_mask))
  out$valid_mask <- IMR::as_incomplete(out$valid_mask)
  out$obs_mask <- IMR::as_incomplete(obs_mask * 1)

  out$y_train <- IMR::as_incomplete(Y * train_mask)
  out$y_valid <- IMR::as_incomplete(Y * out$valid_mask)


  if (!is.null(similarity_rows)) {
    stopifnot(nrow(similarity_rows) == nrow(Y))
    stopifnot(isSymmetric(similarity_rows))
    stopifnot(nrow(similarity_rows) == ncol(similarity_rows))
    out$similarity_rows <- similarity_rows
  } else {
    out$similarity_rows <- NULL
  }
  if (!is.null(similarity_cols)) {
    stopifnot(nrow(similarity_cols) == ncol(Y))
    stopifnot(isSymmetric(similarity_cols))
    stopifnot(nrow(similarity_cols) == ncol(similarity_cols))
    out$similarity_cols <- similarity_cols
  } else {
    out$similarity_cols <- NULL
  }

  if (!is.null(X)) {
    stopifnot(is.matrix(X))
    stopifnot(nrow(X) == nrow(Y))
    xqr <- qr(X)
    #out$X <- X
    out$Xq <- qr.Q(xqr)
    out$Xr <- out$Xr <- qr.R(xqr)
  }
  if (!is.null(Z)) {
    stopifnot(is.matrix(Z))
    stopifnot(nrow(Z) == ncol(Y))
    Zqr <- qr(Z)
    #out$Z <- Z
    out$Zq <- qr.Q(Zqr)
    out$Zr <- out$Zr <- qr.R(Zqr)
  }

  return(out)
}

#' @export
reconstruct <- function(fit, data, trace = TRUE, shared_information=FALSE) {
  out <- list(beta = NA, gamma = NA, M = NA, xbeta = NA, gammaz = NA, estimates = 0)
  # remove this check later
  if(shared_information){
    if(!is.null(fit$beta)) fit$beta <- as.vector(fit$beta)
    if(!is.null(fit$gamma)) fit$gamma <- as.vector(fit$gamma)
  }
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
  if (check_mat(data$Xq) && check_mat(data$Xr)) {
    if(!shared_information && check_mat(fit$beta)){
      if (trace) message("Constructing XBeta ...")
      out$beta <- solve(data$Xr) %*% fit$beta
      out$xbeta <- data$Xq %*% fit$beta
      out$estimates <- out$estimates + out$xbeta
    }else if(shared_information && check_mat(fit$beta, FALSE)){
      if (trace) message("Constructing XBeta ...")
      out$beta <- solve(data$Xr) %*% fit$beta
      out$xbeta <- data$Xq %*% fit$beta
      out$estimates <- out$estimates + out$xbeta  %*% matrix(1, 1, ncol(out$estimates))
    }
  }
  if (check_mat(data$Zq) && check_mat(data$Zr)) {
    if(!shared_information && check_mat(fit$gamma)){
      if (trace) message("Constructing GammaZ ...")
      out$gamma <- fit$gamma %*% solve(t(data$Zr))
      out$gammaz <- fit$gamma %*% t(data$Zq)
      out$estimates <- out$estimates + out$gammaz
    }else if(shared_information && check_mat(fit$gamma, FALSE)){
      if (trace) message("Constructing GammaZ ...")
      out$gamma <- fit$gamma %*% solve(t(data$Zr))
      out$gammaz <- fit$gamma %*% t(data$Zq)
      out$estimates <- out$estimates + matrix(1, nrow(out$estimates), 1) %*% out$gammaz

    }
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
reconstruct_partial <- function(fit, data, target, shared_information = FALSE, trace = FALSE) {
  stopifnot(is_incomplete(target))
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

  if (check_mat(data$Xq) && check_mat(data$Xr)) {
    if (trace) message("Constructing XBeta ...")
    if(shared_information){
      xbeta <- data$Xq %*% fit$beta
      add_to_rows_inplace_cpp(target@x, target@i, xbeta)
    }else
      target@x <- target@x + partial_crossprod(data$Xq, fit$beta, target@i, target@p)
  }
  if (check_mat(data$Zq) && check_mat(data$Zr)) {
    if (trace) message("Constructing GammaZ ...")
    if(shared_information){
      gammaz <- tcrossprod(fit$gamma, data$Zq)
      add_to_cols_inplace_cpp(target@x, target@p, gammaz)
    }else
    target@x <- target@x + partial_crossprod(fit$gamma, data$Zq, target@i, target@p, TRUE)
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
