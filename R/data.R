#' Prepare Data for IMR
#'
#' @export
prepare_data <- function(Y,
                         X = NULL,
                         Z = NULL,
                         similarity_rows = NULL,
                         similarity_cols = NULL,
                         val_prop = 0.2,
                         seed = 2025) {

  out <- list()
  if (!is.null(seed) && is.numeric(seed)) set.seed(seed)

  # ---  Target Matrix Setup ---
  out$Y <- IMR::as.Incomplete(Y)
  out$Y@x <- out$Y@x + 0 # this is to force a copy
  obs_mask <- as.matrix(Y != 0)

  dims <- dim(Y)
  n_total_obs <- sum(obs_mask)

  # --- Train/Validation Split ---
  split_data <- is.numeric(val_prop) && val_prop > 0

  if (split_data) {
    message("Performing train/valid split...")

    valid_mask_mat <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
    train_mask_mat <- obs_mask * (1 - valid_mask_mat)

    out$valid_mask <- IMR::as.Incomplete(valid_mask_mat)
    out$y_train <- IMR::as.Incomplete(Y * train_mask_mat)
    out$y_valid <- IMR::as.Incomplete(Y * valid_mask_mat)

    n_train <- sum(train_mask_mat)
    n_valid <- sum(valid_mask_mat)
  } else {
    n_train <- n_total_obs
    n_valid <- 0
  }

  out$obs_mask <- IMR::as.Incomplete(obs_mask * 1)

  # ---  Similarity Matrices ---
  if (!is.null(similarity_rows) && !is.null(similarity_rows$d) &&
      !is.null(similarity_rows$U)) {

    stopifnot(nrow(similarity_rows$U) == nrow(Y),
              length(similarity_rows$d) == ncol(similarity_rows$U))
    out$similarity_rows <- similarity_rows
  } else {
    out$similarity_rows <- NULL
  }

  if (!is.null(similarity_cols) && !is.null(similarity_cols$d) &&
      !is.null(similarity_cols$U)) {

    stopifnot(nrow(similarity_cols$U) == ncol(Y),
              length(similarity_cols$d) == ncol(similarity_cols$U))
    out$similarity_cols <- similarity_cols
  } else {
    out$similarity_cols <- NULL
  }

  # ---  Covariates ---
  if (!is.null(X)) {
    stopifnot(is.matrix(X), nrow(X) == nrow(Y))
    xqr <- qr(X)
    out$Xq <- qr.Q(xqr)
    out$Xr <- qr.R(xqr)
  }

  if (!is.null(Z)) {
    stopifnot(is.matrix(Z), nrow(Z) == ncol(Y))
    Zqr <- qr(Z)
    out$Zq <- qr.Q(Zqr)
    out$Zr <- qr.R(Zqr)
  }

  # ---  Metadata  ---
  out$meta <- list(
    dimensions = dims,
    sparsity = 1 - (n_total_obs / (dims[1] * dims[2])),
    total_obs = n_total_obs,
    split_data = split_data,
    n_train = n_train,
    n_valid = n_valid,
    val_prop = if (split_data) val_prop else 0,
    has_X = !is.null(X),
    num_X_vars = if (!is.null(X)) ncol(X) else 0,
    has_Z = !is.null(Z),
    num_Z_vars = if (!is.null(Z)) ncol(Z) else 0,
    has_sim_row = !is.null(out$similarity_rows),
    has_sim_col = !is.null(out$similarity_cols)
  )

  structure(out, class = "imr_data")
}
#' @export
print.imr_data <- function(x, ...) {
  m <- x$meta # Alias for cleaner code

  cat("\n== IMR Data Object ==\n")

  # Base Dimensions
  cat(sprintf("Target Matrix (Y): %d rows x %d cols\n", m$dimensions[1], m$dimensions[2]))
  cat(sprintf("Observed Entries:  %d (%.2f%% Sparsity)\n",
              m$total_obs, m$sparsity * 100))

  # Train/Valid Split info
  if (m$split_data) {
    cat(sprintf("  - Training:      %d (%.1f%%)\n", m$n_train, 100 * (1 - m$val_prop)))
    cat(sprintf("  - Validation:    %d (%.1f%%)\n", m$n_valid, 100 * m$val_prop))
  } else {
    cat("  - Training:      Using 100% of data (No validation split)\n")
  }

  # Covariates
  cat("\n-- Covariates --\n")
  cat(sprintf("Row Covariates (X): %s\n",
              if (m$has_X) sprintf("%d variables", m$num_X_vars) else "[None]"))
  cat(sprintf("Col Covariates (Z): %s\n",
              if (m$has_Z) sprintf("%d variables", m$num_Z_vars) else "[None]"))

  # Similarities
  cat("\n-- Similarity Matrices (Decomposed) --\n")
  cat(sprintf("Row Similarity: %s\n", if (m$has_sim_row) "Provided" else "[None]"))
  cat(sprintf("Col Similarity: %s\n", if (m$has_sim_col) "Provided" else "[None]"))

  cat("=====================\n")
  invisible(x)
}

#' @export
prepare_data <- function(Y, X = NULL, Z = NULL,
                         similarity_rows = NULL,
                         similarity_cols = NULL,
                         val_prop = 0.2, seed = 2025) {
  out <- list()
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  out$Y <- IMR::as.Incomplete(Y)
  message("Performing train/valid split")
  obs_mask <- as.matrix(Y != 0)
  out$valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
  train_mask <- IMR::as.Incomplete(obs_mask * (1 - out$valid_mask))
  out$valid_mask <- IMR::as.Incomplete(out$valid_mask)
  out$obs_mask <- IMR::as.Incomplete(obs_mask * 1)

  out$y_train <- IMR::as.Incomplete(Y * train_mask)
  out$y_valid <- IMR::as.Incomplete(Y * out$valid_mask)


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
