#' Prepare Data for IMR
#'
#' @export
imr_data <- function(Y,
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
    stopifnot(
      nrow(similarity_rows$U) == nrow(Y),
      length(similarity_rows$d) == ncol(similarity_rows$U)
    )
    out$similarity_rows <- similarity_rows
  } else {
    out$similarity_rows <- NULL
  }

  if (!is.null(similarity_cols) && !is.null(similarity_cols$d) &&
    !is.null(similarity_cols$U)) {
    stopifnot(
      nrow(similarity_cols$U) == ncol(Y),
      length(similarity_cols$d) == ncol(similarity_cols$U)
    )
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
  # --- Model Structure ----
  out$model <- list(
    with_row_covariates = out$meta$has_X,
    with_col_covariates = out$meta$has_Z,
    with_low_rank_component = TRUE,
    with_row_similarity = out$meta$has_sim_row,
    with_col_similarity = out$meta$has_sim_col
  )

  structure(out, class = "imr_data")
}
#----------------------------------------
update.imr_data <- function(object,
                            with_row_covariates = NULL,
                            with_col_covariates = NULL,
                            with_low_rank_component = NULL,
                            with_row_similarity = NULL,
                            with_col_similarity = NULL,
                            ...) {

  update_flag <- function(obj, flag_val, flag_name, dependency_name) {
    if (!is.null(flag_val)) {
      # Ensure it's a single boolean
      flag_val <- as.logical(flag_val[1])

      if (flag_val && !obj$meta[[dependency_name]]) {
        stop(sprintf("Cannot set '%s = TRUE' because the underlying data ('%s') was not provided.",
                     flag_name, dependency_name), call. = FALSE)
      }
      obj$model[[flag_name]] <- flag_val
    }
    return(obj)
  }

  object <- update_flag(object, with_row_covariates, "with_row_covariates", "has_X")
  object <- update_flag(object, with_col_covariates, "with_col_covariates", "has_Z")
  object <- update_flag(object, with_row_similarity, "with_row_similarity", "has_sim_row")
  object <- update_flag(object, with_col_similarity, "with_col_similarity", "has_sim_col")

  if (!is.null(with_low_rank_component)) {
    object$model$with_low_rank_component <- with_low_rank_component
  }

  return(object)
}


#' @export
print.imr_data <- function(x, ...) {
  m <- x$meta         # What data is available in memory
  a <- x$model # What the solver is instructed to use

  cat("\n== IMR Data Object ==\n")

  # --- 1. Base Dimensions ---
  cat(sprintf("Target Matrix (Y): %d rows x %d cols\n", m$dimensions[1], m$dimensions[2]))
  cat(sprintf("Observed Entries:  %d (%.2f%% Sparsity)\n", m$total_obs, m$sparsity * 100))

  # --- 2. Train/Valid Split info ---
  if (m$split_data) {
    cat(sprintf("  - Training:      %d (%.1f%%)\n", m$n_train, 100 * (1 - m$val_prop)))
    cat(sprintf("  - Validation:    %d (%.1f%%)\n", m$n_valid, 100 * m$val_prop))
  } else {
    cat("  - Training:      Using 100% of data (No validation split)\n")
  }

  # --- 3. Active Submodel Configuration ---
  cat("\n-- Model Configuration --\n")

  # Helper function to align text and format the active/inactive tags
  format_status <- function(has_data, is_active, data_desc) {
    if (!has_data) return(sprintf("[None]"))
    status_tag <- if (is_active) "[ACTIVE]" else ""
    return(sprintf("%-15s %s", data_desc, status_tag))
  }

  # Low-rank is a purely algorithmic component (doesn't depend on external data)
  cat(sprintf("%-20s:  %23s\n", "Low-Rank Matrix (M)",
              if (a$with_low_rank_component) "[ACTIVE]" else ""))

  # Covariates
  cat(sprintf("%-20s: %s\n", "Row Covariates (X)",
              format_status(m$has_X, a$with_row_covariates, sprintf("%d vars", m$num_X_vars))))

  cat(sprintf("%-20s: %s\n", "Col Covariates (Z)",
              format_status(m$has_Z, a$with_col_covariates, sprintf("%d vars", m$num_Z_vars))))

  # Similarities
  cat(sprintf("%-20s: %s\n", "Row Similarity",
              format_status(m$has_sim_row, a$with_row_similarity, "Provided")))

  cat(sprintf("%-20s: %s\n", "Col Similarity",
              format_status(m$has_sim_col, a$with_col_similarity, "Provided")))

  cat("==========================\n")
  invisible(x)
}
# print.imr_data <- function(x, ...) {
#   m <- x$meta # Alias for cleaner code
#
#   cat("\n== IMR Data Object ==\n")
#
#   # Base Dimensions
#   cat(sprintf("Target Matrix (Y): %d rows x %d cols\n", m$dimensions[1], m$dimensions[2]))
#   cat(sprintf(
#     "Observed Entries:  %d (%.2f%% Sparsity)\n",
#     m$total_obs, m$sparsity * 100
#   ))
#
#   # Train/Valid Split info
#   if (m$split_data) {
#     cat(sprintf("  - Training:      %d (%.1f%%)\n", m$n_train, 100 * (1 - m$val_prop)))
#     cat(sprintf("  - Validation:    %d (%.1f%%)\n", m$n_valid, 100 * m$val_prop))
#   } else {
#     cat("  - Training:      Using 100% of data (No validation split)\n")
#   }
#
#   # Covariates
#   cat("\n-- Covariates --\n")
#   cat(sprintf(
#     "Row Covariates (X): %s\n",
#     if (m$has_X) sprintf("%d variables", m$num_X_vars) else "[None]"
#   ))
#   cat(sprintf(
#     "Col Covariates (Z): %s\n",
#     if (m$has_Z) sprintf("%d variables", m$num_Z_vars) else "[None]"
#   ))
#
#   # Similarities
#   cat("\n-- Similarity Matrices (Decomposed) --\n")
#   cat(sprintf("Row Similarity: %s\n", if (m$has_sim_row) "Provided" else "[None]"))
#   cat(sprintf("Col Similarity: %s\n", if (m$has_sim_col) "Provided" else "[None]"))
#
#   cat("=====================\n")
#   invisible(x)
# }
#'


#' @export
reconstruct <- function(fit, data, trace = TRUE) {
  stopifnot(inherits(fit, "imr_fit"))
  stopifnot(inherits(data, "imr_data"))

  out <- list(
    beta = NULL, gamma = NULL, M = NULL,
    xbeta = NULL, gammaz = NULL, estimates = 0
  )

  coefs <- fit$coefficients
  meta <- fit$meta

  # Get dimensions to safely initialize estimates if M is absent
  n <- data$meta$dimensions[1]
  m <- data$meta$dimensions[2]
  out$estimates <- matrix(0, nrow = n, ncol = m)

  # --- reconstruct M ---
  if (!is.null(coefs$u) && !is.null(coefs$d) && !is.null(coefs$v)) {
    if (trace) message("Constructing M ...")
    out$M <- coefs$u %*% (coefs$d * t(coefs$v))
    out$estimates <- out$estimates + out$M
  }

  # --- Reconstructing row covariate effects ---
  if (!is.null(coefs$beta) && !is.null(data$Xq) && !is.null(data$Xr)) {
    if (trace) message("Constructing XBeta ...")

    # Back-transform beta to original scale
    out$beta <- solve(data$Xr, coefs$beta)
    out$xbeta <- data$Xq %*% coefs$beta

    if (meta$shared_effects["beta"]) {
      # Shared beta: xbeta is an n x 1 matrix
      out$estimates <- sweep(out$estimates, 1, as.vector(out$xbeta), "+")
    } else {
      # Unshared beta: xbeta is an n x m matrix.
      out$estimates <- out$estimates + out$xbeta
    }
  }

  # --- Reconstructing column covariate effects ---
  if (!is.null(coefs$gamma) && !is.null(data$Zq) && !is.null(data$Zr)) {
    if (trace) message("Constructing GammaZ ...")

    # Back-transform gamma to original scale
    out$gamma <- coefs$gamma %*% solve(t(data$Zr))
    out$gammaz <- coefs$gamma %*% t(data$Zq)

    if (meta$shared_effects["gamma"]) {
      # Shared gamma: gammaz is a 1 x m matrix.
      out$estimates <- sweep(out$estimates, 2, as.vector(out$gammaz), "+")
    } else {
      # Unshared gamma: gammaz is an n x m matrix.
      out$estimates <- out$estimates + out$gammaz
    }
  }

  # --- Adding Intercepts ---
  if (!is.null(coefs$beta0)) {
    if (trace) message("Constructing row intercepts ...")
    out$estimates <- sweep(out$estimates, 1, as.vector(coefs$beta0), "+")
  }

  if (!is.null(coefs$gamma0)) {
    if (trace) message("Constructing column intercepts ...")
    out$estimates <- sweep(out$estimates, 2, as.vector(coefs$gamma0), "+")
  }

  if (trace) message("done.")
  return(out)
}
#-----------------------------

#' @export
reconstruct_partial <- function(fit, data, target, trace = FALSE) {
  stopifnot(inherits(fit, "imr_fit"))
  stopifnot(inherits(data, "imr_data"))
  stopifnot(IMR::is.Incomplete(target))

  coefs <- fit$coefficients
  meta <- fit$meta

  # --- Reconstruct M ---
  if (!is.null(coefs$u) && !is.null(coefs$d) && !is.null(coefs$v)) {
    if (trace) message("Constructing M ...")
    target@x <- IMR:::partial_crossprod(
      coefs$u,
      sweep(t(coefs$v), 1, coefs$d, "*"),
      target@i, target@p
    )
  } else {
    target@x <- rep(0.0, length(target@i))
  }

  # --- reconstruct row covariate effects ---
  if (!is.null(coefs$beta) && !is.null(data$Xq)) {
    if (trace) message("Constructing XBeta ...")

    if (meta$shared_effects["beta"]) {
      xbeta_vec <- as.numeric(data$Xq %*% coefs$beta)
      IMR:::add_to_rows_inplace_cpp(target@x, target@i, xbeta_vec)
    } else {
      target@x <- target@x + IMR:::partial_crossprod(data$Xq, coefs$beta, target@i, target@p)
    }
  }

  # ---  Column covariate effects ---
  if (!is.null(coefs$gamma) && !is.null(data$Zq)) {
    if (trace) message("Constructing GammaZ ...")

    if (meta$shared_effects["gamma"]) {
      gammaz_vec <- as.numeric(coefs$gamma %*% t(data$Zq))
      IMR:::add_to_cols_inplace_cpp(target@x, target@p, gammaz_vec)
    } else {
      target@x <- target@x + IMR:::partial_crossprod(coefs$gamma, data$Zq, target@i, target@p, vtranspose = TRUE)
    }
  }

  # ---  Intercepts ---
  if (!is.null(coefs$beta0)) {
    if (trace) message("Constructing row intercepts ...")
    IMR:::add_to_rows_inplace_cpp(target@x, target@i, coefs$beta0)
  }

  if (!is.null(coefs$gamma0)) {
    if (trace) message("Constructing column intercepts ...")
    IMR:::add_to_cols_inplace_cpp(target@x, target@p, coefs$gamma0)
  }

  return(target)
}

#----------------------------------------------------------------
#' @export
generate_similarity <- function(x,
                                d = NULL,
                                matern_params = list(smoothness = 1.5, range = 1),
                                rbf_params = list(ell = 1),
                                jitter = 0,
                                invert = FALSE) {
  S <- NULL
  source_type <- "User Matrix"
  params_used <- list()

  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) stop("Input matrix 'x' must be square.")
    S <- x
  } else if (is.character(x)) {
    if (is.null(d)) stop("Distance matrix 'd' is required for kernel generation.")

    if (tolower(x) == "matern") {
      source_type <- "Matern Kernel"
      params_used <- matern_params

      S <- fields::Matern(d,
        smoothness = matern_params$smoothness,
        range = matern_params$range
      )
    } else if (tolower(x) == "rbf") {
      source_type <- "RBF Kernel"
      params_used <- rbf_params
      ell <- rbf_params$ell
      S <- exp(-(d^2) / (2 * ell^2))
    } else {
      stop("Unknown method. 'x' must be a matrix, 'matern', or 'RBF'.")
    }
    if (!invert) {
      warning(paste(
        "Generated a raw Covariance matrix without inversion.",
        "fit_imr() expects the inverse."
      ))
    }
  } else {
    stop("Input 'x' must be a matrix or a character string ('matern', 'RBF').")
  }

  if (invert) {
    S <- tryCatch(
      {
        chol2inv(S)
      },
      error = function(e) {
        stop("Matrix inversion failed (matrix might be singular).")
      }
    )
    source_type <- paste(source_type, "(Inverted)")
  }
  if (is.numeric(jitter) && jitter > 0) {
    S <- S + diag(jitter, nrow(S), ncol(S))
    # source_type <- paste(source_type, "(With Jitter)")
  } else {
    jitter <- 0
  }

  decomp <- eigen(S, symmetric = TRUE)
  evals <- abs(decomp$values)
  cond_num <- max(evals) / min(evals[evals > 0])

  structure(
    list(
      U = decomp$vectors,
      d = decomp$values,
      meta = list(
        source = source_type,
        dim = dim(S),
        params = params_used,
        inverted = invert,
        jitter = jitter,
        condition_number = cond_num
      )
    ),
    class = "imr_similarity"
  )
}

#' @export
print.imr_similarity <- function(x, ...) {
  cat("\n== IMR Similarity Decomposition ==\n")
  cat(sprintf("Source:           %s\n", x$meta$source))
  cat(sprintf("Dimensions:       %d x %d\n", x$meta$dim[1], x$meta$dim[2]))
  cat(sprintf("Jitter value:     %s\n", x$meta$jitter))

  if (length(x$meta$params) > 0) {
    p_names <- names(x$meta$params)
    p_vals <- unlist(x$meta$params)
    cat(sprintf("Parameters:       %s\n", paste(p_names, p_vals, sep = "=", collapse = ", ")))
  }
  cond_fmt <- if (x$meta$condition_number > 1e4) "%.2e" else "%.2f"
  cat(sprintf("Condition Number: %s\n", sprintf(cond_fmt, x$meta$condition_number)))

  cat("Top 5 Eigenvalues:", paste(format(head(x$d, 5), digits = 3), collapse = ", "), "...\n")
  cat("==================================\n")

  invisible(x)
}
