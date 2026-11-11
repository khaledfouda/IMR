#' Prepare Data and Model Structure for IMR
#'
#' @description
#' `imr_data` creates a structured data object for Incomplete Matrix Regression (IMR).
#' It handles target matrix formatting, optional train/validation splitting,
#'  QR decomposition of covariates, and initializes a default model structure
#'  based on the provided inputs.
#'
#' @param Y A numeric matrix (dense or sparse) representing the target/response
#' matrix. This is the only required parameter.
#' Missing values and Zeros are treated as unobserved.
#' @param X A numeric matrix of row covariates. One of the dimensions must match the number of rows of `Y`.
#' The other dimension must be strictly smaller. Defaults to `NULL`.
#' @param Z A numeric matrix of column covariates. One of the dimensions must match the number of columns of `Y`
#'   The other dimension must be strictly smaller. Defaults to `NULL`.
#' @param similarity_rows An `"imr_similarity"` object (see \code{\link{imr_similarity}})
#'   describing the row similarity/information matrix. Defaults to `NULL`.
#' @param similarity_cols An `"imr_similarity"` object (see \code{\link{imr_similarity}})
#'   describing the column similarity/information matrix. Defaults to `NULL`.
#' @param val_prop Numeric. The proportion of observed entries in `Y` to hold out for
#'   validation. If `0.0` (the default), no validation split is performed. Validation is required
#'   for hyperparameter tuning.
#' @param seed Integer. An optional random seed for reproducibility when creating the
#'   train/validation split. Defaults to `NULL`.
#' @param train_mask,valid_mask Optional matrices of 0s and 1s with the
#'   dimensions of `Y`, marking the training and the validation entries. When
#'   `valid_mask` is given, `train_mask` is required as well, the two may not
#'   overlap, and both may mark observed entries only. The split is then taken
#'   from the masks, and `val_prop` and `seed` are not used. Both default to
#'   `NULL`.
#' @param center,scale Logical. If `TRUE` (the default), each column of `X` and
#'   `Z` is centered to mean zero (`center`) and scaled to unit standard
#'   deviation (`scale`). One transformed copy of each matrix is stored, with its
#'   column means and standard deviations. The penalties apply to the
#'   coefficients of the transformed covariates, and `coef()`, `summary()` and
#'   `reconstruct()` report the coefficients on the original scale. A constant
#'   column cannot be scaled and raises an error. Without the matching intercept
#'   in the model, centering changes the fitted model.
#'
#' @details
#' By default, `imr_data` assumes that if you provide `X`, `Z`, or similarity matrices,
#' you want to include them in your IMR model. It sets `$model` logical
#' flags to `TRUE` for any provided data. If you wish to change this model structure
#' later without re-processing the data matrices, use the `update()` method.
#'
#' @return An object of class `"imr_data"`, containing:
#' \itemize{
#'   \item \code{Y}, \code{y_train}, \code{y_valid}: The target matrix and, when
#'   \code{val_prop > 0}, its training/validation splits (as `Incomplete` objects).
#'   \item \code{Xq}, \code{Xr}, \code{Zq}, \code{Zr}: The QR decompositions of
#'   the covariates.
#'   \item \code{similarity_rows}, \code{similarity_cols}: The provided similarity
#'   matrices.
#'   \item \code{meta}: A list of dataset dimensions, sparsity, and variable counts.
#'   \item \code{model}: A list of logical flags defining the current model structure.
#' }
#'
#' @seealso [update.imr_data()], [imr_similarity()]
#' @export
#'
#' @examples
#'
#' # create sample data
#' Y <- matrix(
#'   c(2, NA, 3,
#'   4, .5, NA,
#'   NA, 1, 5), 3, byrow= TRUE
#' )
#'
#' # create covariate matrix of 2 variables
#' X <- matrix(rnorm(3*2), 3, 2)
#'
#' # create the data object
#' data <- imr_data(Y =  Y, X = X, val_prop = 0.0)
#'
#' # update the model to add row intercepts
#' data <- update(data, row_intercept = TRUE)
#'
#' # print the metadata
#' print(data)
#'
imr_data <- function(Y,
                     X = NULL,
                     Z = NULL,
                     similarity_rows = NULL,
                     similarity_cols = NULL,
                     val_prop = 0.0,
                     seed = NULL,
                     train_mask = NULL,
                     valid_mask = NULL,
                     center = TRUE,
                     scale = TRUE) {
  .imr_check_param(center, "bool")
  .imr_check_param(scale, "bool")
  out <- list()
  if (!is.null(seed) && is.numeric(seed)) set.seed(seed)

  # ---  Target Matrix Setup ---
  out$Y <- Y <- as_incomplete(Y)
  out$Y@x <- out$Y@x + 0 # this is to force a copy
  obs_mask <- as.matrix(Y != 0)

  dims <- dim(Y)
  n_total_obs <- sum(obs_mask)

  # every row and column needs at least two observed entries
  n_few_rows <- sum(rowSums(obs_mask) < 2)
  n_few_cols <- sum(colSums(obs_mask) < 2)
  if (n_few_rows > 0 || n_few_cols > 0) {
    stop(sprintf(
      "Every row and column of Y needs at least two observed entries; %d rows and %d columns have fewer.",
      n_few_rows, n_few_cols
    ), call. = FALSE)
  }

  # --- Train/Validation Split ---
  # if a validation mask is provided, it'll be used!
  split_data <- is.numeric(val_prop) && val_prop > 0
  is_binary_matrix <- function(m) {
    is.matrix(m) && isTRUE(all(m == 0 | m == 1))
  }

  if(!is.null(valid_mask)){
    if (!is_binary_matrix(valid_mask)) {
      stop("valid_mask must be a matrix of 0s and 1s or NULL")
    }
    if (!is_binary_matrix(train_mask)) {
      stop("If valid_mask provided then train_mask must be a matrix of 0s and 1s.")
    }
    if (any(train_mask == 1 & valid_mask == 1)) {
      stop("train_mask and valid_mask cannot overlap.")
    }
    if (any((train_mask == 1 | valid_mask == 1) & obs_mask != 1)) {
      stop("train_mask or valid_mask has a 1 where obs_mask does not.")
    }

    out$valid_mask <- as_incomplete(valid_mask)
    out$y_train <- as_incomplete(Y * train_mask)
    out$y_valid <- as_incomplete(Y * valid_mask)

    n_train <- sum(train_mask)
    n_valid <- sum(valid_mask)

    split_data <- TRUE


  } else if (split_data) {
    # message("Performing train/valid split...")

    valid_mask_mat <- mask_train_test_split(obs_mask, val_prop, seed)
    train_mask_mat <- obs_mask * (1 - valid_mask_mat)

    out$valid_mask <- as_incomplete(valid_mask_mat)
    out$y_train <- as_incomplete(Y * train_mask_mat)
    out$y_valid <- as_incomplete(Y * valid_mask_mat)

    n_train <- sum(train_mask_mat)
    n_valid <- sum(valid_mask_mat)
  }else {
    n_train <- n_total_obs
    n_valid <- 0
  }

  # and the training split must keep at least one entry in every row and column
  if (split_data) {
    n_empty_rows <- sum(tabulate(out$y_train@i + 1, dims[1]) == 0)
    n_empty_cols <- sum(diff(out$y_train@p) == 0)
    if (n_empty_rows > 0 || n_empty_cols > 0) {
      stop(sprintf(
        "The training split leaves %d rows and %d columns of Y without observed entries; use another seed, a smaller val_prop or other masks.",
        n_empty_rows, n_empty_cols
      ), call. = FALSE)
    }
  }

  out$obs_mask <- as_incomplete(obs_mask * 1)

  # ---  Similarity Matrices ---
  stopifnot(is.null(similarity_rows) || inherits(similarity_rows, "imr_similarity"))
  stopifnot(is.null(similarity_cols) || inherits(similarity_cols, "imr_similarity"))
  out$similarity_rows <- similarity_rows
  out$similarity_cols <- similarity_cols

  # ---  Covariates ---
  # One copy of each matrix is kept. update.imr_data() below centers and scales it
  # in place and computes its step size (data$LX, data$LZ).
  if (!is.null(X)) {
    stopifnot(is.matrix(X))
    if(nrow(X) != nrow(Y))
      if(ncol(X) == nrow(Y)) X = t(X) else stop("None of the dimensions of X match the number of rows of Y")
    out$X <- X
  }

  if (!is.null(Z)) {
    stopifnot(is.matrix(Z))
    if(nrow(Z) != ncol(Y))
      if(ncol(Z) == ncol(Y)) Z = t(Z) else stop("None of the dimensions of Z match the number of columns of Y")
    out$Z <- Z
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
    names_X_vars = if (!is.null(X)) colnames(X) else NULL,
    has_Z = !is.null(Z),
    num_Z_vars = if (!is.null(Z)) ncol(Z) else 0,
    names_Z_vars = if (!is.null(Z)) colnames(Z) else NULL,
    has_sim_row = !is.null(out$similarity_rows),
    has_sim_col = !is.null(out$similarity_cols)
  )
  # --- Model Structure ----
  out$model <- list(
    row_covariates = out$meta$has_X,
    col_covariates = out$meta$has_Z,
    low_rank_component = TRUE,
    row_similarity = out$meta$has_sim_row,
    col_similarity = out$meta$has_sim_col,
    row_intercept = FALSE,
    col_intercept = FALSE,
    shared_beta = FALSE,
    shared_gamma = FALSE
  )

  # centers and scales X and Z, and computes their step sizes
  update.imr_data(structure(out, class = "imr_data"), center = center, scale = scale)
}
#----------------------------------------
#' Update IMR Model Structure
#'
#' @description
#' Updates the logical model structure of an existing `"imr_data"` object. This allows
#' you to toggle model components (like intercepts, covariates, or low-rank components)
#' on or off without needing to re-process the underlying data matrices.
#'
#' @param object An object of class `"imr_data"`.
#' @param row_covariates Logical. Include row covariates in the model?
#'   (Requires `X` to have been provided to `imr_data`).
#' @param col_covariates Logical. Include column covariates in the model?
#'   (Requires `Z` to have been provided to `imr_data`).
#' @param low_rank_component Logical. Include the low-rank component (\eqn{M=U D V^T})?
#' @param row_similarity Logical. Include row similarity penalties?
#'   (Requires `similarity_rows` to have been provided to `imr_data`).
#' @param col_similarity Logical. Include column similarity penalties?
#'   (Requires `similarity_cols` to have been provided to `imr_data`).
#' @param row_intercept Logical. Include a row intercept/bias term?
#' @param col_intercept Logical. Include a column intercept/bias term?
#' @param shared_beta Logical. Should row covariate effects (`beta`) be shared
#'   across all columns?
#' @param shared_gamma Logical. Should column covariate effects (`gamma`) be shared
#'   across all rows?
#' @param center,scale Logical or `NULL`. `NULL` (the default) keeps the current
#'   centering or scaling of `X` and `Z`. `TRUE` or `FALSE` applies or removes it
#'   on the stored copy of each matrix, and the step sizes are recomputed. See
#'   [imr_data()].
#' @param ... Additional arguments (ignored).
#'
#' @return The modified `"imr_data"` object with updated `$model` flags.
#'
#' @seealso [imr_data()], [imr_similarity()]
#'
#' @inherit imr_data examples
#'
#' @export
#' @method update imr_data
update.imr_data <- function(object,
                            row_covariates = NULL,
                            col_covariates = NULL,
                            low_rank_component = NULL,
                            row_similarity = NULL,
                            col_similarity = NULL,
                            row_intercept = NULL,
                            col_intercept = NULL,
                            shared_beta = NULL,
                            shared_gamma = NULL,
                            center = NULL,
                            scale = NULL,
                            ...) {
  update_flag <- function(obj, flag_val, flag_name, dependency_name) {
    if (!is.null(flag_val)) {
      # Ensure it's a single boolean
      flag_val <- as.logical(flag_val[1])

      if (flag_val && !obj$meta[[dependency_name]]) {
        stop(sprintf(
          "Cannot set '%s = TRUE' because the underlying data ('%s') was not provided.",
          flag_name, dependency_name
        ), call. = FALSE)
      }
      obj$model[[flag_name]] <- flag_val
    }
    return(obj)
  }

  object <- update_flag(object, row_covariates, "row_covariates", "has_X")
  object <- update_flag(object, col_covariates, "col_covariates", "has_Z")
  object <- update_flag(object, row_similarity, "row_similarity", "has_sim_row")
  object <- update_flag(object, col_similarity, "col_similarity", "has_sim_col")

  toggle <- function(obj, val, name) {
    if (!is.null(val)) obj$model[[name]] <- as.logical(val[1])
    return(obj)
  }

  object <- toggle(object, low_rank_component, "low_rank_component")
  object <- toggle(object, row_intercept, "row_intercept")
  object <- toggle(object, col_intercept, "col_intercept")

  if (!is.null(shared_beta)) {
    if (shared_beta && !object$model$row_covariates) {
      warning("Setting shared_beta = TRUE does nothing if row covariates are inactive.")
    }
    object$model$shared_beta <- as.logical(shared_beta[1])
  }

  if (!is.null(shared_gamma)) {
    if (shared_gamma && !object$model$col_covariates) {
      warning("Setting shared_gamma = TRUE does nothing if col covariates are inactive.")
    }
    object$model$shared_gamma <- as.logical(shared_gamma[1])
  }

  # --- Centering and scaling of X and Z ---
  # One copy of each matrix is stored. With its column centers c and scales s
  # (data$x_center, data$x_scale; NULL when not applied), the original matrix is
  # stored * s + c, by column. A change is applied to the stored copy, and only
  # then are the moments and the step size 1 / ||X||_op^2 computed.
  if (!is.null(center)) center <- as.logical(center[1])
  if (!is.null(scale)) scale <- as.logical(scale[1])
  rescale <- function(obj, name) {
    key <- tolower(name)
    c_old <- obj[[paste0(key, "_center")]]
    s_old <- obj[[paste0(key, "_scale")]]
    to_center <- if (is.null(center)) !is.null(c_old) else center
    to_scale <- if (is.null(scale)) !is.null(s_old) else scale
    step <- paste0("L", name)
    if (to_center == !is.null(c_old) && to_scale == !is.null(s_old) && !is.null(obj[[step]]))
      return(obj)

    mat <- obj[[name]]
    s <- if (is.null(s_old)) 1 else s_old
    shift <- if (is.null(c_old)) 0 else c_old
    # centers and scales of the original matrix
    c_new <- if (!to_center) NULL else if (!is.null(c_old)) c_old else s * colMeans(mat)
    s_new <- if (!to_scale) NULL else if (!is.null(s_old)) s_old else {
      constant <- which(apply(mat, 2, function(x) all(x == x[1])))
      if (length(constant) > 0)
        stop(sprintf("Constant columns of %s cannot be scaled: %s.", name, toString(constant)),
             call. = FALSE)
      apply(mat, 2, sd)
    }
    # new stored copy: (mat * s + shift - c_new) / s_new
    divisor <- if (is.null(s_new)) 1 else s_new
    a <- s / divisor
    b <- (shift - if (is.null(c_new)) 0 else c_new) / divisor
    if (any(a != 1)) mat <- sweep(mat, 2, a, "*")
    if (any(b != 0)) mat <- sweep(mat, 2, b, "+")

    obj[[name]] <- mat
    obj[[paste0(key, "_center")]] <- c_new
    obj[[paste0(key, "_scale")]] <- s_new
    obj[[step]] <- 1 / max(eigen(crossprod(mat), symmetric = TRUE, only.values = TRUE)$values)
    obj
  }
  if (object$meta$has_X) object <- rescale(object, "X")
  if (object$meta$has_Z) object <- rescale(object, "Z")

  return(object)
}

#' @title Print an IMR data object
#' @param x An `imr_data` object
#' @param ... Additional arguments to comply with generic function
#' @return The input `imr_data` object `x`, invisibly. Called for its side effect
#'   of printing the data dimensions, split, and model configuration to the console.
#' @seealso [imr_data()] [update.imr_data()] [imr_similarity()]
#' @inherit imr_data examples
#' @export
#' @method print imr_data
print.imr_data <- function(x, ...) {
  m <- x$meta # What data is available in memory
  a <- x$model # What the solver is instructed to use

  cat("\n== IMR Data Object ==\n")

  # --- 1. Base Dimensions ---
  n_cells <- m$dimensions[1] * m$dimensions[2]
  n_miss <- n_cells - m$total_obs
  cat(sprintf("Target Matrix (Y)  : %d x %d  (%d cells)\n", m$dimensions[1], m$dimensions[2], n_cells))
  cat(sprintf("Observed           : %d (%.2f%%)\n", m$total_obs, 100 * m$total_obs / n_cells))
  cat(sprintf("Missing (sparsity) : %d (%.2f%%)\n", n_miss, 100 * m$sparsity))

  # --- 2. Train/Valid Split info (percentages of the observed entries) ---
  cat("\n-- Data Split --\n")
  if (m$split_data) {
    cat(sprintf("Training           : %d (%.1f%% of observed)\n", m$n_train, 100 * m$n_train / m$total_obs))
    cat(sprintf("Validation         : %d (%.1f%% of observed)\n", m$n_valid, 100 * m$n_valid / m$total_obs))
  } else {
    cat("Using 100% of observed data (no validation split)\n")
  }

  # --- 3. Active Submodel Configuration ---
  cat("\n-- Model Configuration --\n")

  # One consistent line per component: [active] / [inactive] / [not provided],
  # optionally followed by a short descriptor.
  cfg <- function(label, provided, active, extra = "") {
    tag <- if (!provided) "[not provided]" else if (active) "[active]" else "[inactive]"
    line <- if (nzchar(extra)) paste0(tag, "  ", extra) else tag
    cat(sprintf("  %-22s : %s\n", label, line))
  }

  cov_extra <- function(provided, active, nvars, shared) {
    if (!provided) return("")
    s <- sprintf("%d var%s", nvars, if (nvars == 1) "" else "s")
    if (active) s <- paste0(s, if (shared) ", shared" else ", per-column")
    s
  }

  # Low-rank is a purely algorithmic component (doesn't depend on external data)
  cfg("Low-Rank Matrix (M)", TRUE, a$low_rank_component)
  cfg("Row Intercepts", TRUE, a$row_intercept)
  cfg("Column Intercepts", TRUE, a$col_intercept)
  cfg("Row Covariates (X)", m$has_X, a$row_covariates,
    cov_extra(m$has_X, a$row_covariates, m$num_X_vars, a$shared_beta))
  cfg("Column Covariates (Z)", m$has_Z, a$col_covariates,
    cov_extra(m$has_Z, a$col_covariates, m$num_Z_vars, a$shared_gamma))
  cfg("Row Similarity", m$has_sim_row, a$row_similarity)
  cfg("Column Similarity", m$has_sim_col, a$col_similarity)

  cat("==========================\n")
  invisible(x)
}


#' @title Reconstruct Model Estimates and Components
#'
#' @description
#' Calculates the estimated response matrix from a fitted IMR model. In addition
#' to the final predicted values, it computes and returns the individual structural
#' components of the model, such as the low-rank matrix, covariate effects, and intercepts.
#'
#' @param fit An object of class `"imr_fit"`, generated by \code{\link{imr_fit}}.
#' @param data An object of class `"imr_data"`, generated by \code{\link{imr_data}}.
#' @param trace Logical. If `TRUE` (the default), prints progress messages for
#'   each reconstruction step.
#'
#' @details
#' During the model fitting process, covariates are made orthogonal using
#' QR decomposition. The `reconstruct`
#' function reverses this process, using the `R` matrices (`Xr` and `Zr`) to
#' back-transform the estimated coefficients (`beta` and `gamma`) to the original
#' scale of the input data.
#'
#' The final `estimates` matrix is computed by summing the active components
#' based on the model structure: the low-rank matrix (\eqn{M = U D V^T}), row covariate
#' effects (\eqn{X \beta}), column covariate effects (\eqn{\Gamma Z}), and any intercepts.
#' The function automatically handles matrix dimension sweeping for shared effects
#'
#' @return A list containing the reconstructed components and final estimates:
#' \itemize{
#'   \item \code{beta}: The row covariate coefficients transformed back to the original scale.
#'   \item \code{gamma}: The column covariate coefficients transformed back to the original scale.
#'   \item \code{M}: The reconstructed low-rank component matrix.
#'   \item \code{xbeta}: The matrix of row covariate effects.
#'   \item \code{gammaz}: The matrix of column covariate effects.
#'   \item \code{beta0}: The estimated row intercepts.
#'   \item \code{gamma0}: The estimated column intercepts.
#'   \item \code{estimates}: The final predicted response matrix (the sum of all active components).
#' }
#'
#' @seealso [reconstruct_partial()]
#'
#' @inherit imr_fit examples
#'
#' @export
reconstruct <- function(fit, data, trace = TRUE) {
  stopifnot(inherits(fit, "imr_fit"))
  stopifnot(inherits(data, "imr_data"))

  out <- list(
    beta = NULL, gamma = NULL, M = NULL,
    beta0 = NULL, gamma0 = NULL,
    xbeta = NULL, gammaz = NULL, estimates = 0
  )

  # coefficients on the original scale of the covariates (see coef.imr_fit()), and on
  # the scale of the stored designs data$X and data$Z
  coefs <- coef(fit)
  fitted <- fit$coefficients
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
  if (!is.null(coefs$beta) && !is.null(data$X)) {
    if (trace) message("Constructing XBeta ...")

    out$beta <- coefs$beta
    # the stored design gives (X - 1 xbar') beta when X is centered; the shift
    # beta' xbar is added back when gamma0 absorbs it (see coef.imr_fit())
    out$xbeta <- data$X %*% fitted$beta
    if (!is.null(fit$scaling$x_center) && !is.null(coefs$gamma0))
      out$xbeta <- sweep(out$xbeta, 2, as.vector(crossprod(coefs$beta, fit$scaling$x_center)), "+")

    if (meta$shared_effects["beta"]) {
      # Shared beta: xbeta is an n x 1 matrix
      out$estimates <- sweep(out$estimates, 1, as.vector(out$xbeta), "+")
    } else {
      # Unshared beta: xbeta is an n x m matrix.
      out$estimates <- out$estimates + out$xbeta
    }
  }

  # --- Reconstructing column covariate effects ---
  if (!is.null(coefs$gamma) && !is.null(data$Z)) {
    if (trace) message("Constructing GammaZ ...")

    out$gamma <- coefs$gamma
    # as above, with the shift Gamma zbar absorbed by beta0
    out$gammaz <- fitted$gamma %*% t(data$Z)
    if (!is.null(fit$scaling$z_center) && !is.null(coefs$beta0))
      out$gammaz <- sweep(out$gammaz, 1, as.vector(coefs$gamma %*% fit$scaling$z_center), "+")

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
    out$beta0 <- coefs$beta0
    out$estimates <- sweep(out$estimates, 1, as.vector(coefs$beta0), "+")
  }

  if (!is.null(coefs$gamma0)) {
    if (trace) message("Constructing column intercepts ...")
    out$gamma0 <- coefs$gamma0
    out$estimates <- sweep(out$estimates, 2, as.vector(coefs$gamma0), "+")
  }

  if (trace) message("done.")
  return(out)
}
#-----------------------------

#' Partially Reconstruct Model Estimates for Specific Indices
#'
#' @description
#' Computes the estimated response matrix from a fitted IMR model, but *only* for
#' a specific subset of matrix entries. This is highly memory-efficient and fast
#' when you only need to evaluate predictions on a test set or validation mask,
#' as it avoids computing the full dense prediction matrix.
#'
#' @param fit An object of class `"imr_fit"`, generated by \code{\link{imr_fit}}.
#' @param data An object of class `"imr_data"`, generated by \code{\link{imr_data}}.
#' @param irow,pcol Integer vectors corresponding to the `@i` and `@p` slots of a
#'    `"CsparseMatrix"` (`dgCMatrix`), denoting the indices of the elements to be
#'    estimated. The returned values correspond to the `@x` slot of that object.
#' @param trace Logical. If `TRUE`, prints progress messages for each reconstruction
#'   step. Defaults to `FALSE`.
#' @param return_matrix Logical. If `TRUE` then it returns a `CsparseMatrix` with
#'    the following slots: `@i=irow`, `@p=pcol`, `@x=estimates`. If `FALSE`, only
#'    the `@x` slot is returned as a vector. Defaults to `FALSE`.
#'
#' @details
#' Unlike the full \code{\link{reconstruct}} function, which allocates and calculates
#' every cell of an \eqn{n \times m} matrix, `reconstruct_partial` relies on optimized C++
#' routines to calculate only the targeted entries.
#'
#' @return Either a `CsparseMatrix` object or a numeric vector (see above).
#'
#' @seealso [reconstruct()]
#' @importFrom methods new
#' @inherit imr_fit examples
#' @export
reconstruct_partial <- function(fit, data, irow, pcol, trace = FALSE, return_matrix = FALSE) {
  stopifnot(inherits(fit, "imr_fit"))
  stopifnot(inherits(data, "imr_data"))

  coefs <- fit$coefficients
  meta <- fit$meta

  # --- Reconstruct M ---
  if (!is.null(coefs$u) && !is.null(coefs$d) && !is.null(coefs$v)) {
    if (trace) message("Constructing M ...")
    target <- partial_crossprod(
      coefs$u,
      sweep(t(coefs$v), 1, coefs$d, "*"),
      irow, pcol
    )
  } else {
    target <- rep(0.0, length(irow))
  }

  # --- reconstruct row covariate effects ---
  if (!is.null(coefs$beta) && !is.null(data$X)) {
    if (trace) message("Constructing XBeta ...")

    if (meta$shared_effects["beta"]) {
      xbeta_vec <- as.numeric(data$X %*% coefs$beta)
      # xbeta_vec <- as.numeric(data$Xq %*% coefs$beta)
      add_to_rows_inplace_cpp(target, irow, xbeta_vec)
    } else {
      target <- target + partial_crossprod(data$X, coefs$beta, irow, pcol)
      #target <- target + partial_crossprod(data$Xq, coefs$beta, irow, pcol)
    }
  }

  # ---  Column covariate effects ---
  if (!is.null(coefs$gamma) && !is.null(data$Z)) {
    if (trace) message("Constructing GammaZ ...")

    if (meta$shared_effects["gamma"]) {
      gammaz_vec <- as.numeric(coefs$gamma %*% t(data$Z))
      # gammaz_vec <- as.numeric(coefs$gamma %*% t(data$Zq))
      add_to_cols_inplace_cpp(target, pcol, gammaz_vec)
    } else {
      target <- target + partial_crossprod(coefs$gamma, data$Z, irow, pcol, vtranspose = TRUE)
      # target <- target + partial_crossprod(coefs$gamma, data$Zq, irow, pcol, vtranspose = TRUE)
    }
  }

  # ---  Intercepts ---
  if (!is.null(coefs$beta0)) {
    if (trace) message("Constructing row intercepts ...")
    add_to_rows_inplace_cpp(target, irow, coefs$beta0)
  }

  if (!is.null(coefs$gamma0)) {
    if (trace) message("Constructing column intercepts ...")
    add_to_cols_inplace_cpp(target, pcol, coefs$gamma0)
  }

  if(return_matrix)
    target <- as_incomplete(new("dgCMatrix",
                  i = irow,
                  p = pcol,
                  x = target,
                  Dim = fit$meta_data$dimensions))

  return(target)
}

#----------------------------------------------------------------

#' @title Generate Similarity or Information Matrix Object
#'
#' @description
#' Creates an `"imr_similarity"` object containing the Eigenvalue-decomposition of a
#' similarity or information matrix. This object is required for incorporating
#' row or column similarities in an IMR model. The function can accept a
#' matrix or generate a covariance matrix from distance data using
#' specified kernels.
#'
#' @param x Either a square numeric matrix, or a character string specifying a
#'   kernel method (`"matern"` or `"rbf"`). If a matrix is provided, it should
#'   be a similarity or information matrix (if `invert = FALSE`), or a covariance
#'   matrix (if `invert = TRUE`).
#' @param normalize Logical. Should the kernel be symmetrically normalized by its
#'   degree matrix (\eqn{S \leftarrow D^{-1/2} S D^{-1/2}}, where \eqn{D} is the
#'   diagonal matrix of row sums)? This assumes a non-negative kernel and is highly
#'   recommended for better performance. See the appendix of the accompanying paper
#'   (Fouda et al., 2026; \url{https://arxiv.org/abs/2606.26325}) for details.
#'   Defaults to `TRUE`.
#' @param invert Logical. Should the kernel be inverted? The IMR model
#'   expects an information matrix. If your input (or generated kernel)
#'   represents a covariance matrix, you must set this to `TRUE` to invert it.
#'   Defaults to `TRUE`.
#' @param jitter Numeric. A small positive value added to the diagonal of the
#'   matrix to improve numerical stability and reduce the condition number.
#'   This is applied before any inversion. A value between `0` and `1` and appropriately away from
#'   the boundaries is recommended. Defaults to `0.2`.
#' @param d A numeric distance matrix. This is strictly required if `x` is set to
#'   `"matern"` or `"rbf"`. Defaults to `NULL`.
#' @param matern_smoothness,matern_range  the parameters for the Matern kernel:
#'   `smoothness` and `range`. Defaults to `(smoothness = 1.5, range = 1)`.
#' @param rbf_ell The length-scale parameter for the Radial Basis Function
#'   (RBF) kernel. Defaults to `1`.
#'
#' @return An object of class `"imr_similarity"`. This is a list containing:
#' \itemize{
#'   \item \code{U}: A matrix of eigenvectors.
#'   \item \code{d}: A vector of  eigenvalues.
#'   \item \code{meta}: A list of metadata.
#' }
#' @seealso [imr_data()], [print.imr_similarity()]
#'
#' @examples
#' # generate 5 random spatial locations
#' coords <- data.frame(
#' x = runif(5, 0, 10),
#' y = runif(5, 0 , 10))
#'
#' # compute the distance matrix
#' distance_matrix <- as.matrix(dist(coords))
#'
#' # generate the similarity object of a 5/2 matern kernel
#' sim <- imr_similarity(x = "matern", d = distance_matrix, matern_smoothness = 1.5)
#'
#' # print the matrix's metadata
#' print(sim)
#'
#' @references
#' Fouda, K., Labbe, A., & Oualkacha, K. (2026).
#' \emph{Incomplete Matrix Regression}.
#' \url{https://arxiv.org/abs/2606.26325}
#'
#' @export
imr_similarity <- function(x,
                           normalize = TRUE,
                           invert = TRUE,
                           jitter = 0.2,
                           d = NULL,
                           matern_smoothness = 1.5,
                           matern_range = 1,
                           rbf_ell = 1) {
  S <- NULL
  source_type <- "User Matrix"
  params_used <- list()

  stopifnot(is.numeric(jitter) && jitter >= 0 && jitter < 1)

  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) stop("Input matrix 'x' must be square.")
    S <- x
  } else if (is.character(x)) {
    if (is.null(d)) stop("Distance matrix 'd' is required for kernel generation.")

    if (tolower(x) == "matern") {
      source_type <- "Matern Kernel"
      params_used <- c(smoothness = matern_smoothness,
                       range      = matern_range)

      S <- fields::Matern(d,
        smoothness = matern_smoothness,
        range = matern_range
      )
    } else if (tolower(x) == "rbf") {
      source_type <- "RBF Kernel"
      params_used <- c(ell = rbf_ell)
      ell <- rbf_ell
      S <- exp(-(d^2) / (2 * ell^2))
    } else {
      stop("Unknown method. 'x' must be a matrix, 'matern', or 'RBF'.")
    }
    # if (!invert) {
    #   warning(paste(
    #     "Generated a raw Covariance matrix without inversion.",
    #     "fit_imr() expects the inverse."
    #   ))
    # }
  } else {
    stop("Input 'x' must be a matrix or a character string ('matern', 'RBF').")
  }

  if(!(normalize && invert && jitter > 0))
    warning(paste("For best practice, normalize the kernel and add a small jitter between 0 and 1",
    "but not too close to the boundaries. fit_imr() expects the inverse of the kernel so make sure",
    "to either set invert=TRUE or provide the kernel inverted. If the kernel is provided inverted,
    you should set normalize=FALSE and jitter=FALSE."))

  if(normalize) {
    #if(!invert)
     # warning("Normalizing without inverting. Normalization should be done to the Kernel.")}
    sqrt_inv_degree <- 1 / sqrt(rowSums(S))
    # S = D^{-1/2} S D^{-1/2}
    S <- sweep(sweep(S, 1, sqrt_inv_degree, "*"), 2, sqrt_inv_degree, "*")
    source_type <- paste(source_type, "(Normalized)")
  }

  if (jitter > 0) {
    #if(!invert)
      #warning("Additing Jitter without inverting. Jitter should be added to the Kernel.")
    S <- S + diag(jitter, nrow(S), ncol(S))
  } else {
    jitter <- 0
  }

  if (invert) {
    S <- inv(S)
    source_type <- paste(source_type, "(Inverted)")
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


#' @title Print an IMR similarity object
#' @param x An `imr_similarity` object
#' @param ... Additional arguments to comply with generic function
#' @return The input `imr_similarity` object `x`, invisibly. Called for its side
#'   effect of printing the decomposition metadata to the console.
#' @seealso [imr_similarity()]
#' @inherit imr_similarity examples
#' @importFrom utils head
#' @export
#' @method print imr_similarity
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

  k <- min(5, length(x$d))
  ellipsis <- if (length(x$d) > k) ", ..." else ""
  cat(sprintf(
    "Top %d Eigenvalues: %s%s\n", k,
    paste(format(head(x$d, k), digits = 3), collapse = ", "), ellipsis
  ))
  cat("==================================\n")

  invisible(x)
}
