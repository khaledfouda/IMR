#' @title Identified components of a fitted IMR model
#'
#' @description
#' An optional step after [imr_fit()]. The fitted matrix
#' \eqn{\hat\Theta} is unique, but its split into covariate effects, intercepts
#' and the low-rank component is not. `imr_identify()` returns the one split
#' of the same \eqn{\hat\Theta} that satisfies the orthogonality conditions
#' below. The model is not fitted again.
#'
#' @param fit An `imr_fit` object, from [imr_fit()].
#' @param data The `imr_data` object the model was fitted on.
#'
#' @details
#' Let \eqn{\tilde X = [1_n, X]} when the column intercepts are in the model
#' (\eqn{\tilde X = X} otherwise), and \eqn{\tilde Z = [1_m, Z]} when the row
#' intercepts are in the model. The intercepts are the matching first row of
#' \eqn{\tilde\beta} and first column of \eqn{\tilde\Gamma}. The identified
#' components satisfy
#' \deqn{\tilde X^\top M = 0, \qquad M \tilde Z = 0, \qquad
#'   \tilde X^\top \tilde\Gamma = 0,}
#' and are computed from the fitted ones as
#' \deqn{\tilde\beta \leftarrow \tilde\beta + (\tilde X^\top \tilde X)^{-1}
#'   \tilde X^\top (\tilde\Gamma \tilde Z^\top + A B^\top),}
#' \deqn{\tilde\Gamma \leftarrow P_{\tilde X}^\perp (\tilde\Gamma + A B^\top
#'   \tilde Z (\tilde Z^\top \tilde Z)^{-1}), \qquad
#'   A \leftarrow P_{\tilde X}^\perp A, \qquad B \leftarrow P_{\tilde Z}^\perp B,}
#' where \eqn{M = A B^\top}. Only matrices with one small dimension are formed,
#' and \eqn{\hat\Theta} itself is never constructed.
#'
#' The designs are the ones stored in `data`. When the covariates are centered
#' and the matching intercept is not in the model, the conditions hold for the
#' centered design, which is the model that was fitted.
#'
#' The identified coefficients are in general not sparse, even when the fitted
#' ones are. Models with shared effects are not supported.
#'
#' @return An `imr_fit` object with the identified coefficients, the same
#'   residuals, and `meta$identified = TRUE`. It can be passed to [coef()],
#'   [summary()], [reconstruct()] and [reconstruct_partial()], which give the
#'   same fitted values as `fit`.
#'
#' @seealso [imr_fit()], [reconstruct()]
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' m <- 10
#' Y <- matrix(rnorm(n * m), n, m)
#' Y[sample(n * m, 30)] <- NA
#' X <- matrix(rnorm(n * 2), n, 2)
#'
#' data <- imr_data(Y = Y, X = X, val_prop = 0.0)
#' data <- update(data, col_intercept = TRUE)
#' fit <- imr_fit(data, rank = 2, lambda_beta = 0.1)
#'
#' identified <- imr_identify(fit, data)
#' # the low-rank component is now orthogonal to [1, X]
#' max(abs(crossprod(cbind(1, data$X), identified$coefficients$u)))
#' # and the fitted values are the same
#' max(abs(reconstruct(fit, data, trace = FALSE)$estimates -
#'           reconstruct(identified, data, trace = FALSE)$estimates))
#'
#' @export
imr_identify <- function(fit, data) {
  stopifnot(inherits(fit, "imr_fit"))
  stopifnot(inherits(data, "imr_data"))
  model <- fit$model
  if ((model$row_covariates && model$shared_beta) ||
      (model$col_covariates && model$shared_gamma))
    stop("imr_identify() does not support shared effects.", call. = FALSE)

  coefs <- fit$coefficients
  n <- fit$meta_data$dimensions[1]
  m <- fit$meta_data$dimensions[2]
  has_m <- length(coefs$d) > 0

  # Augmented designs and coefficients ----
  # xt = [1_n, X] with b = [gamma0'; beta], and zt = [1_m, Z] with
  # g = [beta0, gamma]. A side with neither part is NULL. These are the only
  # copies of the designs, and they are released on return.
  xt <- cbind(if (!is.null(coefs$gamma0)) rep(1, n),
              if (!is.null(coefs$beta)) data$X)
  zt <- cbind(if (!is.null(coefs$beta0)) rep(1, m),
              if (!is.null(coefs$gamma)) data$Z)
  b <- rbind(coefs$gamma0, coefs$beta)
  g <- cbind(coefs$beta0, coefs$gamma)
  qx <- if (!is.null(xt)) qr(xt)
  qz <- if (!is.null(zt)) qr(zt)
  if (!is.null(qx) && qx$rank < ncol(xt))
    stop("The row design [1, X] is rank deficient.", call. = FALSE)
  if (!is.null(qz) && qz$rank < ncol(zt))
    stop("The column design [1, Z] is rank deficient.", call. = FALSE)

  # M = A B', with A = U diag(d) and B = V
  if (has_m) {
    a_mat <- t(t(coefs$u) * coefs$d)
    b_mat <- coefs$v
  }

  # Row side: b + (xt' xt)^{-1} xt' (g zt' + A B') ----
  if (!is.null(qx)) {
    if (!is.null(qz)) b <- b + tcrossprod(qr.coef(qx, g), zt)
    if (has_m) b <- b + tcrossprod(qr.coef(qx, a_mat), b_mat)
  }

  # Column side: P_x^perp (g + A B' zt (zt' zt)^{-1}) ----
  if (!is.null(qz)) {
    if (has_m) g <- g + tcrossprod(a_mat, qr.coef(qz, b_mat))
    if (!is.null(qx)) g <- qr.resid(qx, g)
  }

  # Low-rank component: (P_x^perp A) (P_z^perp B)', back to u, d, v ----
  if (has_m) {
    if (!is.null(qx)) a_mat <- qr.resid(qx, a_mat)
    if (!is.null(qz)) b_mat <- qr.resid(qz, b_mat)
    # A = Ua Da Va' and B Va Da = Uw Dw Vw' give A B' = (Ua Vw) Dw Uw'
    sa <- svd(a_mat)
    sw <- svd(t(t(b_mat %*% sa$v) * sa$d))
    r_eff <- max(1, sum(sw$d > max(n, m) * .Machine$double.eps * sw$d[1]))
    keep <- seq_len(r_eff)
    coefs$u <- (sa$u %*% sw$v)[, keep, drop = FALSE]
    coefs$d <- sw$d[keep]
    coefs$v <- sw$u[, keep, drop = FALSE]
  }

  # Split the augmented coefficients back ----
  if (!is.null(coefs$gamma0)) {
    coefs$gamma0 <- b[1, ]
    if (!is.null(coefs$beta)) b <- b[-1, , drop = FALSE]
  }
  if (!is.null(coefs$beta)) {
    dimnames(b) <- dimnames(coefs$beta)
    coefs$beta <- b
  }
  if (!is.null(coefs$beta0)) {
    coefs$beta0 <- g[, 1]
    if (!is.null(coefs$gamma)) g <- g[, -1, drop = FALSE]
  }
  if (!is.null(coefs$gamma)) {
    dimnames(g) <- dimnames(coefs$gamma)
    coefs$gamma <- g
  }

  fit$coefficients <- coefs
  fit$meta$identified <- TRUE
  # the fitted values, hence the residuals, are unchanged; the share of each
  # component is not
  if (!is.null(fit$meta$sum_squares))
    fit$meta$sum_squares <- .imr_sum_squares(
      c(coefs, list(residuals = fit$residuals)), data
    )
  return(fit)
}
