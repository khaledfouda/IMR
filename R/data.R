#' @export

prepare_data <- function(Y, X = NULL, Z = NULL, val_prop = 0.2, seed = 2025) {
  out <- list( model = list(Xq = NULL, Zq = NULL))
  if ((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  out$model$Y <- out$Y <- IMR::as.Incomplete(Y)
  message("Performing train/valid split")
  obs_mask <- as.matrix(Y != 0)
  out$model$valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
  out$model$train_mask <- IMR::as.Incomplete(obs_mask * (1-out$model$valid_mask))
  out$model$valid_mask <- IMR::as.Incomplete(out$model$valid_mask)
  out$model$y_train <- IMR::as.Incomplete(Y * out$model$train_mask)
  out$model$y_valid <- IMR::as.Incomplete(Y * out$model$valid_mask)
  rm(obs_mask)

  #out$y_train <- as(Y * (1 - valid_mask), "Incomplete")
  #out$y_valid <- as(Y * (valid_mask), "Incomplete")
  #out$valid_mask <- valid_mask

  if (!is.null(X)) {
    stopifnot(is.matrix(X))
    stopifnot(nrow(X) == nrow(Y))
    xqr <- qr(as.matrix(X))
    out$X <- X
    out$model$Xq <- qr.Q(xqr)
    out$Xr <- qr.R(xqr)
  }
  if (!is.null(Z)) {
    stopifnot(is.matrix(Z))
    stopifnot(nrow(Z) == ncol(Y))
    Zqr <- qr(as.matrix(Z))
    out$Z <- Z
    out$model$Zq <- qr.Q(Zqr)
    out$Zr <- qr.R(Zqr)
  }

  return(out)
}
