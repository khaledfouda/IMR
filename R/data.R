#' @export

prepare_data <- function(Y, X=NULL, Z=NULL, val_prop = 0.2, seed = 2025){
  out <- list(Xq = NULL, Zq = NULL)
  if((!is.null(seed)) & is.numeric(seed)) set.seed(seed)
  out$Y <- IMR::as.Incomplete(Y)
  message("Performing train/valid split")
  obs_mask <- as.matrix(Y != 0)
  valid_mask <- IMR:::mask_train_test_split(obs_mask, val_prop, seed)
  out$y_train <- as(Y * (1-valid_mask), "Incomplete")
  out$y_valid <- as(Y * (valid_mask), "Incomplete")
  rm(obs_mask)

  if(! is.null(X)){
    stopifnot(is.matrix(X))
    stopifnot(nrow(X) == nrow(Y))
    xqr <- qr(as.matrix(X))
    out$X <- X
    out$Xq = qr.Q(xqr)
    out$Xr = qr.R(xqr)
  }
  if(! is.null(Z)){
    stopifnot(is.matrix(Z))
    stopifnot(nrow(Z) == ncol(Y))
    Zqr <- qr(as.matrix(Z))
    out$Z <- Z
    out$Zq = qr.Q(Zqr)
    out$Zr = qr.R(Zqr)
  }
  return(out)
}



