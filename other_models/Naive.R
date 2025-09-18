naive.fit <- function(Y, X, return_xbeta = FALSE) {
  start_time <- Sys.time()
  svdH <- reduced_hat_decomp.H(X)
  #Y_naive = as.matrix(Y)
  Y_naive = naive_MC(Y)

  Xbeta <-  svdH$u %*% (svdH$v  %*% Y_naive)
  if (return_xbeta)
    return(Xbeta)
  M <- Y_naive - Xbeta
  #----------------------
  # initialization for beta = X^-1 Y
  # comment for later: shouldn't be X^-1 H Y??
  beta = as.matrix(IMR:::inv(X,FALSE) %*% Xbeta)
  obj = list(
    estimates = Y_naive,
    M = M,
    beta = beta,
    time           = round(
      as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs")
      )),
    total_num_fits = 1
  )
  obj$time_per_fit   = obj$time
  return(obj)
}

reduced_hat_decomp.H <- function(X) {
    # returns only the SVD of the hat matrix.
    qrX = Matrix::qr(X)
    rank = qrX$rank
    Q <- qr.Q(qrX)
    H <- Q %*% t(Q)
    svdH <- tryCatch({
      irlba::irlba(H, nu = rank, nv = rank, tol = 1e-5)
    },
    error = function(e) {
      message(paste("SvdH:", e))
      IMR::opt_svd(H, rank)
    })
    list(u = svdH$u,
         v = svdH$d * t(svdH$v),
         rank = rank)
  }
