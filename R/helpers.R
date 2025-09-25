#' @export
inv <- function(X, is_square = nrow(X) == ncol(X)) {
  if (is_square) {
    # Try to apply solve() and catch any errors indicating singularity
    tryCatch(
      {
        return(solve(X))
      },
      error = function(e) {
        return(MASS::ginv(X))
      }
    )
  } else {
    # Use ginv() for non-square matrices
    inv_X <-
      return(MASS::ginv(X))
  }
}
#-------------------------------------
#' @export
reconstruct <- function(fit, dat, partial=FALSE, trace=TRUE){
  out <- list(beta=NA, gamma=NA, M = NA, xbeta=NA, gammaz=NA, estimates=0)
  check_mat <- function(mat, is_matrix=TRUE){
    if(any(is.na(mat))) return(FALSE)
    if(is_matrix & is.matrix(mat)) return(TRUE)
    if(is_matrix) return(FALSE)
    if(is.vector(mat)) return(TRUE)
    return(FALSE)
  }
  #-----
  if(check_mat(fit$u) & check_mat(fit$d,FALSE) & check_mat(fit$v)){
    if(trace) message("Constructing M ...")
    out$M <- fit$u %*% (fit$d * t(fit$v))
    out$estimates <- out$M
  }
  if(check_mat(fit$beta) & check_mat(dat$X) & check_mat(dat$Xr)){
    if(trace) message("Constructing XBeta ...")
    out$beta <- solve(dat$Xr) %*% fit$beta
    out$xbeta <- dat$X %*% out$beta
    out$estimates <- out$estimates + out$xbeta
  }
  if(check_mat(fit$gamma) & check_mat(dat$Z) & check_mat(dat$Zr)){
    if(trace) message("Constructing GammaZ ...")
    out$gamma <- fit$gamma %*% solve(t(dat$Zr))
    out$gammaz <- out$gamma %*% t(dat$Z)
    out$estimates <- out$estimates + out$gammaz
  }
  if(check_mat(fit$beta0, FALSE)){
    if(trace) message("Constructing row intercepts ...")
    out$estimates <- out$estimates + fit$beta0 %*% matrix(1,1,ncol(out$estimates))
  }
  if(check_mat(fit$gamma0, FALSE)){
    if(trace) message("Constructing column intercepts ...")
    out$estimates <- out$estimates + matrix(1,nrow(out$estimates),1) %*% t(fit$gamma0)
  }
  if(trace) message("done.")
  return(out)
}

#-----------------------------
trim_eig <- function(d, tol=1e-10) d[d>0 & !is.nan(d) & !is.na(d)]

#' @export
mask_train_test_split <-
  function(obs_mask,
           testp = 0.2,
           seed = NULL) {
    # returns a new mask similar to mask with a new train-test sets.
    # testp: proportion of test out of the nonzero cells in mask
    # new_mask:
    #------------------------------------------------------------------
    #  0 -> train    (obs_mask=1) |
    #  1 -> test(valid)(obs_mask=1) |
    #  0 -> missing  (obs_mask=0)
    #-----------------------------------------------------------------
    if(!is.null(seed)) set.seed(seed)
    n_rows <- dim(obs_mask)[1]
    n_cols <- dim(obs_mask)[2]
    # Create a data frame of all matrix indices
    # we only consider non-missing data (ie, with mask_ij=1)
    indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)[obs_mask == 1, ]
    # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
    indices <-  indices[sample(1:nrow(indices)),]
    row.names(indices) <- NULL

    test.idx = sample(1:nrow(indices),
                      size = nrow(indices) * testp,
                      replace = FALSE)
    test.indices = indices[test.idx, ]

    new_mask = matrix(0, nrow = n_rows, ncol = n_cols)
    new_mask[as.matrix(test.indices[, c("row", "col")])] <- 1

    return(new_mask)
  }


#-------------------------
# Do not export
verify_low_rank <- function(M, J, min_eigv = 1e-6) {
  D <- M$d
  JD <- max(sum(D >= min_eigv), 1)
  if (JD >= J) {
    # more singular values than required: downscale
    M$u <- M$u[, seq(J), drop = FALSE]
    M$v <- M$v[, seq(J), drop = FALSE]
    M$d <- D[seq(J)]
  } else {
    # less singular values than required: upscale
    Ja <- J - JD
    M$d <- c(D, rep(D[JD], Ja))
    U <- M$u
    nr <- nrow(U)
    Ua <- matrix(stats::rnorm(nr * Ja), nr, Ja)
    Ua <- Ua - U %*% (t(U) %*% Ua)
    Ua <- svd_small_nc_cpp(Ua)$u
    M$u <- cbind(U, Ua)
    M$v <- cbind(M$v, matrix(0, nrow(M$v), Ja))
  }
  M
}

# Do not export
verify_warm_start <- function(M, J, min_eigv = 1e-6) {
  if (is.null(M)) {
    return(NULL)
  }
  d <- M$d
  if (is.null(d)) {
    return(NULL)
  }
  if (any(d > 0)) {
    if (length(d) == 1) {
      M$u <- matrix(M$u, ncol = 1)
      M$v <- matrix(M$v, ncol = 1)
    }
    verify_low_rank(M, J, min_eigv)
  } else {
    NULL
  }
}
error_metric = list(
  #--- error functions:
  unexplained_variance = function(predicted, true, adjusted = FALSE, k = NA) {
    # SSE / SST
    if(! adjusted){
      return(sum((true - predicted) ^ 2) / sum((true - mean(true)) ^ 2))
    }else{
      n = length(true)
      stopifnot(is.numeric(k))
      return(((sum((true - predicted) ^ 2) /
                 sum((true - mean(true)) ^ 2)) *
                (n - 1) / (n - k - 1)))
    }
  },

  mape = function(predicted, true) {
    mean(abs((as.double(true) - as.double(predicted)) / true), na.rm = TRUE) * 100
  },
  mae = function(predicted, true) {
    mean(abs(as.double(true) - as.double(predicted)), na.rm = TRUE)
  },

  rmse_normalized = function(predicted, true) {
    sqrt(mean((as.double(true) - as.double(predicted)) ^ 2, na.rm = TRUE)) /
      sd(as.double(true), na.rm = TRUE)
  },

  rmse = function(predicted, true) {
    sqrt(mean((as.double(true) - as.double(predicted)) ^ 2, na.rm = TRUE))
  },
  rel.rmse = function(predicted, true) {
    sqrt(sum((true - predicted) ^ 2, na.rm = TRUE)) /
      sqrt(sum(true^ 2, na.rm = TRUE))
  },

  spearman_R2 = function(predicted, true) {
    cor(true, predicted, method = "spearman")
  }
)



# SVD operations. general purpose, selects the optimal function to call
#' @export
opt_svd <-
  function(mat,
           k = NULL, # number of singular values to retain - default:return all
           nr = nrow(mat),
           nc = ncol(mat),
           rthin = nc > 2 * nr,
           cthin = nr > 2 * nc,
           trim = FALSE,
           tol = NULL) {
    if (is.null(k)) {
      if (rthin) {
        return(svd_small_nr_cpp(mat))
      }
      if (cthin) {
        return(svd_small_nc_cpp(mat))
      }
      return(base::svd(mat))
    }
    if (k == min(nr, nc)) {
      return(base::svd(mat))
    }
    if (rthin || cthin || k > 5) {
      return(RSpectra::svds(mat, k))
    }
    return(irlba::irlba(mat, k))
  }


# partial cross product at certain indices and returns a vector
#'
# partial_crossprod <-
#   function(u, v, irow, pcol, vtranpose = FALSE) {
#     dd <- dim(u)
#     nnrow <- as.integer(dd[1])
#     nrank <- dd[2]
#     stopifnot(all(irow < nnrow))
#     storage.mode(u) <- "double"
#     storage.mode(v) <- "double"
#     storage.mode(irow) <- "integer"
#     storage.mode(pcol) <- "integer"
#     nomega <- as.integer(length(irow))
#     r = double(nomega)
#
#     if (vtranpose) {
#       nncol <- as.integer(nrow(v))
#       stopifnot(nrank == ncol(v))
#       pcrossprodt_call(nnrow, nncol, nrank, u, v, irow, pcol, nomega, r)
#     } else {
#       nncol <- as.integer(ncol(v))
#       stopifnot(nrank == nrow(v))
#       pcrossprod_call(nnrow, nncol, nrank, u, v, irow, pcol, nomega, r)
#     }
#     return(r)
#
#     # .Fortran(
#     #   call_fun,
#     #   nnrow,
#     #   nncol,
#     #   nrank,
#     #   u,
#     #   v,
#     #   irow,
#     #   pcol,
#     #   nomega,
#     #   r = double(nomega),
#     #   PACKAGE = "IMR"
#     # )$r
#   }
