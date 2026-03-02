#' Do not export
inv <- function(X, tol = sqrt(.Machine$double.eps)) {
  # Ensure X is a base matrix
  if (!is.matrix(X)) X <- as.matrix(X)

  is_square <- nrow(X) == ncol(X)

  if (is_square) {
    #  Try Cholesky for symmetric matrices
    if (isSymmetric(X)) {
      tryCatch(
        {
          return(chol2inv(chol(X)))
        },
        error = function(e) {
          # not positive-definite
        }
      )
    }
    # standard solve or ginv if it fails
    tryCatch(
      {
        return(solve(X))
      },
      error = function(e) {
        # if singular
        return(MASS::ginv(X, tol = tol))
      }
    )
  }

  # ginv for rectangle matrices
  return(MASS::ginv(X, tol = tol))
}

#-------------------------------------
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
    if (!is.null(seed)) set.seed(seed)
    n_rows <- dim(obs_mask)[1]
    n_cols <- dim(obs_mask)[2]
    # Create a data frame of all matrix indices
    # we only consider non-missing data (ie, with mask_ij=1)
    indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)[obs_mask == 1, ]
    # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
    indices <- indices[sample(1:nrow(indices)), ]
    row.names(indices) <- NULL

    test.idx <- sample(1:nrow(indices),
      size = nrow(indices) * testp,
      replace = FALSE
    )
    test.indices <- indices[test.idx, ]

    new_mask <- matrix(0, nrow = n_rows, ncol = n_cols)
    new_mask[as.matrix(test.indices[, c("row", "col")])] <- 1

    return(new_mask)
  }


#-------------------------
# Do not export
verify_low_rank <- function(M, J, min_eigv = 1e-6) {
  D <- M$d
  # Count valid singular values, cap at length(D)
  JD <- sum(D >= min_eigv)
  # Ensure JD is at least 1 so we don't break matrix operations
  JD <- max(JD, 1)

  if (JD >= J) {
    # We have enough components: just truncate exactly to J
    M$u <- M$u[, seq(J), drop = FALSE]
    M$v <- M$v[, seq(J), drop = FALSE]
    M$d <- D[seq(J)]
  } else {
    # We don't have enough valid components: upscale to exactly J
    # First, keep ONLY the JD valid components
    U <- M$u[, seq(JD), drop = FALSE]
    V <- M$v[, seq(JD), drop = FALSE]
    D_valid <- D[seq(JD)]

    #  Add Ja new components
    Ja <- J - JD
    M$d <- c(D_valid, rep(D_valid[JD], Ja))

    nr <- nrow(U)
    Ua <- matrix(stats::rnorm(nr * Ja), nr, Ja)

    #  orthogonalization
    Ua <- Ua - U %*% (t(U) %*% Ua)

    #  Orthonormalize the new columns
    Ua <- qr.Q(qr(Ua))

    M$u <- cbind(U, Ua)
    M$v <- cbind(V, matrix(0, nrow(V), Ja))
  }

  return(M)
}

# Do not export
verify_warm_start <- function(M, J, min_eigv = 1e-16) {
  if (is.null(M) || is.null(M$d) || is.null(M$u) || is.null(M$v)) {
    warning("warm start verification failed - missing u, d, or v. Reinitializing...")
    return(NULL)
  }

  d <- M$d
  if (any(d >= min_eigv)) {
    # Fix R's habit of dropping dimensions to vectors
    if (length(d) == 1) {
      M$u <- as.matrix(M$u)
      M$v <- as.matrix(M$v)
    }
    return(verify_low_rank(M, J, min_eigv))
  } else {
    warning("warm start verification failed - no significant singular values. Reinitializing...")
    return(NULL)
  }
}


#' @export
error_metrics <- list(
  # Mean Absolute Percentage Error
  mape = function(predicted, actual, na.rm = TRUE) {
    if (na.rm) {
      idx <- complete.cases(predicted, actual) & actual != 0
    } else {
      idx <- actual != 0
    }
    mean(abs((actual[idx] - predicted[idx]) / actual[idx])) * 100
  },

  #  Mean Absolute Error
  mae = function(predicted, actual, na.rm = TRUE) {
    mean(abs(actual - predicted), na.rm = na.rm)
  },

  # RMSE
  rmse = function(predicted, actual, na.rm = TRUE) {
    sqrt(mean((actual - predicted)^2, na.rm = na.rm))
  },

  # 6. Relative RMSE
  rrmse = function(predicted, actual, na.rm = TRUE) {
    if (na.rm) {
      idx <- complete.cases(predicted, actual)
      predicted <- predicted[idx]
      actual <- actual[idx]
    }
    sqrt(sum((actual - predicted)^2)) / sqrt(sum(actual^2))
  },

  # 7. Spearman Correlation (Rho)
  spearman = function(predicted, actual, na.rm = TRUE) {
    use_method <- if (na.rm) "complete.obs" else "everything"
    cor(actual, predicted, method = "spearman", use = use_method)
  }
)

#' Return one of the error metrics above.
#' Do not export
get_metric <- function(metric) {
  if (stringr::str_to_lower(metric) == "rmse") {
    return(error_metrics$rmse)
  }
  if (stringr::str_to_lower(metric) == "rrmse") {
    return(error_metrics$rrmse)
  }
  if (stringr::str_to_lower(metric) == "mae") {
    return(error_metrics$mae)
  }
  if (stringr::str_to_lower(metric) == "mape") {
    return(error_metrics$mape)
  }
  if (stringr::str_to_lower(metric) == "spearman") {
    return(error_metric$spearman)
  }
  stop("Unvalid error metric")
}

#' Evaluate Model Predictions on all metrics
#' @export
evaluate <- function(predicted, actual, metric = "all", na.rm = TRUE) {
  p <- as.numeric(predicted)
  a <- as.numeric(actual)
  if (stringr::str_to_lower(metric) != "all") {
    return(IMR:::get_metric(metric)(p, a, na.rm))
  }
  tibble(
    RMSE            = error_metrics$rmse(p, a, na.rm),
    Rel_RMSE        = error_metrics$rrmse(p, a, na.rm),
    MAE             = error_metrics$mae(p, a, na.rm),
    MAPE            = error_metrics$mape(p, a, na.rm),
    Spearman_Rho    = error_metrics$spearman(p, a, na.rm)
  )
}


# SVD operations. general purpose, selects the optimal function to call
#' @export
svd_opt <- function(mat,
                    k = NULL,
                    tol = NULL) {
  nr <- nrow(mat)
  nc <- ncol(mat)

  # thin thresholds
  rthin <- nc > 2 * nr
  cthin <- nr > 2 * nc

  # If k is requested but equals the max possible rank, treat as full SVD
  if (!is.null(k) && k >= min(nr, nc)) {
    k <- NULL
  }

  # ---  FAST PATH: Thin Matrices or Full SVD ---
  # If the matrix is thin, use fastsvd
  if (rthin || cthin || is.null(k)) {
    if (rthin) {
      dec <- IMR:::svd_small_nr_cpp_fast(mat)
    } else if (cthin) {
      dec <- IMR:::svd_small_nc_cpp_fast(mat)
    } else {
      dec <- base::svd(mat)
    }

    # If a specific k was requested, truncate the full exact SVD
    if (!is.null(k)) {
      dec$u <- dec$u[, seq_len(k), drop = FALSE]
      dec$d <- dec$d[seq_len(k)]
      dec$v <- dec$v[, seq_len(k), drop = FALSE]
    }
  } else {
    #  large matrices with k != NULL
    # irlba is faster for very small k or sparse matrices.
    # RSpectra is faster for larger k or dense matrices.
    if (inherits(mat, "sparseMatrix") || k <= 5) {
      dec <- irlba::irlba(mat, nv = k)
    } else {
      dec <- RSpectra::svds(mat, k)
    }
  }

  # ---  Finally, truncate eigenvalues, if requested
  if (!is.null(tol)) {
    valid_k <- sum(dec$d > tol)

    if (valid_k < length(dec$d)) {
      dec$u <- dec$u[, seq_len(valid_k), drop = FALSE]
      dec$d <- dec$d[seq_len(valid_k)]
      dec$v <- dec$v[, seq_len(valid_k), drop = FALSE]
    }
  }

  return(dec)
}

#
