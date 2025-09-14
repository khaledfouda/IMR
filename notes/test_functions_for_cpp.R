# Helper: pure-R reference that mirrors CSC update on x using p (0-based)
ref_add_to_cols_inplace <- function(x, p, add_per_col, alpha = 1) {
  stopifnot(length(p) >= 1L, length(add_per_col) == (length(p) - 1L))
  nnz <- length(x)
  stopifnot(p[1] == 0L, p[length(p)] == nnz)

  m <- length(add_per_col)
  for (j in seq_len(m)) {
    start <- p[j] + 1L        # convert 0-based to 1-based for R's x
    end   <- p[j + 1L]        # inclusive in R indexing
    if (start <= end) {
      x[start:end] <- x[start:end] + alpha * add_per_col[j]
    }
  }
  x
}

# One-shot test you can run in the console
run_add_to_cols_inplace_cpp_test <- function() {
  # Build a 3x4 matrix with an empty 2nd column; p = c(0,2,2,4,6)
  # Column-wise nonzeros:
  #   C1: [1, 2], C2: [], C3: [3, 4], C4: [5, 6]
  Matrix::Matrix(
    c(1,0,2,
      0,0,0,
      3,4,0,
      0,5,6),
    nrow = 3, ncol = 4, byrow = FALSE, sparse = TRUE
  ) -> M

  add_per_col <- c(10, 20, 30, 0.5)
  alpha <- 2

  # Expected x after update
  expected_x <- ref_add_to_cols_inplace(
    x = M@x, p = M@p, add_per_col = add_per_col, alpha = alpha
  )

  # Try the C++ function; it modifies M@x in place if implemented correctly
  res <- try({
    add_to_cols_inplace_cpp(M@x, M@p, add_per_col, alpha)
    M@x
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    message("C++ call errored: ", conditionMessage(attr(res, "condition")))
    message(
      "This reproduces the bug: the guard `if (pj > 0 && pj < nnz)` rejects p[0] == 0.\n",
      "Fix the checks to enforce 0 <= p[j] <= p[j+1] <= nnz and re-run."
    )
    return(invisible(FALSE))
  }

  # Compare with reference
  equal <- isTRUE(all.equal(res, expected_x, tolerance = 1e-12))
  if (!equal) {
    max_diff <- max(abs(res - expected_x))
    stop(sprintf("Mismatch vs. reference; max abs diff = %.3g", max_diff))
  }

  message("PASS: C++ result matches the reference implementation.")
  invisible(TRUE)
}



run_partial_crossprod_test <- function(n=50,r=10,m=60,prob=0.3){
  mask <- matrix(rbinom(n*m,1,prob),n,m)
  A <- matrix(rnorm(n*r), n, r)
  B <- matrix(rnorm(m*r), r, m)
  mask <- as.Incomplete(mask)
  out1 <- partial_crossprod(A, B, mask@i, mask@p)
  out2 <- partial_crossprod(A, t(B), mask@i, mask@p, vtranspose = T)
  true <- (A %*% B)[as.matrix(mask)==1]
  if(all.equal(out1, true) & all.equal(out2,true)){
    message("PASS: partial_crossprod() test")
  }else(
    message("Mismatch with reference")
  )
}

run_soft_thresh_test <- function(n=50, m=10, a = 0.5){
  A <- matrix(rnorm(n*m), n, m)
  out <- soft_threshold_cpp(A, a)
  true <- sign(A) * pmax(abs(A) - a, 0)
  if(all.equal(out, true) ){
    message("PASS: soft_threshold() test")
  }else(
    message("Mismatch with reference")
  )
}






