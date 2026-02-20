// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

inline Rcpp::NumericVector arma_vec_to_R(const arma::vec& x) {
  Rcpp::NumericVector out(x.n_elem);
  std::memcpy(out.begin(), x.memptr(), x.n_elem * sizeof(double));
  return out;
}

// ---- Map a base matrix or Matrix::dgeMatrix to arma::mat (no copy) ----
static arma::mat dense_view(SEXP mS) {
  if (Rf_isMatrix(mS) && TYPEOF(mS) == REALSXP) {
    IntegerVector dim = Rf_getAttrib(mS, R_DimSymbol);
    return arma::mat(REAL(mS), dim[0], dim[1], /*copy_aux_mem=*/false, /*strict=*/true);
  }
  if (Rf_isS4(mS)) {
    S4 s(mS);
    if (s.is("dgeMatrix")) {
      NumericVector x = s.slot("x");
      IntegerVector Dim = s.slot("Dim");
      return arma::mat(REAL(x), Dim[0], Dim[1], /*copy_aux_mem=*/false, /*strict=*/true);
    }
  }
  stop("svd_small_nc_cpp: expected a base numeric matrix or a Matrix::dgeMatrix.");
}

// [[Rcpp::export]]
Rcpp::List svd_small_nc_cpp(SEXP mS) {
  arma::mat M = dense_view(mS);       // n x p (no copy)
  const arma::uword p = M.n_cols;

  // 1. Gram matrix
  arma::mat B = M.t() * M;            // p x p

  // 2. Eigendecomposition
  arma::vec eval;                     // Guaranteed ascending
  arma::mat evec;
  if (!arma::eig_sym(eval, evec, B)) {
    stop("eig_sym failed on B = t(M) %*% M.");
  }

  // 3. Skip the sorting algorithm! Just flip them backward.
  eval = arma::flipud(eval);
  arma::mat V = arma::fliplr(evec);

  // 4. Safety clamp for floating-point rounding errors before sqrt
  eval.elem(arma::find(eval < 0.0)).zeros();
  arma::vec d = arma::sqrt(eval);

  // 5. Scale V instead of U (Saves n * p divisions!)
  arma::mat V_scaled = V;
  const double eps = std::numeric_limits<double>::epsilon();

  for (arma::uword j = 0; j < p; ++j) {
    if (d(j) > eps) {
      V_scaled.col(j) /= d(j);
    } else {
      V_scaled.col(j).zeros();
    }
  }

  // 6. One clean, fast BLAS Level 3 operation
  arma::mat U = M * V_scaled;

  // 7. Return directly
  Rcpp::NumericVector d_out(d.begin(), d.end());
  return List::create(_["d"] = d_out, _["v"] = V, _["u"] = U);
}


// [[Rcpp::export]]
Rcpp::List svd_small_nr_cpp(SEXP mS) {
  arma::mat M = dense_view(mS);     // n x p (no copy)
  const arma::uword n = M.n_rows;

  // 1. Gram matrix (n x n)
  arma::mat A = M * M.t();

  // 2. Eigendecomposition
  arma::vec eval;      // ascending
  arma::mat evec;      // columns = eigenvectors
  if (!arma::eig_sym(eval, evec, A)) {
    stop("eig_sym failed on M %*% t(M)");
  }

  // 3. Skip the sort! Reverse the arrays to get descending order.
  eval = arma::flipud(eval);
  arma::mat U = arma::fliplr(evec);

  // 4. Safety clamp for floating-point negative zeros
  eval.elem(arma::find(eval < 0.0)).zeros();
  arma::vec d = arma::sqrt(eval);

  // 5. Scale U instead of V (Saves p * n divisions!)
  arma::mat U_scaled = U;
  const double eps = std::numeric_limits<double>::epsilon();
  for (arma::uword j = 0; j < n; ++j) {
    if (d(j) > eps) {
      U_scaled.col(j) /= d(j);
    } else {
      U_scaled.col(j).zeros();
    }
  }

  // 6. One clean, fast BLAS Level 3 operation
  arma::mat V = M.t() * U_scaled;

  // 7. Format output
  Rcpp::NumericVector d_out(d.begin(), d.end());
  return List::create(_["d"] = d_out, _["u"] = U, _["v"] = V);
}
