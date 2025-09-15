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


// Computes B = t(M) %*% M (p×p), eigensolves it, builds U,V,d.
// [[Rcpp::export]]
Rcpp::List svd_small_nc_cpp(SEXP mS) {
  arma::mat M = dense_view(mS);         // n × p (no copy)
  const arma::uword p = M.n_cols;

  // B = crossprod(M)
  arma::mat B = M.t() * M;              // p × p, symmetric PSD

  // eigendecomposition (faster than SVD for symmetric B)
  arma::vec eval;                        // ascending
  arma::mat evec;
  if (!arma::eig_sym(eval, evec, B))
    stop("eig_sym failed on B = t(M) %*% M.");

  // sort descending by eigenvalue (singular values = sqrt(eigenvalues))
  arma::uvec ord = arma::sort_index(eval, "descend");
  arma::vec d = arma::sqrt(eval(ord));
  arma::mat V = evec.cols(ord);         // right singular vectors

  // U = M * V * diag(1/d)  (vectorized; guard zeros)
  arma::mat U = M * V;                  // n × p
  const double eps = std::numeric_limits<double>::epsilon();
  for (arma::uword j = 0; j < p; ++j) {
    double dj = d(j);
    if (dj > eps) U.col(j) /= dj;
    else          U.col(j).zeros();     // rank-deficient: safe/finite
  }
  // ensure that d is an R vector instead of an matrix
  Rcpp::NumericVector d_out(d.begin(), d.end());
  return List::create(_["d"] = d, _["v"] = V, _["u"] = U);
}

// [[Rcpp::export]]
Rcpp::List svd_small_nr_cpp(SEXP mS) {
  arma::mat M = dense_view(mS);     // n × p (no copy)
  //const arma::uword n = M.n_rows;

  // A = M %*% t(M)  (n × n, symmetric PSD)
  arma::mat A = M * M.t();

  // eigendecomposition (cheaper than SVD for symmetric A)
  arma::vec eval;      // ascending
  arma::mat evec;      // columns = eigenvectors
  if (!arma::eig_sym(eval, evec, A))
    stop("eig_sym failed on M %*% t(M)");

  // sort descending (singular values are sqrt(eigenvalues))
  arma::uvec ord = arma::sort_index(eval, "descend");
  arma::vec d = arma::sqrt(eval(ord));   // length n
  arma::mat U = evec.cols(ord);          // n × n (left singular vectors)

  // V = t(M) %*% U %*% diag(1/d)
  arma::mat V = M.t() * U;               // p × n
  const double eps = std::numeric_limits<double>::epsilon();
  for (arma::uword j = 0; j < d.n_elem; ++j) {
    double dj = d(j);
    if (dj > eps) V.col(j) /= dj;        // scale each column by 1/d_j
    else          V.col(j).zeros();      // rank-deficient guard
  }
  // ensure that d is an R vector instead of an matrix
  Rcpp::NumericVector d_out(d.begin(), d.end());
  return List::create(_["d"] = d_out, _["u"] = U, _["v"] = V);
}
