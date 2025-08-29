// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// lasso soft-thresholding
// [[Rcpp::export]]
NumericMatrix soft_threshold_cpp(const NumericMatrix B, const double lambda) {
  const int nrow = B.nrow();
  const int ncol = B.ncol();
  NumericMatrix out(nrow, ncol);

  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      const double val = B(i, j);
      double a = std::abs(val) - lambda;
      if (a < 0.0) a = 0.0;
      out(i, j) = (val > 0.0) ?  a : (val < 0.0 ? -a : 0.0);
    }
  }
  return out;
}

// the following two functions compute row and column means of
// a matrix of class Incomplete (sparse)
// faster than rowMeans and colMeans
// [[Rcpp::export]]
NumericVector row_means_cpp(SEXP yS4, const int n_cols) {
  S4 y(yS4);
  IntegerVector i = y.slot("i");
  IntegerVector p = y.slot("p");
  NumericVector x = y.slot("x");
  IntegerVector dim = y.slot("Dim");
  const int n = dim[0];
  const int m = dim[1];

  NumericVector sums(n); // zero-initialized
  for (int j = 0; j < m; ++j) {
    for (int k = p[j]; k < p[j + 1]; ++k) {
      sums[ i[k] ] += x[k];
    }
  }
  const double denom = static_cast<double>(n_cols);
  for (int r = 0; r < n; ++r) sums[r] /= denom;
  return sums;
}
// [[Rcpp::export]]
NumericVector col_means_cpp(SEXP yS4, const int n_rows) {
  S4 y(yS4);
  IntegerVector p = y.slot("p");
  NumericVector x = y.slot("x");
  const int m = p.size() - 1;

  NumericVector means(m);
  const double denom = static_cast<double>(n_rows);
  for (int j = 0; j < m; ++j) {
    double s = 0.0;
    for (int k = p[j]; k < p[j + 1]; ++k) s += x[k];
    means[j] = s / denom;
  }
  return means;
}

// the following two functions add a vector to the rows or the columns
// of the R sparse matrix. yx = y@x; i = y@i; p = y@p
// the update is done inplace
// [[Rcpp::export]]
void add_to_rows_inplace_cpp(NumericVector yx,
                             const IntegerVector i,
                             const NumericVector add_per_row,
                             double alpha = 1.0)
{
  const R_xlen_t nr = add_per_row.size();
  const R_xlen_t nnz = yx.size();
  if (i.size() != nnz)
    stop("Length mismatch: length(i)=%ld must equal length(yx)=%ld",
         (long)i.size(), (long)nnz);


  for (R_xlen_t k = 0; k < nnz; ++k) {
    int r = i[k];
    if (r >= 0 && r < nr)          // guard
      yx[k] += alpha * add_per_row[r];
    else
      stop("row index out of range");
  }
}
// [[Rcpp::export]]
void add_to_cols_inplace_cpp(NumericVector yx,
                             const IntegerVector p,
                             const NumericVector add_per_col,
                             const double alpha = 1.0) {
  const int m = p.size() - 1;
  const R_xlen_t nnz = yx.size();
  for (int j = 0; j < m; ++j) {
    const int pj = p[j];
    const double v = alpha * add_per_col[j];
    if(pj > 0 && pj < nnz)
      for (int k = pj; k < p[j + 1]; ++k) yx[k] += v;
    else
      stop("Out of bound");
  }
}

// the following functions computes the frobonanci ratio
// between old and new svd decompositions
// [[Rcpp::export]]
double frob_ratio_cpp(const arma::mat& Uold,  const arma::vec& Dsqold, const arma::mat& Vold,
                      const arma::mat& U,     const arma::vec& Dsq,    const arma::mat& V) {
  const arma::uword r = Dsq.n_elem;
  if (U.n_cols != r || V.n_cols != r ||
      Uold.n_cols != r || Vold.n_cols != r ||
      Dsqold.n_elem != r) {
    Rcpp::stop("Dimension mismatch: ranks must agree.");
  }

  // r x r Gram blocks (BLAS-backed)
  arma::mat GU = U.t() * Uold;   // crossprod(U, Uold)
  arma::mat GV = Vold.t() * V;   // crossprod(Vold, V)

  // uvprod = sum_{i,j} Dsq[i] * GU(i,j) * Dsqold[j] * GV(j,i)
  double uvprod = 0.0;
  for (arma::uword i = 0; i < r; ++i) {
    const double di = Dsq(i);
    for (arma::uword j = 0; j < r; ++j) {
      uvprod += di * GU(i, j) * (Dsqold(j) * GV(j, i));
    }
  }

  const double denom = arma::dot(Dsqold, Dsqold);  // sum(Dsqold^2)
  const double num   = denom + arma::dot(Dsq, Dsq) - 2.0 * uvprod;
  const double denom_safe = (denom > 1e-9) ? denom : 1e-9;

  return num / denom_safe;
}

// Build an arma::sp_mat from a Matrix::dgCMatrix (n x m)
  static arma::sp_mat as_spmat_dgc(const S4& y) {
    IntegerVector Dim = y.slot("Dim");
    IntegerVector i   = y.slot("i");
    IntegerVector p   = y.slot("p");
    NumericVector x   = y.slot("x");

    arma::uvec ai = arma::conv_to<arma::uvec>::from( as< std::vector<unsigned int> >(i) );
    arma::uvec ap = arma::conv_to<arma::uvec>::from( as< std::vector<unsigned int> >(p) );
    arma::vec  ax = as<arma::vec>(x);

    return arma::sp_mat(ai, ap, ax, Dim[0], Dim[1]); // n x m
  }

// The following two functions compute the least-squares updates for A and B

// A_mat = (y %*% V + U * diag(Dsq)) %*% diag(D_star)
// with D_star_j = Dsq_j / (Dsq_j + lambda_M)
// [[Rcpp::export]]
arma::mat update_A_cpp(SEXP yS4,              // dgCMatrix (n x m)
                            const arma::mat& V,    // m x J
                            const arma::mat& U,    // n x J
                            const arma::vec& Dsq,  // length J
                            const double lambda_M) {
  S4 y(yS4);
  arma::sp_mat Y = as_spmat_dgc(y);               // n x m
  const arma::uword n = Y.n_rows;
  const arma::uword m = Y.n_cols;

  if (V.n_rows != m) Rcpp::stop("camc_update_A_cpp: V has wrong nrow");
  if (U.n_rows != n) Rcpp::stop("camc_update_A_cpp: U has wrong nrow");
  if (U.n_cols != V.n_cols || U.n_cols != Dsq.n_elem)
    Rcpp::stop("camc_update_A_cpp: J mismatch among U, V, Dsq");

  const arma::uword J = U.n_cols;

  // A = Y %*% V  (n x J)
  arma::mat A = Y * V;

  //  for each j, A[,j] = D_star_j * ( A[,j] + Dsq_j * U[,j] )
  for (arma::uword j = 0; j < J; ++j) {
    const double d   = Dsq(j);
    const double dst = d / (d + lambda_M);
    A.col(j) = dst * (A.col(j) + d * U.col(j));
  }

  return A; // n x J
}
// [[Rcpp::export]]
arma::mat update_B_cpp(SEXP yS4,              // dgCMatrix (n x m)
                            const arma::mat& U,    // n x J
                            const arma::mat& V,    // m x J
                            const arma::vec& Dsq,  // length J
                            const double lambda_M) {
  S4 y(yS4);
  arma::sp_mat Y = as_spmat_dgc(y);     // n x m
  arma::uword n = Y.n_rows, m = Y.n_cols;

  if (U.n_rows != n) stop("U has wrong nrow");
  if (V.n_rows != m) stop("V has wrong nrow");
  if (U.n_cols != V.n_cols || U.n_cols != Dsq.n_elem) stop("J mismatch");

  const arma::uword J = U.n_cols;

  // B = t(Y) %*% U  (m x J)
  arma::mat B = arma::trans(Y) * U;

  //  B[,j] = D_star_j * ( B[,j] + Dsq_j * V[,j] )
  for (arma::uword j = 0; j < J; ++j) {
    const double d   = Dsq(j);
    const double dst = d / (d + lambda_M);
    B.col(j) = dst * (B.col(j) + d * V.col(j));
  }
  return B; // m x J
}
