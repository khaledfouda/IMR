// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix soft_threshold_cpp(const NumericMatrix B, const double lambda) {

  // 1. Create output matrix.
  // Rcpp automatically fills this entirely with 0.0 upon creation!
  NumericMatrix out(B.nrow(), B.ncol());

  // 2. Extract raw C pointers for zero-overhead access
  const double* p_B = REAL(B);
  double* p_out = REAL(out);

  // Total number of elements in the matrix
  R_xlen_t n = B.length();

  // 3. The Flat Loop
  for (R_xlen_t k = 0; k < n; ++k) {
    double val = p_B[k];

    // Simple, highly predictable branches.
    if (val > lambda) {
      p_out[k] = val - lambda;
    } else if (val < -lambda) {
      p_out[k] = val + lambda;
    }
    // Note: We do NOTHING if it's between -lambda and lambda.
    // p_out[k] is already 0.0, so we save the CPU a write operation!
  }

  return out;
}

// the following two functions compute row and column means of
// a matrix of class Incomplete (sparse)
// faster than rowMeans and colMeans
// [[Rcpp::export]]
NumericVector row_means_cpp(const NumericVector x,
                                 const IntegerVector i,
                                 const int n_rows) {

  // 1. Initialize output (automatically filled with 0.0)
  NumericVector means(n_rows);
  IntegerVector counts(n_rows);

  // 2. Extract raw C pointers for zero-overhead access
  double* p_means = REAL(means);
  int* p_counts = INTEGER(counts);
  const double* p_x = REAL(x);
  const int* p_i = INTEGER(i);

  R_xlen_t nnz = x.size();

  // 3. Flat loop: Add each value to its corresponding row sum
  for (R_xlen_t k = 0; k < nnz; ++k) {
    int row_idx = p_i[k];
    p_means[ row_idx ] += p_x[k];
    p_counts[row_idx] += 1;
  }

  // 4. Divide by total columns to get the mean
  // const double denom = static_cast<double>(n_cols);
  // for (int r = 0; r < n_rows; ++r) {
  //   p_means[r] /= denom;
  // }
  for (int r = 0; r < n_rows; ++r) {
    if (p_counts[r] > 0) {
      p_means[r] /= p_counts[r];
    } else {
      p_means[r] = NA_REAL; // Or 0.0, depending on how you want to handle entirely empty rows
    }
  }

  return means;
}

// [[Rcpp::export]]
NumericVector col_means_cpp(const NumericVector x,
                                 const IntegerVector p,
                                 const int n_cols) {

  // 1. Initialize output
  NumericVector means(n_cols);

  // 2. Extract raw C pointers
  double* p_means = REAL(means);
  const double* p_x = REAL(x);
  const int* p_p = INTEGER(p);

  //const double denom = static_cast<double>(n_rows);

  // 3. Column-centric loop
  for (int j = 0; j < n_cols; ++j) {
    double s = 0.0;
    int start = p_p[j];
    int end = p_p[j + 1];

    // Sum all non-zero elements in column j
    // for (int k = start; k < end; ++k) {
    //   s += p_x[k];
    // }
    //
    // p_means[j] = s / denom;
    // The exact number of observed values in this column!
    int count = end - start;

    if (count > 0) {
      // Sum all non-zero elements in column j
      for (int k = start; k < end; ++k) {
        s += p_x[k];
      }
      // Divide by the actual count of observed values
      p_means[j] = s / static_cast<double>(count);
    } else {
      // Handle the edge case where a column is completely empty
      p_means[j] = NA_REAL;
    }
  }

  return means;
}


// the following two functions add a vector to the rows or the columns
// of the R sparse matrix. yx = y@x; i = y@i; p = y@p
// the update is done inplace
// [[Rcpp::export]]
void add_to_rows_inplace_cpp(NumericVector yx,
                             const IntegerVector i,
                             const NumericVector add_per_row)
{
  // [later add this for speed]

  double* p_yx = REAL(yx);
  const int* p_i = INTEGER(i);
  const double* p_add = REAL(add_per_row);

  R_xlen_t n = yx.size();


  for (R_xlen_t k = 0; k < n; ++k) {
      p_yx[k] += p_add[p_i[k]];  //add_per_row[i[k]];
  }
}

// [[Rcpp::export]]
void add_to_cols_inplace_cpp(NumericVector yx,
                             const IntegerVector p,
                             const NumericVector add_per_col) {
  // const int m = p.size() - 1;
  //
  // for (int j = 0; j < m; ++j) {
  //   for (int k = p[j]; k < p[j+1]; ++k)
  //     yx[k] += add_per_col[j];
  // }
  double* p_yx = REAL(yx);
  const int* p_p = INTEGER(p);
  const double* p_add = REAL(add_per_col);

  int m = p.size() - 1;

  for (int j = 0; j < m; ++j) {
    // p_p[j] to p_p[j+1] gives the start and end indices for column j
    int start = p_p[j];
    int end = p_p[j+1];
    double val_to_add = p_add[j];

    for (int k = start; k < end; ++k) {
      p_yx[k] += val_to_add;
    }
  }
}




// [[Rcpp::export]]
double frob_ratio_cpp(const arma::mat& Uold,  const arma::vec& Dsqold, const arma::mat& Vold,
                      const arma::mat& U,     const arma::vec& Dsq,    const arma::mat& V) {

  const arma::uword r = Dsq.n_elem;
  const arma::uword ro = Dsqold.n_elem;

  // 1. BLAS-backed matrix multiplications (The heavy lifting)
  arma::mat GU = U.t() * Uold;   // size r x ro
  arma::mat GV = Vold.t() * V;   // size ro x r

  // 2. Compute trace directly (Loop order optimized for Column-Major access)
  double uvprod = 0.0;
  for (arma::uword j = 0; j < ro; ++j) {
    const double dj = Dsqold(j);
    for (arma::uword i = 0; i < r; ++i) {
      // Inner loop 'i' means we move down the columns of GU, which is contiguous in RAM.
      uvprod += Dsq(i) * GU(i, j) * dj * GV(j, i);
    }
  }

  // 3. Compute final sequence
  const double denom = arma::dot(Dsqold, Dsqold);
  const double num   = denom + arma::dot(Dsq, Dsq) - 2.0 * uvprod;

  return num / std::max(denom, 1e-9);
}

// --------------------------------------------------------------------------
// Huber helpers
//
// rho_c(u) = u^2/2 for |u| <= c
//          = c |u| - c^2/2 otherwise
// psi_c(u) = rho_c'(u) = max(-c, min(c,u)) # we apply this to the residuals
//
// `update_huber_c_cpp` computes the tuning constant, c.
// `huber_clip_cpp` apply psi_c element-wise to an input vector.
// ------------------------------------------------------------------

// Inter-quartile range
// R equivalence: stats::IQR(x, type=7) = quantile(x, 3/4, type=7) - quantile(x, 1/4, type=7)
static double iqr_type7_(std::vector<double>&  v) {

  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(v.size());
  if(n < 3) return 0; // must have more than 3 observations for peace of mind

  // from stats::IQR docs we have k = 1 + (n-1)p
  // we subtract 1 because cpp is 0-based. so we compute k = (n-1)p
  const double k_lo = 0.25 * static_cast<double>(n-1);
  const double k_hi = 0.75 * static_cast<double>(n-1);
  // j is the floor of k. j = floor(np+m)
  const std::ptrdiff_t j_lo = static_cast<std::ptrdiff_t>(std::floor(k_lo));
  const std::ptrdiff_t j_hi = static_cast<std::ptrdiff_t>(std::floor(k_hi));
  // g = k - j in [0,1) (g is gamma in stats::IQR definition)
  const double g_lo = k_lo - static_cast<double>(j_lo);
  const double g_hi = k_hi - static_cast<double>(j_hi);

  // Q(p) = (1-g) v[i] + g v[i+1]
  double a_lo, b_lo, a_hi, b_hi;
  // compute upper quartile
  // partition the vector and leaves the j_hi smallest elements in [0, j_hi).
  // partial sort. we only need j_hi, h_lo+1 to be in place.
  std::nth_element(v.begin(), v.begin() + j_hi, v.end());
  a_hi = v[j_hi];
  b_hi = (g_hi > 0.0)
    ? *std::min_element(v.begin() + j_hi + 1, v.end())
      : a_hi;

  // lower quartile inside [0, j_hi + 1)
  std::nth_element(v.begin(), v.begin() + j_lo + 1, v.begin() + j_hi + 1);
  b_lo = v[j_lo + 1];
  a_lo = *std::max_element(v.begin(), v.begin() + j_lo + 1);

  // return (1-g)*a + g*b
  const double q_lo = (g_lo > 0.0 && b_lo != a_lo)
    ? (1.0 - g_lo) * a_lo + g_lo * b_lo : a_lo;
  const double q_hi = (g_hi > 0.0 && b_hi != a_hi)
    ? (1.0 - g_hi) * a_hi + g_hi * b_hi : a_hi;

  return q_hi - q_lo;

}

// update Huber c .. fill description later

// [[Rcpp::export]]
double update_huber_c_cpp(const NumericVector yx,
                          const double huber_shift,
                          const double c_old,
                          const int max_sample = 100000) {

  if (ISNAN(huber_shift) || ISNAN(c_old))
    Rcpp::stop("update_huber_c_cpp: `huber_shift` and `c_old` must not be NA/NaN.");

  const R_xlen_t n = yx.size();
  const double* p_x = REAL(yx);

  //std::vector<double> v(p_x, p_x + n);

  // nth_element permutes its input, so a copy is mandatory; the copy is a memcpy
  std::vector<double> v;

  if (max_sample <= 0 || n <= static_cast<R_xlen_t>(max_sample)) {
    // use all data
    v.assign(p_x, p_x + n);

  } else {

    // work on a sample instead
    // Sampling is with replacement for a lower computational cost
    // it shouldn't be a problem for very large vectors which is what sampling is intended for
    const R_xlen_t m = static_cast<R_xlen_t>(max_sample);
    const double dn = static_cast<double>(n);
    const double phi_inv = 0.6180339887498948482;

    v.resize(static_cast<std::size_t>(m));

    for (R_xlen_t t = 0; t < m; ++t) {
      double u = static_cast<double>(t) * phi_inv;
      u -= std::floor(u);
      R_xlen_t idx = static_cast<R_xlen_t>(u* dn);
      if (idx >= n) idx = n - 1;
      v[static_cast<std::size_t>(t)] = p_x[idx];
    }
  }

  const double d = iqr_type7_(v) / 1.349;

  const double cand = huber_shift * d;
  return (cand < c_old) ? cand : c_old;
}

// clip residual vector to the Huber chosen constant
// [[Rcpp::export]]
NumericVector huber_clip_cpp(const NumericVector yx, const double huber_c) {

  if (ISNAN(huber_c)) Rcpp::stop("huber_clip_cpp: `huber_c` must not be NA/NaN.");
  if (huber_c < 0.0)  Rcpp::stop("huber_clip_cpp: `huber_c` must be non-negative.");

  const R_xlen_t n = yx.size();
  NumericVector out(Rcpp::no_init(n));

  const double* p_x   = REAL(yx);
  double*       p_out = REAL(out);
  const double  neg_c = -huber_c;

  for (R_xlen_t k = 0; k < n; ++k) {
    const double val = p_x[k];
    // Both comparisons are false for NA/NaN, so `val` passes through intact.
    p_out[k] = (val > huber_c) ? huber_c : ((val < neg_c) ? neg_c : val);
  }

  return out;
}

// [[Rcpp::export]]
void huber_clip_into_cpp(const NumericVector yx,
                         const double huber_c,
                         NumericVector out) {

  if (ISNAN(huber_c)) Rcpp::stop("huber_clip_into_cpp: `huber_c` must not be NA/NaN.");
  if (huber_c < 0.0)  Rcpp::stop("huber_clip_into_cpp: `huber_c` must be non-negative.");
  if (out.size() != yx.size())
    Rcpp::stop("huber_clip_into_cpp: `out` and `yx` must have the same length.");

  const double* p_x   = REAL(yx);
  double*       p_out = REAL(out);
  const double  neg_c = -huber_c;
  const R_xlen_t n    = yx.size();

  for (R_xlen_t k = 0; k < n; ++k) {
    const double val = p_x[k];
    p_out[k] = (val > huber_c) ? huber_c : ((val < neg_c) ? neg_c : val);
  }
}

// [[Rcpp::export]]
void huber_clip_inplace_cpp(NumericVector yx, const double huber_c) {

  if (ISNAN(huber_c)) Rcpp::stop("huber_clip_inplace_cpp: `huber_c` must not be NA/NaN.");
  if (huber_c < 0.0)  Rcpp::stop("huber_clip_inplace_cpp: `huber_c` must be non-negative.");

  double*      p_yx  = REAL(yx);
  const double neg_c = -huber_c;
  const R_xlen_t n   = yx.size();

  for (R_xlen_t k = 0; k < n; ++k) {
    const double val = p_yx[k];
    p_yx[k] = (val > huber_c) ? huber_c : ((val < neg_c) ? neg_c : val);
  }
}




// The following two functions compute the least-squares updates for A and B


static arma::sp_mat as_spmat_dgc(const S4& y) {
  IntegerVector Dim = y.slot("Dim");

  // Cast directly from R memory to arma::ivec, then convert to arma::uvec
  arma::uvec ai = arma::conv_to<arma::uvec>::from(as<arma::ivec>(y.slot("i")));
  arma::uvec ap = arma::conv_to<arma::uvec>::from(as<arma::ivec>(y.slot("p")));
  arma::vec  ax = as<arma::vec>(y.slot("x"));

  return arma::sp_mat(ai, ap, ax, Dim[0], Dim[1]);
}

// [[Rcpp::export]]
arma::mat update_A_cpp(SEXP yS4,
                       const arma::mat& U,
                       const arma::mat& V,
                       const arma::vec& Dsq,
                       const double lambda_M) {

  // 1. Build sparse matrix efficiently
  arma::sp_mat Y = as_spmat_dgc(yS4);

  // 2. Let Armadillo's highly-optimized engine do the math
  arma::mat A = Y * V + U.each_row() % Dsq.t();

  // 3. Scale columns
  for (arma::uword j = 0; j < U.n_cols; ++j) {
    A.col(j) *= (1.0 / (1.0 + lambda_M / Dsq(j)));
  }

  return A;
}



// [[Rcpp::export]]
arma::mat update_A_sim_cpp(SEXP yS4,              // dgCMatrix (n x m)
                                const arma::mat& U,    // n x J
                                const arma::mat& V,    // m x J
                                const arma::vec& Dsq,  // length J
                                const arma::mat& Ur,   // n x n
                                const arma::vec& dr) { // length n

  // 1. Fast sparse construction (zero std::vector copies)
  arma::sp_mat Y = as_spmat_dgc(yS4);

  // 2. Construct W
  arma::mat W = Y * V;
  W = W.each_row() % Dsq.t() + U.each_row() % (Dsq % Dsq).t();

  // 3. Pre-multiply by t(Ur) (Level 3 BLAS Matrix-Matrix multiply)
  arma::mat W_tilde = Ur.t() * W;

  // 4. Create an intermediate matrix to hold the scaled columns
  // arma::no_zeros tells C++ not to waste time filling this with 0.0s since we
  // are going to immediately overwrite it anyway.
  arma::mat A_inner(U.n_rows, U.n_cols, arma::fill::none);

  // 5. The Loop: ONLY do the fast element-wise math here (O(n) operations)
  for (arma::uword j = 0; j < U.n_cols; ++j) {
    const double d = Dsq(j);
    arma::vec a = 1.0 / (dr + d);

    A_inner.col(j) = a % W_tilde.col(j);
  }

  // 6. Final transformation: ONE big Matrix-Matrix multiply (Level 3 BLAS)
  return Ur * A_inner;
}


// [[Rcpp::export]]
arma::mat update_B_cpp(SEXP yS4,              // dgCMatrix (n x m)
                       const arma::mat& U,    // n x J
                       const arma::mat& V,    // m x J
                       const arma::vec& Dsq,  // length J
                       const double lambda_M) {

  // 1. Build sparse matrix efficiently with zero std::vector copies
  arma::sp_mat Y = as_spmat_dgc(yS4);     // n x m

  // 2. Let Armadillo's optimized engine handle the transposed sparse multiplication
  // B = t(Y) %*% U (m x J) + V * diag(Dsq)
  arma::mat B = Y.t() * U + V.each_row() % Dsq.t();

  // 3. Apply the simplified diagonal scaling
  for (arma::uword j = 0; j < U.n_cols; ++j) {
    B.col(j) *= (1.0 / (1.0 + lambda_M / Dsq(j)));
  }

  return B; // m x J
}



// [[Rcpp::export]]
arma::mat update_B_sim_cpp(SEXP yS4,              // dgCMatrix (n x m)
                           const arma::mat& U,    // n x J
                           const arma::mat& V,    // m x J
                           const arma::vec& Dsq,  // length J
                           const arma::mat& Uc,   // n x n
                           const arma::vec& dc) { // length n

  // 1. Fast sparse construction (zero std::vector copies)
  arma::sp_mat Y = as_spmat_dgc(yS4);

  // 2. Construct W
  arma::mat W = Y.t() * U;
  W = W.each_row() % Dsq.t() + V.each_row() % (Dsq % Dsq).t();

  // 3. Pre-multiply by t(Ur) (Level 3 BLAS Matrix-Matrix multiply)
  arma::mat W_tilde = Uc.t() * W;

  // 4. Create an intermediate matrix to hold the scaled columns
  // arma::no_zeros tells C++ not to waste time filling this with 0.0s since we
  // are going to immediately overwrite it anyway.
  arma::mat B_inner(V.n_rows, V.n_cols, arma::fill::none);

  // 5. The Loop: ONLY do the fast element-wise math here (O(n) operations)
  for (arma::uword j = 0; j < U.n_cols; ++j) {
    const double d = Dsq(j);
    arma::vec b = 1.0 / (dc + d);

    B_inner.col(j) = b % W_tilde.col(j);
  }

  // 6. Final transformation: ONE big Matrix-Matrix multiply (Level 3 BLAS)
  return Uc * B_inner;
}





