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
                                 const int n_rows,
                                 const int n_cols) {

  // 1. Initialize output (automatically filled with 0.0)
  NumericVector means(n_rows);

  // 2. Extract raw C pointers for zero-overhead access
  double* p_means = REAL(means);
  const double* p_x = REAL(x);
  const int* p_i = INTEGER(i);

  R_xlen_t nnz = x.size();

  // 3. Flat loop: Add each value to its corresponding row sum
  for (R_xlen_t k = 0; k < nnz; ++k) {
    p_means[ p_i[k] ] += p_x[k];
  }

  // 4. Divide by total columns to get the mean
  const double denom = static_cast<double>(n_cols);
  for (int r = 0; r < n_rows; ++r) {
    p_means[r] /= denom;
  }

  return means;
}

// [[Rcpp::export]]
NumericVector col_means_cpp(const NumericVector x,
                                 const IntegerVector p,
                                 const int n_rows,
                                 const int n_cols) {

  // 1. Initialize output
  NumericVector means(n_cols);

  // 2. Extract raw C pointers
  double* p_means = REAL(means);
  const double* p_x = REAL(x);
  const int* p_p = INTEGER(p);

  const double denom = static_cast<double>(n_rows);

  // 3. Column-centric loop
  for (int j = 0; j < n_cols; ++j) {
    double s = 0.0;
    int start = p_p[j];
    int end = p_p[j + 1];

    // Sum all non-zero elements in column j
    for (int k = start; k < end; ++k) {
      s += p_x[k];
    }

    p_means[j] = s / denom;
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

//
// // Build an arma::sp_mat from a Matrix::dgCMatrix (n x m)
//   static arma::sp_mat as_spmat_dgc(const S4& y) {
//     IntegerVector Dim = y.slot("Dim");
//     IntegerVector i   = y.slot("i");
//     IntegerVector p   = y.slot("p");
//     NumericVector x   = y.slot("x");
//
//     arma::uvec ai = arma::conv_to<arma::uvec>::from( as< std::vector<unsigned int> >(i) );
//     arma::uvec ap = arma::conv_to<arma::uvec>::from( as< std::vector<unsigned int> >(p) );
//     arma::vec  ax = as<arma::vec>(x);
//
//     return arma::sp_mat(ai, ap, ax, Dim[0], Dim[1]); // n x m
//   }
//
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

//
// arma::mat update_A_sim_cpp(SEXP yS4,              // dgCMatrix (n x m)
//                            const arma::mat& U,    // n x J
//                            const arma::mat& V,    // m x J
//                            const arma::vec& Dsq,  // length J
//                            const double lambda_M,
//                            const arma::mat& Ur,   // n x n
//                            const arma::vec& dr) { // length n
//   // y: sparse matrix
//   S4 y(yS4);
//   arma::sp_mat Y = as_spmat_dgc(y);     // n x m
//   arma::mat A = Ur.t() * ( Y * V + U.each_row() % Dsq.t());  // n x J
//   A = A.each_row() % Dsq.t();
//
//   for (arma::uword j = 0; j < U.n_cols; ++j) {
//     const double d = Dsq(j);
//     arma::vec a = 1 / (dr + d );
//     A.col(j) = Ur * (a % A.col(j));
//   }
//
//   return A; // n x J
// }



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
                        const arma::mat& Uc,   // m x m
                        const arma::vec& dc) { // length m
  // y: sparse matrix
  S4 y(yS4);
  arma::sp_mat Y = as_spmat_dgc(y);     // n x m
  arma::uword  J = U.n_cols;
  // B <- crossprod(Uc,(crossprod(Y, U) + V %*% Dsq))
  arma::mat B = Uc.t() * ( arma::trans(Y) * U + V.each_row() % Dsq.t());  // m x J

  for (arma::uword j = 0; j < J; ++j) {
    const double d = Dsq(j);
    arma::vec a = d / (dc + d );
    B.col(j) = Uc * (a % B.col(j));
  }

  return B; // m x J
}



