// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
extern "C" {
#include <R_ext/RS.h>  // F77_CALL/F77_NAME
}

// Fortran symbols
extern "C" {
  void F77_NAME(pcrossprod)(int*,int*,int*,double*,double*,int*,int*,int*,double*);
  void F77_NAME(pcrossprodt)(int*,int*,int*,double*,double*,int*,int*,int*,double*);
}



// slower than R function
// [[Rcpp::export]]
Rcpp::NumericVector partial_crossprod(Rcpp::NumericMatrix u,
                                          Rcpp::NumericMatrix v,
                                          Rcpp::IntegerVector irow,
                                          Rcpp::IntegerVector pcol,
                                          bool vtranspose = false)
{
  int nrow  = u.nrow();
  int nrank  = u.ncol();
  int nomega = irow.size();

  if (nrow <= 0 || nrank <= 0)
    Rcpp::stop("u must be a non-empty numeric matrix");


  int ncol = vtranspose ? v.nrow() : v.ncol();
  if ( (vtranspose && nrank != v.ncol()) || (!vtranspose && nrank != v.nrow()) )
    Rcpp::stop(vtranspose ?
                 "When vtranspose=TRUE, ncol(u) must equal ncol(v)." :
                 "When vtranspose=FALSE, ncol(u) must equal nrow(v).");



  Rcpp::NumericVector r(nomega);

  if (vtranspose) {
    F77_CALL(pcrossprodt)(&nrow, &ncol, &nrank,
             u.begin(), v.begin(),
             irow.begin(), pcol.begin(),
             &nomega, r.begin());
  } else {
    F77_CALL(pcrossprod )(&nrow, &ncol, &nrank,
             u.begin(), v.begin(),
             irow.begin(), pcol.begin(),
             &nomega, r.begin());
  }
  return r;
}

// [[Rcpp::export]]
Rcpp::NumericVector partial_crossprod_fast(Rcpp::NumericMatrix u,
                                      Rcpp::NumericMatrix v,
                                      Rcpp::IntegerVector irow,
                                      Rcpp::IntegerVector pcol,
                                      bool vtranspose = false)
{
  int nrow  = u.nrow();
  int nrank  = u.ncol();
  int nomega = irow.size();

  if (nrow <= 0 || nrank <= 0)
    Rcpp::stop("u must be a non-empty numeric matrix");


  int ncol = vtranspose ? v.nrow() : v.ncol();
  if ( (vtranspose && nrank != v.ncol()) || (!vtranspose && nrank != v.nrow()) )
    Rcpp::stop(vtranspose ?
                 "When vtranspose=TRUE, ncol(u) must equal ncol(v)." :
                 "When vtranspose=FALSE, ncol(u) must equal nrow(v).");



  Rcpp::NumericVector r(nomega);

  if (vtranspose) {
    F77_CALL(pcrossprodt)(&nrow, &ncol, &nrank,
             u.begin(), v.begin(),
             irow.begin(), pcol.begin(),
             &nomega, r.begin());
  } else {
    F77_CALL(pcrossprod )(&nrow, &ncol, &nrank,
             u.begin(), v.begin(),
             irow.begin(), pcol.begin(),
             &nomega, r.begin());
  }
  return r;
}
