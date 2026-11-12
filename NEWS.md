# IMR (development version)

* New `imr()`: prepares the data, computes the grid limits, tunes and fits the
  model in one call, with the most used arguments of `imr_data()`, `update()`,
  `imr_set_grid_limits()` and `imr_tune()`. It returns the tuning result
  together with the data object and the grid. When the grid fixes every tuning
  parameter, the model is fitted once, without a validation split.
* New validation metric `"medae"`, the median absolute error, in
  `imr_tune_grid()`, `get_metric()` and `evaluate()`. It is computed only when
  requested by name: the default table of `evaluate()` is unchanged.
* `evaluate()` gains `mask`: only the entries where `mask == 1` are evaluated.
  When `actual` is an incomplete matrix (see `as_incomplete()`), only its
  observed entries are evaluated; previously its unobserved entries were
  compared as zeros.
* New `imr_identify()`: an optional step after `imr_fit()` that returns the
  decomposition of the same fitted matrix in which the low-rank component is
  orthogonal to the row and column covariates (and to the intercepts), and the
  column-covariate effects are orthogonal to the row covariates. The result is
  an `imr_fit` object, so `coef()`, `summary()`, `reconstruct()` and
  `reconstruct_partial()` apply to it. Models with shared effects are not
  supported.
* `imr_data()` gains `center` and `scale` (both `TRUE`). The covariates are
  centered and scaled to unit standard deviation instead of being orthonormalized
  by QR, and the covariate blocks are updated by a proximal-gradient step with
  step size $1/\|X\|_{op}^2$ (resp. $1/\|Z\|_{op}^2$). One copy of `X` and `Z` is
  stored, with their column centers and scales. `update()` gains `center` and
  `scale` to change these afterwards. The Lasso penalties apply to the
  standardized coefficients; `coef()`, `summary()` and `reconstruct()` report the
  coefficients on the original scale. This changes default results.
* `imr_data()` now stops when a row or column of `Y` has fewer than two observed
  entries, or when the training split leaves one without any.
* `imr_fit()` gains `return_objective`. When `TRUE`, `fit$objective` holds the
  value of the objective after each iteration.
* Fits without a low-rank component no longer stop after one iteration, and the
  stopping rule now includes the intercepts. This changes the least-squares
  initialization and the `lambda_m` grid maximum from `imr_set_grid_limits()`.
* `imr_set_grid_limits()` now verifies a bound 5% above its KKT value, so a
  correct bound is no longer doubled.
* `imr_fit()` no longer iterates when `lambda_m` is at or above its critical
  value, the largest singular value of the residuals of the model without the
  low-rank term. With the default start (`ls_initial = TRUE`) it returns that
  model, with `d = 0` and zero iterations. Such a fit used to need about
  `thresh^(-1/2)` iterations, or to stop at `maxit` with "Did not converge",
  and it returned a small nonzero `d`. In `imr_tune()`, the fit at the
  `lambda_m` grid maximum therefore has rank 0, which can change the selected
  model. Models with similarity matrices are not affected. A warm start whose
  `d` is all zero is now replaced by the default start without a warning.
* `imr_tune()` now passes `nuclear_log_scale` on when it also tunes
  `lambda_beta` or `lambda_gamma`.
* Clearer error message when an internal eigendecomposition fails.
* `imr_convergence(ls_initial = FALSE)` was documented and printed as a random
  initialization. It is not random: the coefficients and intercepts start at
  zero and the low-rank component at the truncated SVD of `Y`. The help page
  and `print()` now say so; the fits are unchanged.

# IMR 0.1.2
I added a 'src/Makevars' file to link the BLAS, LAPACK, and Fortran libraries.

# IMR 0.1.1

* Initial Submission to CRAN.

# IMR 1.0.0

* Article explaining the framework is uploaded to arXiv and is cited in (`inst/CITATION`).
* Improved console output of the `print()` and `summary()` methods.
* Added normalization to the similarity matrices as described in the Appendix of the article.
