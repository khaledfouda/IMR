
# IMR

<!-- badges: start -->
[![R-CMD-check](https://github.com/khaledfouda/IMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/khaledfouda/IMR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Citation

If you use IMR, please cite:

> Fouda, K., Labbe, A., & Oualkacha, K. (2026). *Incomplete Matrix Regression*. arXiv preprint arXiv:2606.26325.


The IMR package provides a framework for both matrix completion and regression on response matrices with missing values. Let $\boldsymbol{Y} \in \Re^{n\times m}$ denote the observed incomplete matrix, where missing values are designated by either NA or zero. The estimator for each entry of the matrix is specified by any combination of the following

$$\boldsymbol{\hat Y_{ij}} = \boldsymbol{\hat\beta_{oi}} + \boldsymbol{\hat\Gamma_{oj}} + \boldsymbol{X_i\hat\beta_j} + \boldsymbol{\hat\Gamma_i Z^{'}_j} + \boldsymbol{\hat M_{ij}}$$

where $\boldsymbol{X}\in \Re^{n\times p}$ and  $\boldsymbol{Z}\in \Re^{m\times q}$ are row ($p$ predictors) and column ($q$ predictors) covariate matrices, respectively. The vectors $\boldsymbol{\hat\beta_{o}}\in \Re^{n}$ and $\boldsymbol{\hat\Gamma_{o}}\in \Re^{m}$ represent the row-level and column-level intercepts. The term $\boldsymbol{\hat\beta}$ denotes the row covariate coefficients, which may be structured as either an $p\times m$ matrix (one coefficient of each (predictor,column) pair, or a $p$-dimensional vector (one coefficient for each predictor), forcing coefficients to be equal across all columns for each covariate. Similarly, the column covariate coefficients, denoted by $\boldsymbol{\hat\Gamma}$, can be either an $n \times q$ matrix or a $q$-dimensional vector where all rows share the same coefficient for each covariate. To avoid having too many parameters, we impose Lasso ($L_1$) penalties on the row and column covariate coefficients. Finally, $\boldsymbol{M}$ is a rank-r (r is a hyper-parameter) low-rank matrix subject to a nuclear norm penalty. Together, they yield the following penalty structure:

$$\mathrm{Penalty} = \lambda_\beta {\|\boldsymbol{\beta}\|}_{1} + \lambda_\Gamma {\|\boldsymbol{\Gamma}\|}_{1} + \lambda_m {\|\boldsymbol{S_r^{1/2}MS_c^{1/2}}\|}_{*},$$

where $\boldsymbol{S_r} \in \Re^{n\times n}$ and $\boldsymbol{S_c} \in \Re^{m \times m}$ are similarity (or information) matrices that describe the correlation structure among the rows and columns of the response matrix, respectively. In the absence of a known correlation structure, this penalty term reduces to the standard nuclear norm, ${\|\boldsymbol{M}\|}_{*}$. We have 4 penalty parameters: $(\lambda_\beta,  \lambda_\Gamma, \lambda_m, r)$. We provide a method to estimate those parameters. 

As we said above, we can use any combination of the model components to define our estimator. Examples of these combinations include:

$\boldsymbol{\hat Y_{ij}} = \boldsymbol{\hat\beta_{oi}} + \boldsymbol{X_i\hat\beta}$ (where $\boldsymbol{\hat\beta}$ is a p-dimensional vector).

$\boldsymbol{\hat Y_{ij}} = \boldsymbol{\hat\beta_{oi}} + \boldsymbol{\hat\Gamma_{oj}} +\boldsymbol{\hat M_{ij}}$ with $\boldsymbol{S_r}$ set to a the inverse of a Matern kernel and $S_c$ left unspecified.

$\boldsymbol{\hat Y_{ij}} = \boldsymbol{\hat M_{ij}}$, where neither $\boldsymbol{S_r}$ or $\boldsymbol{S_c}$ is specified. This corresponds exactly to the Soft-Impute model [@hastie2015].




## Installation

You can install the development version of IMR from [GitHub](https://github.com/) with:


``` r
# install.packages("remotes")
remotes::install_github("khaledfouda/IMR", build_vignettes = TRUE)
# or
# install.packages("pak")
# pak::pak("khaledfouda/IMR")
```

## Example

To illustrate the standard workflow within the package, assume the matrices $\boldsymbol{Y}$ and $\boldsymbol{X}$ are defined as above, and the objective is to fit example 1 above to obtain the complete estimated matrix $\boldsymbol{\hat Y}$. Then, the implementation is as follows:
``` r
library(IMR)
# set the hyperparameter value.
lambda_beta <- 0.02
# load the data example (see ?IMR::Bixi_sample for more information)
Bixi <- IMR::Bixi_sample
# create the data object
data <- imr_data(Y = Bixi$Y, X = Bixi$X)
# update the model structure to fit example 1
data <- update(data, row_covariates = TRUE, # turn  XBeta on (on by default when X is provided)
               shared_beta = TRUE, # make beta a p-dimensional vector (off by default)
               low_rank_component = FALSE, # turn M off (on by default)
               row_intercept = TRUE) # turn row intercepts on (off by default).
# fit the model
fit  <- imr_fit(data, lambda_beta = lambda_beta )
# obtain \hat{Y}
Y_hat <- reconstruct(fit, data)$estimates

```

For a detailed overview of the package and step-by-step example, 
please read the vignette. You can access it from within R by running:

``` r
vignette("IMR", "IMR")
```


