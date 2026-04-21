#' Bixi sample data
#' A subset of data from Bixi (https://bixi.com), a docked bike-sharing service
#' in Montreal, Canada.
#' We use the data compiled by Lei et al., 2025, which contains normalized daily
#' departure counts for each of 579 stations over 193 days (April 15 to
#'  October 27, 2019). The data is accompanied by side information.
#'  We select two temporal variables (temperature and precipitation) and
#' one spatial variable (park area). We also take a sample of 100 stations and
#' the first 100 consecutive days. The data also contains two distance matrices:
#' a spatial distance matrix derived from the geographic coordinates of the stations
#'  (columns) and a temporal distance matrix representing the elapsed days (rows)
#'
#' @format ## `bixi_sample`
#' A list of four matrices:
#' \describe{
#'    \item{Y}{a 100 x 100 matrix of normalized departure counts. Missing values
#'        are set to zero. Each row is a day and each column is a station.}
#'    \item{test}{a 100 x 100 test matrix where test indices are nonzero. All
#'        nonzero elements in `test` are zero in `Y` and vice-versa.}
#'    \item{X}{a 100 x 2 matrix of two row covariates: temperature and precipitation.
#'        Both are normalized.}
#'    \item{Z}{a 100 x 1 matrix of one column covariate: park area.}
#'    \item{spatial_distance}{a 100 x 100 distance matrix between stations}
#'    \item{temporal_distance}{a 100 x 100 distance matrix for days}
#' }
#'
#' @references
#' Lei, M., Labbe, A., & Sun, L. (2025). Scalable Spatiotemporally Varying
#' Coefficient Modelling with Bayesian Kernelized Tensor Regression.
#' \emph{Bayesian Analysis}, 20(3). \doi{10.1214/24-BA1428}.
"bixi_sample"
