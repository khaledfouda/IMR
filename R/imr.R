#' @title Prepare the data, tune and fit an IMR model in one call
#'
#' @description
#' `imr()` runs the usual sequence of calls with their most used arguments:
#' [imr_data()], [update.imr_data()], [imr_set_grid_limits()] and [imr_tune()].
#' For the other arguments of these functions, call them one by one.
#'
#' @inheritParams imr_data
#' @inheritParams imr_convergence
#' @inheritParams imr_tune
#' @param row_intercept,col_intercept Logical. Include the row intercepts and
#'   the column intercepts in the model? Both default to `FALSE`.
#' @param val_prop Numeric, strictly between 0 and 1. The proportion of the
#'   observed entries of `Y` held out for validation. Defaults to `0.2`.
#' @param grid An object of class `"imr_tune_grid"`. The maxima left at `NA`
#'   are computed by [imr_set_grid_limits()]. Defaults to `imr_tune_grid()`.
#' @param metric Character. The error metric computed on the validation set;
#'   see [imr_tune_grid()]. Defaults to the metric of `grid`.
#' @param seed Integer. An optional random seed, passed to [imr_data()] and to
#'   [imr_tune()]. Defaults to `NULL`.
#' @param verbose Integer. `0` is silent. `1` (the default) reports each step
#'   and the time it took. The functions called receive `verbose - 1`, so
#'   larger values also show their output. Above `10`, every fit prints its
#'   iterations.
#'
#' @details
#' The steps are: the data object with its validation split ([imr_data()]),
#' the intercepts ([update.imr_data()]), the maxima of the grid
#' ([imr_set_grid_limits()], on the training split), and the tuning
#' ([imr_tune()]), which ends with the fit on all the observed entries when
#' `final_fit = TRUE`. Every fit uses
#' `imr_convergence(maxit = maxit, thresh = thresh)`.
#'
#' Nothing is tuned when `grid` fixes every tuning parameter of the model, that
#' is, when the minimum equals the maximum for the rank, for \eqn{\lambda_M},
#' and for \eqn{\lambda_\beta} and \eqn{\lambda_\Gamma} when `X` and `Z` are
#' given. No validation split is made then, and the model is fitted once on all
#' the observed entries by [imr_fit()]. The arguments `val_prop`, `metric`,
#' `fast_nuclear`, `nuclear_log_scale`, `final_fit` and `n_cores` are not used.
#'
#' @return A list. When the model is tuned, it is the list returned by
#'   [imr_tune()], whose element `fit` is the `"imr_fit"` object, with two more
#'   elements:
#' \itemize{
#'   \item \code{data}: the `"imr_data"` object the model was fitted on, which
#'     [reconstruct()], [reconstruct_partial()] and [imr_identify()] require.
#'   \item \code{grid}: the grid, with its computed maxima.
#' }
#'   When nothing is tuned, the list holds `fit`, `data` and `grid` only.
#'
#' @seealso [imr_data()], [imr_tune_grid()], [imr_set_grid_limits()],
#'   [imr_tune()], [imr_fit()]
#'
#' @examples
#' # three observed entries in every row and column, so that each keeps at
#' # least one in the training split
#' Y <- matrix(
#'   c(NA, 1, 3, 4,
#'     4, NA, 2, 4,
#'     2, 1, NA, 3,
#'     5, 2, 4, NA), 4, byrow = TRUE
#' )
#'
#' # tune lambda_m and the rank of the model Y = M, then fit on all the entries
#' grid <- imr_tune_grid(nuclear = c(0, NA, 5, 2), rank = c(2, 5, 1, 2))
#' out <- imr(Y, grid = grid, n_cores = 1, seed = 1)
#' out$params
#'
#' # the data object is returned for the functions that need it
#' estimates <- reconstruct(out$fit, out$data, trace = FALSE)$estimates
#'
#' # a grid that fixes both parameters: the model is fitted once
#' out <- imr(Y, grid = imr_tune_grid(nuclear = 0.1, rank = 2))
#'
#' @export
imr <- function(Y,
                X = NULL,
                Z = NULL,
                row_intercept = FALSE,
                col_intercept = FALSE,
                val_prop = 0.2,
                grid = imr_tune_grid(),
                metric = grid$metric,
                maxit = 600,
                thresh = 1e-5,
                fast_nuclear = TRUE,
                nuclear_log_scale = TRUE,
                final_fit = TRUE,
                n_cores = 4,
                seed = NULL,
                verbose = 1) {
  stopifnot(inherits(grid, "imr_tune_grid"))
  # nothing to tune when every parameter of the model is fixed in the grid
  fixed <- function(p) identical(p$min, p$max)
  tune <- !(fixed(grid$nuclear) && fixed(grid$rank) &&
              (is.null(X) || fixed(grid$beta)) &&
              (is.null(Z) || fixed(grid$gamma)))
  if (tune) {
    .imr_check_param(val_prop, "numeric", 0, 1, min_inclusive = FALSE,
                     max_inclusive = FALSE)
    get_metric(metric) # stops on an unknown metric
    grid$metric <- toupper(metric)
  }
  convergence <- imr_convergence(maxit = maxit, thresh = thresh,
                                 trace = verbose > 10, ls_initial = TRUE)
  # runs one step; announces it and reports its time when verbose > 0
  step <- function(label, expr) {
    if (verbose <= 0) return(expr)
    message("imr: ", label)
    start <- Sys.time()
    value <- expr
    message(sprintf("imr: done in %.2f s",
                    as.numeric(Sys.time() - start, units = "secs")))
    return(value)
  }

  data <- step("preparing the data", update.imr_data(
    imr_data(Y, X, Z, val_prop = if (tune) val_prop else 0, seed = seed),
    row_intercept = row_intercept, col_intercept = col_intercept
  ))
  if (!tune) {
    fit <- step("fitting the model", imr_fit(
      data, rank = grid$rank$min, lambda_m = grid$nuclear$min,
      lambda_beta = grid$beta$min, lambda_gamma = grid$gamma$min,
      convergence = convergence
    ))
    return(list(fit = fit, data = data, grid = grid))
  }
  grid <- step("computing the grid limits", imr_set_grid_limits(
    data, grid, convergence = convergence, training = TRUE,
    verbose = verbose - 1
  ))
  out <- step("tuning", imr_tune(
    data, grid, final_fit = final_fit, fast_nuclear = fast_nuclear,
    convergence = convergence, verbose = verbose - 1, n_cores = n_cores,
    seed = seed, nuclear_log_scale = nuclear_log_scale
  ))
  out$data <- data
  out$grid <- grid
  return(out)
}
