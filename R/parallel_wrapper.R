#' Run a 2-parameter grid in parallel
#'
#' Evaluate a function over all combinations of two tuning parameters
#' (`lambda_beta` and `lambda_gamma`) using the future ecosystem for parallel
#' execution.
#'
#' This helper is currently specialised to grids over `lambda_beta` and
#' `lambda_gamma`.
#'
#' @param grid Named list of length 2 with elements `lambda_beta` and
#'   `lambda_gamma`. Each element should be an atomic vector of candidate values.
#' @param f Function to be evaluated on each parameter combination. It must
#'   accept arguments named `lambda_beta` and `lambda_gamma`, along with any
#'   additional arguments passed through `...`.
#' @param combine How to combine the individual results. `"list"` (default)
#'   returns a list; `"rbind"` row-binds data frame results and adds the
#'   tuning-parameter columns.
#' @param .seed Value passed to `future.apply::future_lapply()` as
#'   `future.seed`. Set to `TRUE` for reproducible parallel RNG streams.
#' @param .packages Character vector of package names to attach on workers;
#'   passed to `future.apply::future_lapply()` as `future.packages`.
#' @param .progress Logical. If `TRUE` and the {progressr} package is installed,
#'   a progress bar is shown while evaluating the grid.
#' @param ... Additional arguments forwarded to `f`.
#'
#' @details
#' The function constructs the full Cartesian product of the supplied
#' `lambda_beta` and `lambda_gamma` values, evaluates `f()` for each
#' combination in parallel via `future.apply::future_lapply()`, and then
#' combines the results.
#'
#' When `combine = "list"`, the result is a list of length equal to the number
#' of parameter combinations. The original grid of combinations is attached as
#' the `"grid"` attribute.
#'
#' When `combine = "rbind"`, every call to `f()` must return a data frame. Two
#' extra columns corresponding to the grid parameters are added to each result
#' before row-binding into a single data frame.
#'
#' The function requires the {future.apply} package, and optionally {progressr}
#' when `.progress = TRUE`.
#'
#' @return
#' A list or data frame, depending on `combine`:
#'
#' * If `combine = "list"`, a list of results, with an attached `"grid"`
#'   attribute containing the parameter combinations.
#' * If `combine = "rbind"`, a single data frame created by row-binding all
#'   results and adding the tuning-parameter columns.
#'
#' @examples
#' \dontrun{
#' grid <- list(
#'   lambda_beta  = c(0.1, 1),
#'   lambda_gamma = c(0, 0.5)
#' )
#'
#' f <- function(lambda_beta, lambda_gamma) {
#'   data.frame(
#'     lambda_beta  = lambda_beta,
#'     lambda_gamma = lambda_gamma,
#'     score        = lambda_beta + lambda_gamma
#'   )
#' }
#'
#' # Return a data frame with parameter columns attached
#' parallel_grid(grid, f, combine = "rbind")
#' }
#'
#' @export
parallel_grid <- function(grid,
                          f,
                          combine = "list",
                          .seed = TRUE,
                          .packages = NULL,
                          .progress = FALSE,
                          ...) {
  combine <- match.arg(combine, c("list", "rbind"))

  if (!is.list(grid) || length(grid) != 2L || any(!nzchar(names(grid)))) {
    stop("`grid` must be a named list of length 2 as list(lambda_beta = ..., lambda_gamma = ...).")
  }

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future.apply' to enable parallel execution.")
  }

  stopifnot(names(grid) %in% c("lambda_beta", "lambda_gamma"))

  # Build full Cartesian product of the two parameters
  combos <- do.call(
    expand.grid,
    c(grid, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  )

  n_tasks <- nrow(combos)
  stopifnot(n_tasks > 0L)

  # Turn each row into a named list of arguments for f()
  tasks <- lapply(
    X = seq_len(n_tasks),
    FUN = function(i) as.list(combos[i, , drop = FALSE])
  )

  param_names <- names(grid)

  # Worker wrapper: evaluate f(args, ...) and optionally report progress via `p`
  worker <- function(args) {
    on.exit(if (!is.null(p)) p(), add = TRUE)
    do.call(f, c(args, list(...)))
  }

  # Helper: run the grid in parallel without progress reporting
  run <- function() {
    future.apply::future_lapply(
      X = tasks,
      FUN = worker,
      future.seed = .seed,
      future.packages = .packages
    )
  }

  # Evaluate the grid, optionally with a progress bar
  res <- if (.progress && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(along = tasks)
      future.apply::future_lapply(
        X = tasks,
        FUN = worker,
        future.seed = .seed,
        future.packages = .packages
      )
    })
  } else {
    run()
  }

  if (combine == "list") {
    attr(res, "grid") <- combos
    return(res)
  }

  # rbind path: ensure each result is a data frame, and add parameter columns
  are_df <- vapply(res, is.data.frame, logical(1))
  if (!all(are_df)) {
    bad <- which(!are_df)[1]
    stop(sprintf(
      "combine = 'rbind' requires data.frame results; element %d is a %s.",
      bad, class(res[[bad]])[1]
    ))
  }

  res_with_ids <- Map(
    f = function(df, args) {
      df[[param_names[1]]] <- args[[param_names[1]]]
      df[[param_names[2]]] <- args[[param_names[2]]]
      df[c(param_names, setdiff(names(df), param_names))]
    },
    res,
    tasks
  )

  do.call(rbind, res_with_ids)
}
# ================================================================================
#' Parallel search over a 1D parameter grid
#'
#' Evaluate a function over a one-dimensional grid of parameter values using the
#' future ecosystem for parallel execution. The function `f` is expected to take
#' a single argument `parameter` (the grid value) plus any additional arguments
#' passed through `...`, and return a list with elements `error` and `fit`.
#'
#' @param grid Atomic vector of candidate values for the tuning parameter.
#' @param f Function to be evaluated on each parameter value. It must have the
#'   signature `f(parameter, ...)` and return a list with components:
#'   \describe{
#'     \item{error}{Numeric scalar giving the validation error for that fit.}
#'     \item{fit}{An arbitrary R object (typically a model fit).}
#'   }
#' @param .seed Value passed to `future.apply::future_lapply()` as
#'   `future.seed`. Set to `TRUE` for reproducible parallel RNG streams.
#' @param .packages Character vector of package names to attach on workers;
#'   passed to `future.apply::future_lapply()` as `future.packages`.
#' @param .progress Logical. If `TRUE` and the {progressr} package is installed,
#'   a progress bar is shown while evaluating the grid.
#' @param ... Additional arguments forwarded to `f`.
#'
#' @details
#' For each value in `grid`, the function calls:
#'
#' \preformatted{
#'   f(parameter = grid[i], ...)
#' }
#'
#' and expects a list with elements named `error` and `fit`. The `error`
#' values are collected into a data frame alongside the corresponding
#' parameter values, and the parameter with the smallest error is selected as
#' the best.
#'
#' Parallelisation is handled via `future.apply::future_lapply()`. A progress
#' bar is optionally shown using {progressr} when `.progress = TRUE`.
#'
#' @return
#' A list with components:
#'
#' \describe{
#'   \item{results}{A data frame with columns `parameter` and `error`, giving
#'     the grid and corresponding validation errors.}
#'   \item{best_parameter}{The parameter value with the smallest error.}
#'   \item{best_error}{The smallest error value.}
#'   \item{best_fit}{The `fit` object corresponding to `best_parameter`.}
#' }
#'
#' @examples
#' \dontrun{
#' grid <- seq(0, 1, length.out = 5)
#'
#' f <- function(parameter, n = 100) {
#'   # Dummy example: quadratic loss around 0.4
#'   err <- (parameter - 0.4)^2 + rnorm(1, sd = 0.01)
#'   fit <- list(param = parameter, n = n)
#'   list(error = err, fit = fit)
#' }
#'
#' res <- parallel_grid_1d(
#'   grid = grid,
#'   f = f,
#'   .progress = FALSE
#' )
#'
#' res$results
#' res$best_parameter
#' res$best_error
#' str(res$best_fit)
#' }
#'
#' @export
parallel_grid_1d <- function(grid,
                             f,
                             .seed = TRUE,
                             .packages = NULL,
                             .progress = FALSE,
                             ...) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future.apply' to enable parallel execution.")
  }

  if (!is.atomic(grid) || length(grid) == 0L) {
    stop("`grid` must be a non-empty atomic vector of candidate values.")
  }

  parameter_values <- grid
  n_tasks <- length(parameter_values)

  # Progressor placeholder, used only when .progress = TRUE
  p <- NULL

  # Worker: evaluate f(parameter = value, ...) and optionally update progress
  worker <- function(i) {
    if (!is.null(p)) {
      on.exit(p(), add = TRUE)
    }

    value <- parameter_values[[i]]
    res <- f(value, ...)

    if (!is.list(res) || !all(c("error", "fit") %in% names(res))) {
      stop("`f` must return a list with elements named 'error' and 'fit'.")
    }

    res
  }

  # Helper: run the grid in parallel
  run <- function() {
    future.apply::future_lapply(
      X = seq_len(n_tasks),
      FUN = worker,
      future.seed = .seed,
      future.packages = .packages
    )
  }

  # Evaluate the grid, optionally with a progress bar
  res <- if (.progress && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p <<- progressr::progressor(along = seq_len(n_tasks))
      run()
    })
  } else {
    run()
  }

  # Extract errors and fits
  errors <- vapply(
    X = res,
    FUN = function(x) x$error,
    FUN.VALUE = numeric(1)
  )

  if (all(is.na(errors))) {
    stop("All `error` values returned by `f` are NA; cannot select a best parameter.")
  }

  fits <- lapply(res, function(x) x$fit)

  # Build results data frame
  results <- data.frame(
    parameter = parameter_values,
    error = errors,
    stringsAsFactors = FALSE
  )

  # Select best parameter (smallest error)
  best_idx <- which.min(errors)
  best_parameter <- parameter_values[best_idx]
  best_error <- errors[best_idx]
  best_fit <- fits[[best_idx]]

  list(
    results = results,
    best_parameter = best_parameter,
    best_error = best_error,
    best_fit = best_fit
  )
}


#' @export
parallel_grid_1d_adaptive <- function(param_min,
                                      param_max,
                                      step_sizes,
                                      f,
                                      .seed = TRUE,
                                      .packages = NULL,
                                      .progress = FALSE,
                                      .trace = FALSE,
                                      ...) {
  # 1. check that step_sizes is numeric, positive, and decreasing
  if (!is.numeric(step_sizes) || length(step_sizes) == 0L) {
    stop("`step_sizes` must be a non-empty numeric vector or list of numeric scalars.")
  }

  if (any(step_sizes <= 0)) {
    stop("All `step_sizes` must be strictly positive.")
  }

  # We expect coarse-to-fine: e.g., 0.1, 0.01, 0.001 (decreasing)
  if (!isTRUE(all(diff(step_sizes) < 0))) {
    warning(
      "`step_sizes` is not strictly decreasing; ",
      "for coarse-to-fine search, use larger steps first, then smaller."
    )
  }


  original_min <- param_min
  original_max <- param_max

  for (step in step_sizes) {

    grid <- seq(param_min, param_max, step)
    if (length(grid) == 0L) {
      stop(
        "No grid points generated for step size ", step,
        " with current bounds [", param_min, ", ", param_max, "]."
      )
    }

    result <- parallel_grid_1d(grid, f, .seed, .packages, .progress, ...)

    # adjust min and max for next grid. they should go from  param+step <-> param-step
    param_min <- max(original_min, result$best_parameter - step)
    param_max <- min(original_max, result$best_parameter + step)

    if (.trace) {
      message(
        sprintf(
          "step = %g, best_parameter = %.6f, best_error = %.6f",
          step,
          result$best_parameter,
          result$best_error
        )
      )
    }
  }

  return(result)
}


