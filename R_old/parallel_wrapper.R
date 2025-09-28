#' @export
parallel_grid <- function(grid,
                          f,
                          combine = "list",
                          # id_prefix = ".param_",
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
  tasks <- lapply(seq_len(n_tasks), function(i) as.list(combos[i, , drop = FALSE]))
  param_names <- names(grid)
  #id_cols <- paste0(id_prefix, param_names)

  # Worker wrapper: evaluate f(args, ...) and progress
  worker <- function(args) {
    on.exit(if (!is.null(p)) p(), add = TRUE)
    do.call(f, c(args, list(...)))
  }

  # progress bar
  run <- function() {
    future.apply::future_lapply(
      X = tasks,
      FUN = worker,
      future.seed = .seed,
      future.packages = .packages
    )
  }

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

  # rbind path: ensure each result is a data frame, and add ID columns
  are_df <- vapply(res, is.data.frame, logical(1))
  if (!all(are_df)) {
    bad <- which(!are_df)[1]
    stop(sprintf("combine = 'rbind' requires data.frame results; element %d is a %s.",
                 bad, class(res[[bad]])[1]))
  }

  res_with_ids <- Map(function(df, args) {
    # Avoid clobbering user columns by using a prefix
    df[[param_names[1]]] <- args[[param_names[1]]]
    df[[param_names[2]]] <- args[[param_names[2]]]
    # Put IDs first
    df[c(param_names, setdiff(names(df), param_names))]
  }, res, tasks)

  # if (requireNamespace("data.table", quietly = TRUE)) {
  #   return(data.table::rbindlist(res_with_ids, use.names = TRUE, fill = TRUE))
  # }
  # if (requireNamespace("dplyr", quietly = TRUE)) {
  #   return(dplyr::bind_rows(res_with_ids))
  # }

  do.call(rbind, res_with_ids)
}
