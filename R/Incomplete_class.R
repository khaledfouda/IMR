#' @export
#' @importClassesFrom softImpute Incomplete
as.Incomplete <- function(x) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "CsparseMatrix")
  na <- is.na(x@x)
  if (any(na)) {
    x@x[na] <- 0
    x <- Matrix::drop0(x)
  }
  methods::new("Incomplete", x)
}

#' @export
is.Incomplete <- function(x) inherits(x, "Incomplete")

#' @export
setAs("matrix", "Incomplete", function(from) as.Incomplete(from))

#' @export
setAs("Matrix", "Incomplete", function(from) as.Incomplete(from))

#' @export
#' @exportS3Method base::as.matrix
as.matrix.Incomplete <- function(x, ...) {
  # Unwrap and force into a standard base R dense matrix
  base::as.matrix(methods::as(x, "CsparseMatrix"))
}
#' @export
#' @exportS3Method base::t
t.Incomplete <- function(x) {
  # Unwrap, transpose, and re-wrap as Incomplete
  as.Incomplete(t(methods::as(x, "CsparseMatrix")))
}
