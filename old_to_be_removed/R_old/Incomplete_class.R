#' @export
#' @importClassesFrom softImpute Incomplete
#' @export
as_incomplete <- function(x, ...) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "CsparseMatrix")
  na <- is.na(x@x)
  if(any(na)){
    x@x[na] <- 0
    x <- Matrix::drop0(x)
  }
  x
}
#' @export
is_incomplete <- function(x, ...) inherits(x, "dgCMatrix")
#' @export
setMethod("as.matrix", "imr_incomplete", as_incomplete)
setAs("matrix", "imr_incomplete", function(from) as_incomplete(from))
