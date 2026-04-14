#' @export
#' @importClassesFrom softImpute Incomplete
#' @export
as.Incomplete <- function(x, ...) {
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
is.Incomplete <- function(x, ...) inherits(x, "dgCMatrix")
#' @export
setMethod("as.matrix", "Incomplete", as.Incomplete)
setAs("matrix", "Incomplete", function(from) as.Incomplete(from))
