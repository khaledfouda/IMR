#' @export
#' @importClassesFrom softImpute Incomplete
as.Incomplete <- function(x) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "dgCMatrix")
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
setMethod("as.matrix", "Incomplete", as.Incomplete)

#' @export
setAs("matrix", "Incomplete", function(from) as.Incomplete(from))

#' @export
setAs("Matrix", "Incomplete", function(from) as.Incomplete(from))

#  to allow converting an Incomplete back to a dense matrix
#' @export
setMethod("as.matrix", "Incomplete", function(x) {
  as.matrix(as(x, "dgCMatrix"))
})
