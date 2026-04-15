#' @export
#' @importClassesFrom Matrix dgCMatrix
setClass("Incomplete", "dgCMatrix")

#' @export
as.Incomplete <- function(x) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "dgCMatrix")
  na <- is.na(x@x)
  if (any(na)) {
    x@x[na] <- 0
    x <- Matrix::drop0(x)
  }
  new("Incomplete", x)
}

#' @export
is.Incomplete <- function(x) inherits(x, "Incomplete")

#' @export
setAs("matrix", "Incomplete", function(from) as.Incomplete(from))

#' @export
setAs("Matrix", "Incomplete", function(from) as.Incomplete(from))



#' @export
#' @method as.matrix Incomplete
as.matrix.Incomplete <- function(x, ...) {
  out <- as.matrix(as(x, "dgCMatrix"))
  return(out)
}

#' @export
setMethod("as.matrix", "Incomplete", as.matrix.Incomplete)

#' @export
#' @importMethodsFrom Matrix t
t.Incomplete <- function(x) {
  as.Incomplete(t(as(x, "dgCMatrix")))
}
