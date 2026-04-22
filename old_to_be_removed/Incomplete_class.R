#' @export
#' @importClassesFrom Matrix dgCMatrix
setClass("imr_incomplete", "dgCMatrix")

#' @export
as_incomplete <- function(x) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "dgCMatrix")
  na <- is.na(x@x)
  if (any(na)) {
    x@x[na] <- 0
    x <- Matrix::drop0(x)
  }
  new("imr_incomplete", x)
}

#' @export
is_incomplete <- function(x) inherits(x, "imr_incomplete")

#' @export
setAs("matrix", "imr_incomplete", function(from) as_incomplete(from))

#' @export
setAs("Matrix", "imr_incomplete", function(from) as_incomplete(from))



#' @export
#' @importMethodsFrom Matrix t
t.Incomplete <- function(x) {
  as_incomplete(t(as(x, "dgCMatrix")))
}




# #' @export
# #' @method as.matrix Incomplete
# as.matrix.Incomplete <- function(x, ...) {
#   out <- as.matrix(as(x, "dgCMatrix"))
#   return(out)
# }
# #' @export
# setMethod("as.matrix", "imr_incomplete", as.matrix.Incomplete)
