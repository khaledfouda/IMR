#' @export
#' @importClassesFrom softImpute Incomplete
as.matrix.Incomplete <- function(x, ...) {
  stopifnot(is.matrix(x))
  x[x == 0] <- NA
  x <- as(x, "Incomplete")
  x
}
setMethod("as.matrix", "Incomplete", as.matrix.Incomplete)
