#' @title Clean and convert any matrix to a sparse object
#' @return A standard `dgCMatrix` with both NAs and 0s are treated as missing values.
#' @export
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix drop0
#' @importMethodsFrom Matrix t
#' @examples
#' # create sample data
#' Y <- matrix(
#'   c(2, NA, 0,
#'   4, .5, NA,
#'   NA, NA, 0), 3, byrow= TRUE
#' )
#' # make it sparse with both NAs and 0s dropped
#' Y <- as_incomplete(Y)
#' # verify
#' print(is_incomplete(Y))
as_incomplete <- function(x) {
  stopifnot(inherits(x, c("matrix", "Matrix")))
  x <- as(x, "dgCMatrix")
  # replace NAs with zeros then drop them.
  na <- is.na(x@x)
  if (any(na)) {
    x@x[na] <- 0
    x <- Matrix::drop0(x)
  }
  # add a check attributed
  attr(x, "imr_cleaned") <- TRUE
  return(x)
}

#' @title Has the matrix be processed by `as_incomplete()`?
#' @return Boolean
#' @export
#' @seealso [as_incomplete()]
#' @inherit as_incomplete examples
is_incomplete <- function(x) isTRUE(attr(x, "imr_cleaned"))
