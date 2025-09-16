#' @export
naive_MC <- function(mat) {
  # Input assumes that empty cells are filled with NA so they don't affect the average
  #mat[mat == 0] <- NA
  # Calculate row and column means, excluding NAs
  if(! IMR::is.Incomplete(mat)) mat <- IMR::as.Incomplete(mat)
  # have this instead and accept sparse matrix.
  row_means <- IMR:::row_means_cpp(mat, ncol(mat))
  col_means <- IMR:::col_means_cpp(mat, nrow(mat))
  # row_means <- rowMeans(mat, na.rm = TRUE)
  # col_means <- colMeans(mat, na.rm = TRUE)

  # in the case that some rows/columns don't have obeservations we fill them with overall mean
  # if (any(is.na(col_means)) | any(is.na(row_means))) {
  #   total_mean <- mean(c(row_means, col_means), na.rm = TRUE)
  #   col_means <- tidyr::replace_na(col_means, total_mean)
  #   row_means <- tidyr::replace_na(row_means, total_mean)
  # }


  # Expand row and column means to match the matrix dimensions
  # expanded_row_means <-
  #   matrix(
  #     row_means[rep(1:length(row_means), each = ncol(mat))],
  #     nrow = nrow(mat),
  #     ncol = ncol(mat),
  #     byrow = TRUE
  #   )
  # expanded_col_means <-
  #   matrix(col_means[rep(1:length(col_means), times = nrow(mat))], nrow = nrow(mat), ncol = ncol(mat))
  # avg_means <- (expanded_row_means + expanded_col_means) / 2
  # mat[na_positions] <- avg_means[na_positions]

   na_positions <- mat == 0
  ij <- Matrix::which(na_positions, arr.ind=TRUE)
  mat[ij] <- (row_means[ij[,1]] + col_means[ij[, 2]]) / 2
  mat <- Matrix::drop0(mat)

  return(mat)
}
