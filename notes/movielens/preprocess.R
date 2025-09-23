
prepare_ml_1m_data <- function(min_obs_per_col = 10,
                               increase_missing = FALSE,
                               prop_miss = .98,
                               seed = 2025){

  if(is.numeric(seed)) set.seed(seed)
  require(Matrix)
  require(dplyr)
  require(magrittr)
  setwd("~/Research/CASMC/MovieLens/")
  load("data/ml-1m/Movie_Y.Rdata") # Y
  load("data/ml-1m/Movie_X.Rdata") # X
  load("data/ml-1m/Movie_Q.Rdata") # Query (testing set for evaluating the performance)
  query <- as.data.frame(query)
  colnames(query) <- c("row_id", "column_id", "value")
  #=====================================================================
  query_to_matrix <- function(query,
                              dims=c(max(query$row_id), max(query$column_id))){
    Matrix::sparseMatrix(
      i    = query$row_id,
      j    = query$column_id,
      x    = query$value,
      dims = dims
    )
  }

  matrix_to_query <- function(M){
    idx <- which(as.matrix(M!=0), arr.ind=TRUE)
    data.frame(row_id=idx[,1], column_id=idx[,2], value=M[idx]) %>%
      arrange(row_id, column_id)
  }
  summ_Y <- function(Y){
    message("Proportion of missing: ",round(mean(Y==0),3))
    col_counts <- colSums(Y!=0)
    row_counts <- rowSums(Y!=0)
    message(
      "Columns – ",
      "min obs: ", min(col_counts),
      "  max obs: ", max(col_counts),
      "\n",
      "Rows    – ",
      "min obs: ", min(row_counts),
      "  max obs: ", max(row_counts),
      "\n"
    )
  }
  #=====================================================================
  # I will remove columns with less than min_obs_per_col observations.
  query.mat <- query_to_matrix(query, c(dim(Y)[1],dim(Y)[2]))
  stopifnot(all(matrix_to_query(query.mat) == query))
  message("Dim of Y (before): [" , dim(Y)[1], ",", dim(Y)[2],"]")


  summ_Y(Y)
  message("Keeping only columns with at least ", min_obs_per_col,
          " observations.")
  columns_to_keep = base::which(colSums(Y!=0) >= min_obs_per_col)
  Y <- Y[, columns_to_keep]
  message("Dim of Y (after): [" , dim(Y)[1], ",", dim(Y)[2],"]")
  query.mat <- query.mat[, columns_to_keep]
  query <- matrix_to_query(query.mat)
  summ_Y(Y)
  keyword = paste0("_c_",min_obs_per_col,"_")
  message("Saving data with keyword: ", keyword)
  saveRDS(Y, paste0("data/ml-1m/Movie_Y",keyword,".Rdata"))
  saveRDS(query, paste0("data/ml-1m/Movie_Q",keyword,".Rdata"))
  #=========
  if(increase_missing)
  {

    message("Increase the proportion of missing data to ",prop_miss)

    # we increase percentage of missing to 98%
    n_rows     <- nrow(Y)
    n_cols     <- ncol(Y)
    total_cells <- n_rows * n_cols

    current_zero_count <- sum(Y == 0, na.rm = TRUE)
    target_zero_count  <- ceiling(prop_miss * total_cells)

    # how many new zeros we need
    n_to_add <- target_zero_count - current_zero_count
    nonzero_idx <- which(Y != 0, arr.ind = FALSE)

    if (length(nonzero_idx) < n_to_add) {
      stop("Not enough non-zero entries to reach 99% zeros.")
    }

    selected_idx <- sample(nonzero_idx, size = n_to_add)

    # 4. Decode to (row, col) and record original values
    #    arrayInd() turns linear indices back into row/col pairs
    rc_pairs <- arrayInd(selected_idx, .dim = dim(Y))

    # add the newly missing data into the test set
    query %<>% rbind(
      data.frame(
        row_id        = rc_pairs[, 1],
        column_id     = rc_pairs[, 2],
        value = Y[selected_idx]
      )) %>%
      arrange(row_id, column_id)

    Y[selected_idx] <- 0
    new_zero_prop <- mean(Y == 0, na.rm = TRUE)
    message(sprintf(
      "Now %.2f%% of Y are zeros.",
      new_zero_prop * 100
    ))
    summ_Y(Y)
    #====================================================================
    message("We rerun the previous block to remove any new columns with",
            " less than ", min_obs_per_col, " observations.")
    #=====================================================================
    # I will remove columns with less than min_obs_per_col observations.
    query.mat <- query_to_matrix(query, c(dim(Y)[1],dim(Y)[2]))
    stopifnot(all(matrix_to_query(query.mat) == query))
    message("Dim of Y (before): [" , dim(Y)[1], ",", dim(Y)[2],"]")

    columns_to_keep = base::which(colSums(Y!=0) >= min_obs_per_col)
    Y <- Y[, columns_to_keep]
    message("Dim of Y (after): [" , dim(Y)[1], ",", dim(Y)[2],"]")
    query.mat <- query.mat[, columns_to_keep]
    query <- matrix_to_query(query.mat)
    summ_Y(Y)
    keyword = paste0("_c_",min_obs_per_col,"_", round(100*prop_miss),"_")
    message("Saving data with keyword: ", keyword)
    saveRDS(Y, paste0("data/ml-1m/Movie_Y", keyword, ".Rdata"))
    saveRDS(query, paste0("data/ml-1m/Movie_Q",keyword,".Rdata"))
    #=========
  }
  message("Finally, we save the data as .dat for Python fit")
  obs_ind <- which(Y!=0, arr.ind=TRUE)
  py.Y <- data.frame(userID=obs_ind[,1], movieID=obs_ind[,2], rating=Y[obs_ind])
  colnames(query) <- c("userID", "movieID", "rating")
  write.table(py.Y,
              paste0("data/ml-1m/Movie_Y",keyword, ".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
  write.table(query,
              paste0("data/ml-1m/Movie_test",keyword,".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}
