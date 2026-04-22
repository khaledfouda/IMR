
prepare_ml_1m_data <- function(min_obs_per_col = 10,
                               increase_missing = FALSE,
                               prop_miss = .98,
                               seed = 2025){

  if(is.numeric(seed)) set.seed(seed)
  require(Matrix)
  require(dplyr)
  require(magrittr)

  load("article/movielens/data/Movie_Y.Rdata") # Y
  load("article/movielens/data/Movie_X.Rdata") # X
  load("article/movielens/data/Movie_Q.Rdata") # Query (testing set for evaluating the performance)
  query <- as.data.frame(query)
  colnames(query) <- c("row_id", "column_id", "value")
  # remove interactions from data
  X <- X[,1:4]
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
  saveRDS(Y, paste0("article/movielens/data/Movie_Y",keyword,".Rdata"))
  saveRDS(query, paste0("article/movielens/data/Movie_Q",keyword,".Rdata"))
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
    saveRDS(Y, paste0("article/movielens/data/Movie_Y", keyword, ".Rdata"))
    saveRDS(query, paste0("article/movielens/data/Movie_Q",keyword,".Rdata"))
    #=========
  }
  message("Finally, we save the data as .dat for Python fit")
  obs_ind <- which(Y!=0, arr.ind=TRUE)
  py.Y <- data.frame(userID=obs_ind[,1], movieID=obs_ind[,2], rating=Y[obs_ind])
  colnames(query) <- c("userID", "movieID", "rating")
  write.table(py.Y,
              paste0("article/movielens/data/Movie_Y",keyword, ".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
  write.table(query,
              paste0("article/movielens/data/Movie_test",keyword,".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}
#-------------------------------------------------------------------------------

# ===== part 1: loading and preparing the data ===============
# ===== out:   X, Z, Y, query, test.idx, test.truths, obs_mask
# =============================================================
load_movielens1m <- function() {
  data <- list()
  load("article/movielens/data/Movie_X.Rdata") # X
  load("article/movielens/data/Movie_Y.Rdata", verbose = T)
  X <- X[, 1:4] # keep only main-effects
  data$X <- X
  data$Y <- Y
  input_tag <- "_c_0_"
  # Y <-     readRDS(paste0("article/movielens/data/Movie_Y",input_tag,".Rdata"))
  query <- readRDS(paste0("article/movielens/data/Movie_Q", input_tag, ".Rdata"))
  data$query <- query
  source("other_models/Ma25.R")
  source("other_models/SoftImpute_cv.R")
  source("article/movielens/preprocess.R")
  source("article/movielens/Ma25_fit.R")
  # ========================================================
  # prepare test set and X-QR
  data$test.idx <- cbind(data$query[, 1], data$query[, 2])
  data$test.truths <- data$query[, 3]
  data$Y <- IMR::as_incomplete(data$Y)
  data$obs_mask <- as.matrix((data$Y != 0) * 1)
  mean(data$obs_mask == 1)
  mean(data$obs_mask == 0)
  # ====================================================
  # prepare Z < genres >
  data$Z <- data.table::fread(
    file = "article/movielens/data/movies_Z.dat",
    sep = NULL,
    encoding = "Latin-1",
    header = FALSE
  ) %>%
    tidyr::separate(
      V1,
      into = c("movie_id", "title", "genres"),
      sep = "::"
    )

  genre_labels <- c(
    "Action", "Adventure", "Animation", "Children's", "Comedy", "Crime",
    "Documentary", "Drama", "Fantasy", "Film-Noir", "Horror", "Musical",
    "Mystery", "Romance", "Sci-Fi", "Thriller", "War", "Western"
  )
  data$Z %>%
    transmute(genre = as.vector(strsplit(data$Z$genres, "|", fixed = TRUE))) ->
    genres
  genres <- genres[[1L]]
  i <- rep(seq_along(genres), lengths(genres))
  j <- match(unlist(genres, use.names = FALSE), genre_labels)
  keep <- !is.na(j)
  genre_sparse <- Matrix::sparseMatrix(
    i = i[keep], j = j[keep], x = 1L,
    dims = c(length(genres), length(genre_labels)),
    dimnames = list(NULL, genre_labels)
  )
  genre_df <- as.data.frame(as.matrix(genre_sparse), check.names = FALSE)

  data$Z <- cbind(data$Z[, 1:2], genre_df)
  data$Z %<>% mutate(movie_id = as.numeric(movie_id))
  movies_no_genre <- (1:3952)[!(1:3952 %in% data$Z$movie_id)]
  extra.rows <- data.frame(movie_id = movies_no_genre, title = "")
  for (genre in genre_labels) {
    extra.rows[genre] <- 0
  }
  rbind(data$Z, extra.rows) %>%
    arrange(movie_id) %>%
    dplyr::select(-movie_id, -title) %>%
    as.matrix() ->
    data$Z
  return(data)
}
# =========

prepare_output_movielens <- function(
    model_name,
    X,
    Z = NA,
    estim.test,
    estim.train,
    obs.test,
    obs.train,
    time = NA,
    beta.estim  = NA,
    gamma.estim = NA,
    M.estim     = NA,
    rank.M      = NA,
    test_error  = IMR::error_metric$rmse,
    time_per_fit = NA,
    total_num_fits = NA
) {
  # Core metrics
  results <- list(
    model = model_name,
    time = time,
    time_per_fit = time_per_fit,
    total_num_fits = total_num_fits,
    error.test  = test_error(estim.test, obs.test),
    corr.test   = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    #rank_M      = tryCatch(
    #  qr(M.estim)$rank,
    #  error = function(e) NA
    #),
    rank_M = rank.M,
    rank_beta   = tryCatch(
      qr(beta.estim)$rank,
      error = function(e) NA
    ),
    rank_gamma   = tryCatch(
      qr(gamma.estim)$rank,
      error = function(e) NA
    ),
    sparsity_beta    = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    ),
    sparsity_gamma    = tryCatch(
      sum(gamma.estim == 0) / length(gamma.estim),
      error = function(e) NA
    )
  )


  # Covariate coefficient summaries
  results$cov_summaries_rows <- tryCatch({
    apply(beta.estim, 1, summary) |>
      as.data.frame() |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        prop_non_zero = apply(beta.estim, 1, function(x)
          sum(x != 0) / length(x)
        )
      ) |>
      `rownames<-`(colnames(X))
  }, error = function(e) NA)

  results$cov_summaries_cols <- tryCatch({
    apply(gamma.estim, 2, summary) |>
      as.data.frame() |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        prop_non_zero = apply(gamma.estim, 2, function(x)
          sum(x != 0) / length(x)
        )
      ) |>
      `rownames<-`(colnames(Z))
  }, error = function(e) NA)

  results
}
#---------
