
prepare_ml_1m_data <- function(min_obs_per_col = 10,
                               increase_missing = FALSE,
                               prop_miss = .98,
                               seed = 2025){

  if(is.numeric(seed)) set.seed(seed)
  require(Matrix)
  require(dplyr)
  require(magrittr)

  load("article_results/movielens/data/Movie_Y.Rdata") # Y
  load("article_results/movielens/data/Movie_X.Rdata") # X
  load("article_results/movielens/data/Movie_Q.Rdata") # Query (testing set for evaluating the performance)
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
  saveRDS(Y, paste0("article_results/movielens/data/Movie_Y",keyword,".Rdata"))
  saveRDS(query, paste0("article_results/movielens/data/Movie_Q",keyword,".Rdata"))
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
    saveRDS(Y, paste0("article_results/movielens/data/Movie_Y", keyword, ".Rdata"))
    saveRDS(query, paste0("article_results/movielens/data/Movie_Q",keyword,".Rdata"))
    #=========
  }
  message("Finally, we save the data as .dat for Python fit")
  obs_ind <- which(Y!=0, arr.ind=TRUE)
  py.Y <- data.frame(userID=obs_ind[,1], movieID=obs_ind[,2], rating=Y[obs_ind])
  colnames(query) <- c("userID", "movieID", "rating")
  write.table(py.Y,
              paste0("article_results/movielens/data/Movie_Y",keyword, ".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
  write.table(query,
              paste0("article_results/movielens/data/Movie_test",keyword,".dat"),
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}
#-------------------------------------------------------------------------------



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
prepare_results_for_analysis <- function(){
  require(tidyverse)
  devtools::load_all()
  require(magrittr)

  #===== part 1: loading and preparing the data ===============
  #===== out:   X, Z, Y, query, test.idx, test.truths, obs_mask
  #=============================================================

  load("article_results/movielens/data/Movie_X.Rdata") #X
  load("article_results/movielens/data/Movie_Y.Rdata",verbose = T)
  X <- X[,1:4] # keep only main-effects
  input_tag = "_c_0_"
  #Y <-     readRDS(paste0("article_results/movielens/data/Movie_Y",input_tag,".Rdata"))
  query <- readRDS(paste0("article_results/movielens/data/Movie_Q",input_tag,".Rdata"))
  source("other_models/Ma25.R")
  source("other_models/SoftImpute_cv.R")
  source("article_results/movielens/preprocess.R")
  source("article_results/movielens/Ma25_fit.R")
  #========================================================
  # prepare test set and X-QR
  test.idx   <- cbind(query[, 1], query[, 2])
  test.truths <- query[, 3]
  Y <- IMR::as.Incomplete(Y)
  obs_mask <- as.matrix((Y!=0) * 1)
  mean(obs_mask==1)
  mean(obs_mask==0)
  # ====================================================
  # prepare Z < genres >
  Z <- data.table::fread(
    file = "article_results/movielens/data/movies_Z.dat",
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
  Z %>%
    transmute(genre = as.vector(strsplit(Z$genres, "|", fixed=TRUE))) ->
    genres
  genres = genres[[1L]]
  i <- rep(seq_along(genres), lengths(genres))
  j <- match(unlist(genres, use.names = FALSE), genre_labels)
  keep <- !is.na(j)
  genre_sparse <- Matrix::sparseMatrix(
    i = i[keep], j = j[keep], x = 1L,
    dims = c(length(genres), length(genre_labels)),
    dimnames = list(NULL, genre_labels)
  )
  genre_df <- as.data.frame(as.matrix(genre_sparse), check.names = FALSE)

  Z <- cbind(Z[,1:2], genre_df)
  Z %<>% mutate(movie_id = as.numeric(movie_id))
  movies_no_genre <- (1:3952)[!(1:3952 %in% Z$movie_id)]
  extra.rows <- data.frame(movie_id = movies_no_genre, title = "")
  for(genre in genre_labels)
    extra.rows[genre] <- 0
  rbind(Z, extra.rows) %>%
    arrange(movie_id) %>%
    dplyr::select(-movie_id, -title) %>%
    as.matrix() ->
    Z
  #=========
  dat <- IMR::prepare_data(Y, X, Z, 0.2, seed = 2025)
  fit.imr1 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_intercept.rds")
  fit.imr2 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_rows_nonsp.rds")
  fit.imr3 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_rows_cols.rds")
  fit.si   <- readRDS("article_results/movielens/data/saved_models/SI_fit.rds")
  fit.ma   <- readRDS("article_results/movielens/data/saved_models/Ma_fit.rds")



  out <- list()
  out[[1]] <- IMR::reconstruct(fit.imr1$fit, dat)
  out[[2]] <- IMR::reconstruct(fit.imr2$fit, dat)
  out[[3]] <- IMR::reconstruct(fit.imr3$fit, dat)
  out[[4]]   <- IMR::reconstruct(fit.si$fit, dat)
  # out[[5]] <- IMR::reconstruct(fit.ma$fit, dat)
  out[[5]] <- list()


  fit.ma$rank_M = fit.ma$fit$fit[[1]]$rank
  ffit.ma <- fit.ma$fit$fit[[1]]
  ffit.ma$M <- ffit.ma$L %*% t(ffit.ma$R)
  ffit.ma$xbeta <- cbind(1, X) %*% t(ffit.ma$B)
  ffit.ma$estimates <- ffit.ma$M + ffit.ma$xbeta
  ffit.ma$beta <- t(ffit.ma$B[,-1])
  fit.ma$fit <- out[[5]] <- ffit.ma

  res <- list()
  fits <- list(fit.imr1, fit.imr2, fit.imr3, fit.si, fit.ma)
  models <- c("IMR-Intercept", "IMR-X", "IMR-XZ", "SoftImpute", "Ma")

  for(i in 1:5){
    print(models[[i]])
    res[[i]] <- prepare_output_movielens(
      models[[i]],
      time = fits[[i]]$time,
      X = dat$X,
      Z = dat$Z,
      beta.estim  = out[[i]]$beta,
      gamma.estim = out[[i]]$gamma,
      estim.test = out[[i]]$estimates[test.idx],
      estim.train = as.Incomplete(out[[i]]$estimates * obs_mask)@x,
      obs.test = test.truths,
      obs.train = dat$Y[dat$Y!=0],
      M.estim = out[[i]]$M,
      rank.M = fits[[i]]$rank_M
    )
  }
  res[[5]]$rank_beta = 5
  res[[6]] <- list(model = "Glocal-K",
                   time        = 52.36,
                   error.test  = 0.8516,
                   corr.test   = 0.6278,
                   error.train = 0.7018,
                   rank_M      = 63,
                   rank_beta =NA,
                   rank_gamma = NA,
                   sparsity_beta = NA,
                   sparsity_gamma =NA)
  # 1- results table:
  do.call(rbind, lapply(res, function(x) if(length(x)==14) x[c(1:2,5:12)] else x)) %>%
    as.data.frame() %>%
    mutate(across(c(time, error.test, corr.test, error.train, sparsity_beta,
                    sparsity_gamma, rank_M, rank_beta, rank_gamma), as.numeric)) %>%
    mutate(across(where(is.numeric), ~ round(.x,3))) %>%
    mutate(rank_total = rank_M + ifelse(is.na(rank_beta), 0, rank_beta) +
             ifelse(is.na(rank_gamma), 0, rank_gamma)) %>%
    # dplyr::select(-time_per_fit, -total_num_fits) %>%
    arrange(error.test) ->
    res.df

  return(list(dat=dat, fits=fits, res=res.df, res.list=res, out=out))
}

