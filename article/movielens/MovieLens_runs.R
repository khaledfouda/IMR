require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")
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
  data$Y <- IMR::as.Incomplete(data$Y)
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
data <- load_movielens1m()
seed = 2025

model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0.2)
print(model_data)

model_combn <- data.frame(
  row_covariates = c(F,T,T),
  col_covariates = c(F,F,T),
  intercepts = c(T,T,T),
)
#-----

convergence <- imr_convergence(maxit=600, thresh=1e-5, trace=FALSE, ls_initial = TRUE)


grid <- imr_tune_grid(beta = c(0, NA, 20),
                      gamma = c(0, NA, 20),
                      laplace = c(0, NA, 20, 3),
                      rank = c(10, 30, 1, 2))

model_data <- update(model_data,
                     intercept_row = TRUE,
                     intercept_col = TRUE,
                     col_intercepts = TRUE,
                     row_covariates = TRUE,
                     col_covariates = TRUE); model_data

# notes: max for model IMR-I: laplace = 120 but make it 45-16
#       max for model IMR-IX: beta = 2,  laplace = 135
#       max for model IMR-IXZ: beta: 1.9, laplace = 109, gamma = 4
grid$beta$max = 2
grid <- imr_set_grid_limits(model_data, grid,default_rank = 10,
                            bisection_iter = 5,
                            convergence=convergence, verbose=2)
grid$laplace$max = 45

start <- Sys.time()
bench::bench_time(fitimr2 <- IMR::imr_tune(model_data, grid, final_fit = TRUE,
                                          fast_laplace = FALSE,
                                          laplace_log_scale = TRUE,
                                          convergence=convergence, n_cores=7,
                                          seed = seed, verbose=2)) -> time.imr
time <- Sys.time() - start

print(fitimr$fit)
summary(fitimr$fit)

convergence2 <- imr_convergence(maxit=2000, thresh=1e-7, trace=FALSE)
start <- Sys.time()
bench::bench_time(fitimr <- IMR::imr_fit(model_data, rank = fitimr$params$rank,
                                         lambda_m = fitimr$params$lambda_laplace,
                                          convergence=convergence2)) -> time.imr
time <- Sys.time() - start


print(fitimr)
summary(fitimr)




datrec <- reconstruct(fitimr, model_data)
prepare_output_movielens(
  "IMR-I",
  time = time.imr[2],
  X = data$X,
  Z = data$Z,
  beta.estim = datrec$beta,
  gamma.estim = datrec$gamma,
  estim.test = datrec$estimates[data$test.idx],
  estim.train = as.Incomplete(datrec$estimates * data$obs_mask)@x,
  obs.test = data$test.truths,
  test_error = IMR:::error_metrics$rmse,
  obs.train = data$Y[data$Y != 0],
  M.estim = datrec$M,
  rank.M = 12 # fitimr$rank_M
)
