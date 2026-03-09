
require(tidyverse)
require(magrittr)
fit_IMR_movielens <- function(input_tag = "_c_0_",
                                   seed = 2025,
                              intercept_row = FALSE,
                              intercept_col = FALSE,
                                   hpar = IMR::get_imr_default_hparams()){

  if(is.numeric(seed)) set.seed(seed)
  load("article_results/movielens/data/Movie_X.Rdata") #X
  load("article_results/movielens/data/Movie_Y.Rdata",verbose = T)
  Y0 <- Y
  X <- X[,1:4]
  input_tag = "_c_0_"
  Y <-     readRDS(paste0("article_results/movielens/data/Movie_Y",input_tag,".Rdata"))
  query <- readRDS(paste0("article_results/movielens/data/Movie_Q",input_tag,".Rdata"))
  #========================================================
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
    select(-movie_id, -title) %>%
    as.matrix() ->
    Z

  #========================================================
  # prepare data
  idx   <- cbind(query[, 1], query[, 2])
  truths <- query[, 3]
  Y <- IMR::as.Incomplete(Y)
  obs_mask <- as.matrix((Y!=0) * 1)
  mean(obs_mask==1)
  mean(obs_mask==0)
  # we need X to be orthonormal
  Xqr <- qr(as.matrix(X))
  Zqr <- qr(as.matrix(Z))

  dat <- list(
    Y = Y,
    Y.mat = as.matrix(Y),
    obs_mask = as.matrix(obs_mask),
    Xq = qr.Q(Xqr),
    Xr = qr.R(Xqr),
    Zq = qr.Q(Zqr),
    Zr = qr.R(Zqr),
    X  = X,
    Z = Z
  )
  #====================================================
  # fitting CAMC
  hpar$beta$lambda_max <- 2
  hpar$gamma$lambda_max <- 2
  hpar$M$n.lambda <- 40
  hpar$beta$n.lambda <- 10
  hpar$gamma$n.lambda <- 10

  future::plan(future::sequential)
  future::plan(future::multisession, workers = 9)

  fit.imr <- IMR::imr.cv(
    Y = Y,
    X = dat$Xq,
    Z = dat$Zq,
    intercept_row = T,
    intercept_col = T,
    hpar = hpar,
    verbose = 1,
    max_cores = 9,
    seed = 2025
  )
  out <- IMR:::reconstruct(fit.imr$fit, dat)
  prepare_output_movielens(
    "CAMC",
    time       = fit.imr$time,
    X           = X,
    estim.test  = out$estimates[idx],
    estim.train = out$estimates[dat$obs_mask==1],
    obs.test    = truths,
    obs.train   = dat$Y[dat$obs_mask==1],
    beta.estim  = out$beta,
    M.estim     = out$M,
    rank.M      = fit.imr$rank_M,
    total_num_fits = fit.imr$total_num_fits,
    time_per_fit   = fit.imr$time_per_fit
  ) -> results.camc; results.camc

colMeans(out$gamma) %>% sort()
  M <- readRDS(paste0("article_results/movielens/data/saved_models/Ma_fit",input_tag,".rds"))

  rhat <- M$rank_est$est['h'] # get the estimated rank
  M$fit[[1]]$rmse # training error
  fit.ma <- M$fit[[1]]
  fit.ma$M <- fit.ma$L %*% t(fit.ma$R)
  fit.ma$xbeta <- cbind(1, X) %*% t(fit.ma$B)
  fit.ma$estimates <- fit.ma$M + fit.ma$xbeta
  fit.ma$beta <- fit.ma$B[,-1]

  prepare_output_movielens(
    "Ma",
    time       = NA,
    X           = X,
    estim.test  = fit.ma$estimates[idx],
    estim.train = fit.ma$estimates[dat$obs_mask==1],
    obs.test    = truths,
    obs.train   = dat$Y[dat$obs_mask==1],
    beta.estim  = fit.ma$B,
    M.estim     = fit.ma$M,
    rank.M      = rhat
  ) -> results.ma
  results.ma


  fit.camc <- CAMC_cv_large(
    Y = dat$Y,
    X = dat$Xq,
    obs_mask = dat$obs_mask,
    n_folds = n_folds,
    hpar = hpar,
    verbose = 1,
    max_cores = 8,
    seed = 2025
  )
  saveRDS(fit.camc, paste0("data/saved_models_1m/CAMC_fit",input_tag,".rds"))


  #====== fit Mao

  fit.Mao <- Mao.cv(
    Y               = dat$Y.mat,
    X               = X,
    W               = dat$obs_mask,
    n_folds         = n_folds,
    lambda_1_grid   = seq(0,   1,   length = 20),
    lambda_2_grid   = seq(0.9, 0.1, length = 20),
    alpha_grid      = 1,
    seed            = seed,
    numCores        = 1,
    n1n2_optimized  = TRUE,
    test_error      = utils$error_metric$rmse,
    theta_estimator = Mao_weights$uniform,
    sequential      = FALSE
  )
  saveRDS(fit.Mao, paste0("data/saved_models_1m/Mao_fit",input_tag,".rds"))

}

#=======
prepare_output_movielens <- function(
    model_name,
    time,
    X,
    estim.test,
    estim.train,
    obs.test,
    obs.train,
    beta.estim  = NA,
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
    sparsity    = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    )
  )


  # Covariate coefficient summaries
  results$cov_summaries <- tryCatch({
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

  results
}


#=========



prop_miss = .96
min_obs_per_col = 0#10
increase_missing = FALSE
keyword = ifelse(increase_missing,
                 paste0("_c_",min_obs_per_col,"_", round(100*prop_miss),"_"),
                 paste0("_c_",min_obs_per_col,"_"))


prepare_ml_1m_data(min_obs_per_col =  min_obs_per_col,
                   increase_missing = increase_missing,
                   prop_miss = prop_miss,
                   seed = 2025)
M <- fit_MA25_movielens("", seed = 2025)
