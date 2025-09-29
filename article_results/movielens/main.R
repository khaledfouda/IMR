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
#============================================================
# part 2: fitting the models then saving them to disk =======
#============================================================
# fit IMR
hpar <- IMR::get_imr_default_hparams()
hpar$beta$lambda_max <- 1
hpar$gamma$lambda_max <- 1
hpar$M$n.lambda <- 60
hpar$beta$n.lambda <- 20
hpar$gamma$n.lambda <- 20

IMR:::initialize_parallel_workers(9)


# fit all model variations: [imr+x, imr+xz, imr+intercept, ma, si]
#----------------------
# 1. intercept only >
dat <- IMR::prepare_data(Y, NULL, NULL, 0.2, seed = 2025)



bench::bench_time(fit.imr1 <- IMR:::imr.cv(
  dat$model,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  fast.cv = T,
  separate_tuning = TRUE,
  seed = 2025
)) -> time.imr

fit.imr1$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr1, paste0("article_results/movielens/data/saved_models/",
                        "IMR_fit_intercept.rds"))

# 2. row covariates >
dat <- IMR::prepare_data(Y, X, NULL, 0.2, seed = 2025)
bench::bench_time(fit.imr2 <- IMR:::nlrr.cv(
  dat$model,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  seed = 2025
)) -> time.imr

fit.imr2$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr2, paste0("article_results/movielens/data/saved_models/",
                        "IMR_fit_rows.rds"))


# 3. row+col covariates >
dat <- IMR::prepare_data(Y, X, Z, 0.2, seed = 2025)
bench::bench_time(fit.imr3 <- IMR::imr.cv(
  dat$model,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  seed = 2025,
  separate_tuning = FALSE,
)) -> time.imr

fit.imr3$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr3, paste0("article_results/movielens/data/saved_models/",
                        "IMR_fit_rows_cols.rds"))

#--- soft-Impute
bench::bench_time(fit.si <-
                    simpute.cv(dat$Y,
                               seed = 2025,
                               lambda_max = 100,
                               n.lambda = hpar$M$n.lambda,
                               trace = T)) -> time.si
fit.si$time <- round(lubridate::time_length(time.si[2],  "minute"),2)
saveRDS(fit.si, "article_results/movielens/data/saved_models/SI_fit.rds")
#------------ Ma
M = fit_MA25_movielens("", 2025)
#--------------------------------------------------------------
#==============================================================
#================ part 3 - load the models and show the results
#==============================================================
dat <- IMR::prepare_data(Y, X, Z, 0.2, seed = 2025)
fit.imr1 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_intercept.rds")
fit.imr2 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_rows.rds")
fit.imr3 <- readRDS("article_results/movielens/data/saved_models/IMR_fit_rows_cols.rds")
fit.si   <- readRDS("article_results/movielens/data/saved_models/SI_fit.rds")
fit.ma   <- readRDS("article_results/movielens/data/saved_models/Ma_fit.rds")



out <- list()
out[[1]] <- IMR::reconstruct(fit.imr1$fit, dat)
out[[2]] <- IMR::reconstruct(fit.imr2$fit, dat)
out[[3]] <- IMR::reconstruct(fit.imr3$fit, dat)
out[[4]]   <- IMR::reconstruct(fit.si$fit, dat)
out[[5]] <- IMR::reconstruct(fit.ma$fit, dat)
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


best_min_cols <- c("time", "error.test", "error.train", "rank_M" )
best_max_cols <- c("corr.test")

best_idx_min <- lapply(best_min_cols, function(v)
  which.min(replace(res.df[[v]], is.na(res.df[[v]]), Inf))
) |> rlang::set_names(best_min_cols)

best_idx_max <- lapply(best_max_cols, function(v)
  which.max(replace(res.df[[v]], is.na(res.df[[v]]) | is.nan(res.df[[v]]), -Inf))
) |> rlang::set_names(best_max_cols)
fmt_num <- function(x, digits) ifelse(is.na(x) | is.nan(x), "—", sprintf(paste0("%.", digits, "f"), x))

disp <- res.df %>%
  mutate(
    time          = fmt_num(time, 2),
    error.test    = fmt_num(error.test, 3),
    corr.test     = fmt_num(corr.test, 3),
    error.train   = fmt_num(error.train, 3),
    sparsity_beta = fmt_num(sparsity_beta, 3),
    sparsity_gamma= fmt_num(sparsity_gamma, 3),
    rank_M        = ifelse(is.na(rank_M), "—", as.character(rank_M)),
    rank_beta     = ifelse(is.na(rank_beta), "—", as.character(rank_beta)),
    rank_gamma    = ifelse(is.na(rank_gamma), "—", as.character(rank_gamma)),
    rank_total    = as.character(rank_total)
  ) %>% dplyr::select( - rank_total)
for (v in names(best_idx_min)) {
  idx <- best_idx_min[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}
for (v in names(best_idx_max)) {
  idx <- best_idx_max[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}

require(kableExtra)
col_names <- c(
  "Model", "Time (min)",
  "RMSE", "correlation", "RMSE",
  "$M$", "$\\beta$", "$\\Gamma$",
  "$\\beta$", "$\\Gamma$"
)

kbl(disp,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    escape   = FALSE,
    col.names = col_names,
    caption  = paste0("Performance comparison on the MovieLens 1M dataset. ",
    "Best values per column are bolded and IMR models are shaded.")
) |>
  add_header_above(c(" " = 2, "Test"=2, "Train"=1, " "=5)) |>
  add_header_above(c(" " = 2, "Performance" = 3, "Rank estimation" = 3, "Sparsity" = 2)) |>
  kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8) |>
  row_spec(which(grepl("^IMR", disp$model)), bold = FALSE, background = "#f7f7f7") |>
  column_spec(1, width = "3.2cm")# |>
  # footnote(
  #   general = "Sparsity is the fraction of zeros (higher is sparser). Dashes denote unavailable metrics.",
  #   threeparttable = TRUE
  # )
############################################################################
#================== part 4 > insights ? =====================================
#############################################################################
fit <- out[[2]]
G <- fit$beta[1,]
which.max(G)
which.min(G)

#==============================================================================







