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

# 1- results table:
res.df <- data.frame()
do.call(rbind, lapply(res, function(x) x[1:12])) %>%
  as.data.frame() %>%
  mutate(across(c(time, error.test, corr.test, error.train, sparsity_beta,
                  sparsity_gamma), as.numeric)) %>%
  mutate(across(where(is.numeric), ~ round(.x,3))) %>%
  dplyr::select(-time_per_fit, -total_num_fits) %>%
  arrange(error.test)



#==============================================================================
IMR::initialize_parallel_workers(9)
mod.dat <- IMR::prepare_data(Y, X, NULL, 0.2, seed = 2025)

devtools::load_all()
bench::bench_time(fit.imr1 <- IMR:::imr.cv(
  mod.dat$model,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  fast.cv = F,
  separate_tuning = TRUE,
  seed = 2025
))
r


dat <- mod.dat

fit <- fit.imr
fit$u <- fit.imr$fit$u
fit$v <- fit.imr$fit$v
fit$d <- fit.imr$fit$d

outsi <- IMR::reconstruct(fit.si3$fit, dat)
prepare_output_movielens(
  "SoftImpute",
  time = fit.si$time,
  X = dat$X,
  estim.test = outsi$estimates[test.idx],
  estim.train = as.Incomplete(out$estimates * obs_mask)@x,
  obs.test = test.truths,
  obs.train = dat$Y[dat$Y!=0],
  M.estim = out$M,
  rank.M = fit.si$rank_M
)

fit.imr <- fit.imr1; out <- IMR:::reconstruct(fit.imr$fit, dat)
prepare_output_movielens(
  "IMR",
  time       = fit.imr$time,
  X           = dat$X,
  Z           = dat$Z,
  estim.test  = out$estimates[test.idx],
  estim.train = out$estimates[as.matrix(dat$Y!=0)],
  obs.test    = test.truths,
  obs.train   = dat$Y[dat$Y!=0],
  beta.estim  = out$beta,
  gamma.estim = out$gamma,
  M.estim     = out$M,
  rank.M      = fit.imr$rank_M,
  total_num_fits = fit.imr$total_num_fits,
  time_per_fit   = fit.imr$time_per_fit
) -> results.imr; results.imr

results.imr$cov_summaries_cols %>%
  as.data.frame() %>%
  arrange(desc(Mean))


IMR::as.Incomplete(out$estimates * obs_mask)@x[1:15]
IMR::as.Incomplete(out1$estimates * obs_mask)@x[1:15]


out$M[1:5,1:5]
fit.imr1$fit$beta[,1:10]

fit.imr11$rank_M

# fit MA
bench::bench_time(softImpute::softImpute(dat$Y,
                                         rank.max = lambdas$r[s],
                                         lambda = lambdas$M[s])) -> sit


M <- readRDS(paste0("article_results/movielens/data/saved_models/Ma_fit.rds"))

rhat <- M$rank_est$est['h'] # get the estimated rank


M <- M$fit
M$fit[[1]]$rmse # training error
fit.ma <- M$fit[[1]]
fit.ma$M <- fit.ma$L %*% t(fit.ma$R)
fit.ma$xbeta <- cbind(1, X) %*% t(fit.ma$B)
fit.ma$estimates <- fit.ma$M + fit.ma$xbeta
fit.ma$beta <- fit.ma$B[,-1]
fit.ma$estimates[1:5,1:5]
fit.ma$beta[1:5,]
out[[5]]$beta[,1:5] %>% t()
out[[5]]$estimates[1:5,1:5]

prepare_output_movielens(
  "Ma",
  time       = NULL,
  X           = X,
  estim.test  = fit.ma$estimates[test.idx],
  estim.train = fit.ma$estimates[dat$obs_mask==1],
  obs.test    = test.truths,
  obs.train   = dat$Y[dat$obs_mask==1],
  beta.estim  = fit.ma$B,
  M.estim     = fit.ma$M,
  rank.M      = rhat
) -> results.ma
results.ma; results.ma








