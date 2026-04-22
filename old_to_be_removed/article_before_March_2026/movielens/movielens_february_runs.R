require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")
#===== part 1: loading and preparing the data ===============
#===== out:   X, Z, Y, query, test.idx, test.truths, obs_mask
#=============================================================
load_movielens1m <- function(){
data <- list()
load("article/movielens/data/Movie_X.Rdata") #X
load("article/movielens/data/Movie_Y.Rdata",verbose = T)
X <- X[,1:4] # keep only main-effects
data$X <- X
data$Y <- Y
input_tag = "_c_0_"
#Y <-     readRDS(paste0("article/movielens/data/Movie_Y",input_tag,".Rdata"))
query <- readRDS(paste0("article/movielens/data/Movie_Q",input_tag,".Rdata"))
data$query <- query
source("other_models/Ma25.R")
source("other_models/SoftImpute_cv.R")
source("article/movielens/preprocess.R")
source("article/movielens/Ma25_fit.R")
#========================================================
# prepare test set and X-QR
data$test.idx   <- cbind(data$query[, 1], data$query[, 2])
data$test.truths <- data$query[, 3]
data$Y <- IMR::as_incomplete(data$Y)
data$obs_mask <- as.matrix((data$Y!=0) * 1)
mean(data$obs_mask==1)
mean(data$obs_mask==0)
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
  transmute(genre = as.vector(strsplit(data$Z$genres, "|", fixed=TRUE))) ->
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

data$Z <- cbind(data$Z[,1:2], genre_df)
data$Z %<>% mutate(movie_id = as.numeric(movie_id))
movies_no_genre <- (1:3952)[!(1:3952 %in% data$Z$movie_id)]
extra.rows <- data.frame(movie_id = movies_no_genre, title = "")
for(genre in genre_labels)
  extra.rows[genre] <- 0
rbind(data$Z, extra.rows) %>%
  arrange(movie_id) %>%
  dplyr::select(-movie_id, -title) %>%
  as.matrix() ->
  data$Z
return(data)
}
#=========
data <- load_movielens1m()
fit.imr1 <- readRDS("article/movielens/data/saved_models/IMR_fit_intercept.rds")
f1 <- fit.imr1$fit


start <- Sys.time()
{
data <- load_movielens1m()
data$mdat <- IMR::prepare_data(data$Y, data$X, data$Z, val_prop = 0.2, seed = 2025)
bench::bench_time(fitimr <- IMR::imr.fit(
  Y = data$Y,
  X = data$mdat$Xq,
  Z = data$mdat$Zq,
  r = 14,
  lambda_m = 17.88, #20.79,#14.53,#
  lambda_beta = .21053,#.03158,
  lambda_gamma = .3157892,
  row_intercept = T,
  col_intercept = T,
  trace = F,
  thresh = 1e-6,
  maxit = 3000,
  ls_initial = F
)) -> time.imr;print(time.imr)
f2 <- fitimr

saveRDS(fitimr, "article/movielens/data/february_fitimr1.rds")

data <- load_movielens1m()
data$mdat <- IMR::prepare_data(data$Y, data$X, data$Z, val_prop = 0.2, seed = 2025)
dat2 <- data$mdat;
dat2$X <- data$X; dat2$Z=data$Z
datrec <- IMR::reconstruct(fitimr, dat2)

prepare_output_movielens(
  "IMR-I",
  time = 3,#time.imr[2],
  X = dat2$X,
  Z = dat2$Z,
  beta.estim  = datrec$beta,
  gamma.estim = datrec$gamma,
  estim.test = datrec$estimates[data$test.idx],
  estim.train = as_incomplete(datrec$estimates * data$obs_mask)@x,
  obs.test = data$test.truths,
  test_error = IMR:::error_metric$rmse,
  obs.train = data$Y[data$Y!=0],
  M.estim = datrec$M,
  rank.M = 12#fitimr$rank_M
)
}

IMR:::error_metric$rmse(datrec$estimates[data$test.idx], data$test.truths)
IMR:::error_metric$rmse(as_incomplete(datrec$estimates * data$obs_mask)@x,data$Y[data$Y!=0])
#========================================

fit.imr1 <- readRDS("article/movielens/data/february_fitimr1.rds")
fit.imr2 <- readRDS("article/movielens/data/february_fitimr2.rds")
fit.imr3 <- readRDS("article/movielens/data/february_fitimr3.rds")
fit.imr1$time = 3.63
fit.imr1$lambda_laplace = 17.8827
fit.imr1$rank_M = length(fit.imr1$d > 0)
fit.imr1$lambda_beta = .21053
fit.imr1$lambda_gamma = .3157895
#
fit.imr2$time = 2.06
fit.imr2$lambda_laplace = 14.53
fit.imr2$rank_M = length(fit.imr2$d > 0)
fit.imr2$lambda_beta = .03158
fit.imr2$lambda_gamma = NA
#
fit.imr3$time = 1.52
fit.imr3$lambda_laplace = 20.79
fit.imr3$rank_M = length(fit.imr3$d > 0)
fit.imr3$lambda_beta = NA
fit.imr3$lambda_gamma = NA
#
fit.si   <- readRDS("article/movielens/data/saved_models/SI_fit.rds")
fit.ma   <- readRDS("article/movielens/data/saved_models/Ma_fit.rds")

data <- load_movielens1m()
data$mdat <- IMR::prepare_data(data$Y, data$X, data$Z, val_prop = 0.2, seed = 2025)
data$mdat$X <- data$X; data$mdat=data$Z

fit.si$time_old <- fit.si$time
as.numeric(bench::bench_time(fit.si$fit2 <-
                    softImpute::softImpute(as.matrix(data$Y), 9, 0, thresh=1e-6, maxit=3000))[2])->
  fit.si$time

out <- list()
out[[1]] <- IMR::reconstruct(fit.imr1, data$mdat)
out[[2]] <- IMR::reconstruct(fit.imr2, data$mdat)
out[[3]] <- IMR::reconstruct(fit.imr3, data$mdat)
out[[4]] <- IMR::reconstruct(fit.si$fit, data$mdat)

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
models <- c("IMR-IXZ", "IMR-IX", "IMR-I", "SoftImpute", "Ma")

for(i in 1:5){
  print(models[[i]])
  res[[i]] <- prepare_output_movielens(
    models[[i]],
    time = fits[[i]]$time,
    X = data$X,
    Z = data$Z,
    beta.estim  = out[[i]]$beta,
    gamma.estim = out[[i]]$gamma,
    estim.test = out[[i]]$estimates[data$test.idx],
    estim.train = as_incomplete(out[[i]]$estimates * data$obs_mask)@x,
    obs.test = data$test.truths,
    test_error = IMR:::error_metric$rmse,
    obs.train = data$Y[data$Y!=0],
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

out <- list(dat=data, fits=fits, res=res.df, res.list=res, out=out)





