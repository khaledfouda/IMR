require(tidyverse)
devtools::load_all()
require(magrittr)
load("article_results/movielens/data/Movie_X.Rdata") #X
load("article_results/movielens/data/Movie_Y.Rdata",verbose = T)
X <- X[,1:4] # keep only main-effects
input_tag = "_c_0_"
#Y <-     readRDS(paste0("article_results/movielens/data/Movie_Y",input_tag,".Rdata"))
query <- readRDS(paste0("article_results/movielens/data/Movie_Q",input_tag,".Rdata"))
source("other_models/Ma25.R")
source("other_models/SoftImpute_cv.R")
#========================================================
# prepare test set and X-QR
idx   <- cbind(query[, 1], query[, 2])
truths <- query[, 3]
Y <- IMR::as.Incomplete(Y)
obs_mask <- as.matrix((Y!=0) * 1)
mean(obs_mask==1)
mean(obs_mask==0)
# we need X to be orthonormal
Xqr <- qr(as.matrix(X))

dat <- list(
  Y = Y,
  Y.mat = as.matrix(Y),
  obs_mask = as.matrix(obs_mask),
  Xq = qr.Q(Xqr),
  Xr = qr.R(Xqr),
  X  = X
)
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
Zqr <- qr(as.matrix(Z))
dat$Z <- Z
dat$Zq = qr.Q(Zqr)
dat$Zr = qr.R(Zqr)
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
bench::bench_time(fit.imr1 <- IMR::imr.cv(
  Y = Y,
  #X = dat$Xq,
  #Z = dat$Zq,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  seed = 2025
)) -> time.imr

fit.imr1$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr1, paste0("article_results/movielens/data/saved_models/",
                        "IMR_fit_intercept.rds"))

# 2. row covariates >
bench::bench_time(fit.imr2 <- IMR::imr.cv(
  Y = Y,
  X = dat$Xq,
  #Z = dat$Zq,
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
bench::bench_time(fit.imr3 <- IMR::imr.cv(
  Y = Y,
  X = dat$Xq,
  Z = dat$Zq,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar,
  verbose = 1,
  seed = 2025
)) -> time.imr

fit.imr3$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr3, paste0("article_results/movielens/data/saved_models/",
                        "IMR_fit_rows_cols.rds"))

#--- soft-Impute
bench::bench_time(fit.si <-
                    simpute.cv(dat$Y,
                               seed = 2025,
                               n.lambda = hpar$M$n.lambda,
                               trace = T)) -> time.si
fit.si$time <- round(lubridate::time_length(time.si[2],  "minute"),2)
saveRDS(fit.si, "article_results/movielens/data/saved_models/SI_fit.rds")
#------------ Ma
M = fit_MA25_movielens("", 2025)
#--------------------------------------------------------------


out <- IMR:::reconstruct(fit.imr$fit, dat)
prepare_output_movielens(
  "IMR",
  time       = fit.imr$time,
  X           = dat$X,
  Z           = dat$Z,
  estim.test  = out$estimates[idx],
  estim.train = out$estimates[dat$obs_mask==1],
  obs.test    = truths,
  obs.train   = dat$Y[dat$obs_mask==1],
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


# fit MA
bench::bench_time(softImpute::softImpute(dat$Y,
                                         rank.max = lambdas$r[s],
                                         lambda = lambdas$M[s])) -> sit


M <- readRDS(paste0("article_results/movielens/data/saved_models/Ma_fit.rds"))

rhat <- M$rank_est$est['h'] # get the estimated rank
M$fit[[1]]$rmse # training error
fit.ma <- M$fit[[1]]
fit.ma$M <- fit.ma$L %*% t(fit.ma$R)
fit.ma$xbeta <- cbind(1, X) %*% t(fit.ma$B)
fit.ma$estimates <- fit.ma$M + fit.ma$xbeta
fit.ma$beta <- fit.ma$B[,-1]

prepare_output_movielens(
  "Ma",
  time       = M$,
  X           = X,
  estim.test  = fit.ma$estimates[idx],
  estim.train = fit.ma$estimates[dat$obs_mask==1],
  obs.test    = truths,
  obs.train   = dat$Y[dat$obs_mask==1],
  beta.estim  = fit.ma$B,
  M.estim     = fit.ma$M,
  rank.M      = rhat
) -> results.ma
results.ma; results.ma








