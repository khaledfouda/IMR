library(devtools)
# devtools::uninstall(); devtools::install()
# clean_dll(); Rcpp::compileAttributes(); document()
devtools::load_all()
library(tidyverse)
library(magrittr)
source("./article_results/bixi/data_bixi.R")
source("./other_models/SoftImpute_cv.R")
source("./article_results/bixi/fit_bktr.R")
source("./article_results/bixi/helpers_bixi.R")
#------------------------------------------------

miss = 0.5
seed = 2026

dat <- prepare_bixi_data(miss, "25Sep", seed = seed,
                         val_prop = 0.2,
                         x_keep = c("x_max_temp_f",
                                    "x_total_precip_mm","x_holiday" ))
#--------------------------------------------------------
#-- change similarity matrices
require(BKTR)
bkdat <- BixiData$new()

file_prefix <- paste0(
  "./article_results/bixi/data/splits/",
  round(100*miss),
  "percent_",
  timestamp,
  "_"
)
train.df <- readRDS(paste0(file_prefix, "train.rds"))

train.df %<>% rename(location=column, time = row)
#train.df <- setkey(as.data.table(train.df), location, time)

bkdat$temporal_positions_df %<>%
  filter(time %in% train.df$time)
bkdat$spatial_positions_df %<>%
  filter(location %in% train.df$location)

spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 3)
temporal_kernel = BKTR::KernelSE$new()
spatial_kernel$set_positions(bkdat$spatial_positions_df)
temporal_kernel$set_positions(bkdat$temporal_positions_df)

distance.mat = spatial_kernel$distance_matrix %>% as.matrix()
min(distance.mat[upper.tri(distance.mat)])
d <- distance.mat
dvec <- d[upper.tri(d)]
dvec <- dvec[dvec > 0]
range0 <- median(dvec)
K <- fields::Matern(d, smoothness = 3/2, range = range0)
#tau2 <- 1e-3
#K <- K + diag(tau2, nrow(K))
dat$modd$similarity_cols <- chol2inv(K)
#dat$modd$similarity_rows <- diag(1, 193, 193)


distance.mat = temporal_kernel$distance_matrix %>% as.matrix()
min(distance.mat[upper.tri(distance.mat)])
d <- distance.mat
dvec <- d[upper.tri(d)]
dvec <- dvec[dvec > 0]
range0 <- median(dvec)
K <- fields::Matern(d, smoothness = 3/2, range = range0)
#tau2 <- 1e-3
#K <- K + diag(tau2, nrow(K))
dat$modd$similarity_rows <- chol2inv(K)



#--------------


hparam <- IMR::get_imr_default_hparams(dat$modd$similarity_rows,
                                       dat$modd$similarity_cols, 0, 0)
hparam$laplace$step_sizes <- c(1, 0.1, 0.01)
hparam$laplace$min <- 0
hparam$laplace$max <- 5
hparam$rank$n_streaks <- hparam$laplace$n_streaks <- 1
hparam$rank$max <- 15
hparam$beta$n_streaks <- hparam$gamma$n_streaks <- 2
hparam$beta$max <- 0.31
hparam$gamma$max <- 0.2
hparam$beta$step_sizes <- 1e-7
hparam$gamma$step_sizes <- 0.2 / 300


initialize_parallel_workers(9)

bench::bench_time(res0 <- IMR:::imr.cv_laplace(dat$modd,final_fit = T,
                                             trace=2, hpar=hparam, intercept_row = T,
                                             intercept_col = T, ls_initial = T,
                                             seed = seed, warm_start = NULL, maxit=600,
                                             shared_information = T, thresh=1e-5,
                                             num_cores = 0)) -> timee
round(lubridate::time_length(timee[2], "seconds"), 2)
s0 <- output_wrapper_bixi(res0$best_fit, dat,T)
s0$res$error.test
vestim <- IMR::reconstruct_partial(res0$best_fit, dat$modd, dat$modd$y_valid, T)

error_metric$rmse(dat$modd$y_valid@x, vestim@x)



rec <- IMR::reconstruct(res$best_fit, dat$modd,shared_information = T)

rec$beta
rec$gamma
res$best_fit$beta0


res.bktr <- BKTR_Bixi_Wrapper(dat, "25Sep", 2025, miss,return_fit = T)

res.bktr$results$error.test


hparam$laplacian_row <- hparam$laplacian_col <- NULL
dat$modd$similarity_rows <- dat$modd$similarity_cols <- NULL

bench::bench_time(res <- IMR:::imr.cv(dat$modd,
                                              trace=2, hpar=hparam, intercept_row = T,
                                              intercept_col = T, ls_initial = T,
                                              seed = seed, warm_start = NULL, maxit=900,
                                              thresh = 1e-6,
                                              shared_information = T,
                                              num_cores = 0)) -> timee
round(lubridate::time_length(timee[2], "minutes"), 2)
s <- output_wrapper_bixi(res$best_fit, dat,T)
s$res$error.test

res$best_fit$gamma
res$best_fit$beta
s$rec$beta
s$rec$gamma


dat$modd$similarity_cols[1:5,1:5]


distance.mat = temporal_kernel$distance_matrix

distance.mat = spatial_kernel$distance_matrix %>% as.matrix()

distance.mat = fields::rdist(as.matrix(1:579))
kernel =   fields::Matern(distance.mat, smoothness = 3/2)

svd(K)$d

min(distance.mat[upper.tri(distance.mat)])
d <- distance.mat
dvec <- d[upper.tri(d)]
dvec <- dvec[dvec > 0]
range0 <- median(dvec)

K <- fields::Matern(d, smoothness = 3/2, range = range0)
tau2 <- 1e-3
K <- K + diag(tau2, nrow(K))
solve(K)[1:5,1:5]
Ki <- chol2inv(chol(K))
Ki[1:5,1:5]

e <- eigen(K, symmetric=TRUE, only.values=TRUE)$values
range(e)
kappa(K)
max(abs(Ki - t(Ki)))


distance.mat[1:10,1:10]
solve(kernel)[1:10,1:10]


e <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
range(e)
kappa(K)

spatial_kernel[upper.tri(spatial_kernel)] = t(spatial_kernel)[upper.tri(spatial_kernel)]

L = chol(K)

solve(L)[1:10,1:10]


mask_train_valid_split <- function(
    obs_mask,
    valid_p = 0.2,
    seed = NULL,
    min_train_row = 1,
    min_train_col = 1
) {
  stopifnot(is.matrix(obs_mask))
  stopifnot(is.numeric(valid_p), valid_p >= 0, valid_p <= 1)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  obs_mask <- obs_mask != 0
  n_rows <- nrow(obs_mask)
  n_cols <- ncol(obs_mask)

  obs_ij <- which(obs_mask, arr.ind = TRUE)
  n_obs <- nrow(obs_ij)
  if (n_obs == 0) {
    return(matrix(0L, nrow = n_rows, ncol = n_cols))
  }

  target_valid <- floor(n_obs * valid_p)

  row_left <- rowSums(obs_mask)
  col_left <- colSums(obs_mask)

  # Random order over observed cells
  ord <- sample.int(n_obs)
  obs_ij <- obs_ij[ord, , drop = FALSE]

  valid_linear <- integer(0)

  for (k in seq_len(n_obs)) {
    if (length(valid_linear) >= target_valid) break

    r <- obs_ij[k, 1]
    c <- obs_ij[k, 2]

    # Only hold out if the remaining training counts stay >= minimums
    if ((row_left[r] - 1L) >= min_train_row && (col_left[c] - 1L) >= min_train_col) {
      valid_linear <- c(valid_linear, r + (c - 1L) * n_rows)
      row_left[r] <- row_left[r] - 1L
      col_left[c] <- col_left[c] - 1L
    }
  }

  if (length(valid_linear) < target_valid) {
    warning(
      sprintf(
        "Could only select %d/%d validation cells (valid_p=%.3f) under min_train_row=%d, min_train_col=%d.",
        length(valid_linear), target_valid, valid_p, min_train_row, min_train_col
      ),
      call. = FALSE
    )
  }

  valid_mask <- matrix(0L, nrow = n_rows, ncol = n_cols)
  valid_mask[valid_linear] <- 1L
  valid_mask
}

obs <- as.matrix(dat$modd$obs_mask)
valid1 = mask_train_valid_split(as.matrix(dat$modd$obs_mask), 0.05, 2025, 5, 5)
apply( obs * (1-valid2), 1, sum) %>% summary()
apply( valid2, 1, sum) %>% summary()

valid2 = IMR::mask_train_test_split(as.matrix(dat$modd$obs_mask), 0.05, 2025)

