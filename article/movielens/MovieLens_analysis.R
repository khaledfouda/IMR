data <- load_movielens1m()
seed = 2025

model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0.2)
print(model_data)

model_combn <- data.frame(
  row_covariates = c(F,T,T),
  col_covariates = c(F,F,T),
  intercepts = c(T,T,T)
)


imrI <- readRDS("article/movielens/data/saved_models/March_IMR_I_fit.rds")
imrIX <- readRDS("article/movielens/data/saved_models/March_IMR_IX_fit.rds")
imrIXZ <- readRDS("article/movielens/data/saved_models/March_IMR_IXZ_fit.rds")
fit.si   <- readRDS("article/movielens/data/saved_models/SI_fit.rds")
fit.ma   <- readRDS("article/movielens/data/saved_models/Ma_fit.rds")
imrI$rank_M <- imrI$meta$rank
imrIX$rank_M <- imrIX$meta$rank
imrIXZ$rank_M <- imrIXZ$meta$rank
#-----------------------------------------------------------------------

imrIXZ

out <- list()
out[[1]] <- IMR::reconstruct(imrI, model_data)
out[[2]] <- IMR::reconstruct(imrIX, model_data)
out[[3]] <- IMR::reconstruct(imrIXZ, model_data)
out[[4]]   <- IMR::reconstruct(structure(list(coefficients=fit.si$fit),class= "imr_fit"),
                               model_data)
out[[5]] <- list()


fit.ma$rank_M = fit.ma$fit$fit[[1]]$rank
ffit.ma <- fit.ma$fit$fit[[1]]
ffit.ma$M <- ffit.ma$L %*% t(ffit.ma$R)
ffit.ma$xbeta <- cbind(1, data$X) %*% t(ffit.ma$B)
ffit.ma$estimates <- ffit.ma$M + ffit.ma$xbeta
ffit.ma$beta <- t(ffit.ma$B[,-1])
fit.ma$fit <- out[[5]] <- ffit.ma

res <- list()
fits <- list(imrI, imrIX, imrIXZ, fit.si, fit.ma)
models <- c("IMR-Intercept", "IMR-X", "IMR-XZ", "SoftImpute", "Ma")



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
    estim.train = as.Incomplete(out[[i]]$estimates * data$obs_mask)@x,
    obs.test = data$test.truths,
    obs.train = model_data$Y[model_data$Y!=0],
    M.estim = out[[i]]$M,
    test_error = IMR:::error_metrics$rmse,
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

dat <- data
res.list <- res
res <- res.df
out <- list(dat=dat, fits=fits, res.df=res.df, res=res, out=out)
