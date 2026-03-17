# first, we read the data and prepare it for the analysis (table and figure)
data <- load_movielens1m()
seed <- 2025
model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0)
print(model_data)

imr_i <- readRDS("article/movielens/data/saved_models/March_IMR_I_fit_1e6.rds")
imr_ixz <- readRDS("article/movielens/data/saved_models/March_IMR_IXZ_fit_1e6.rds")
fit_si <- readRDS("article/movielens/data/saved_models/SI_fit.rds")
fit_ma <- readRDS("article/movielens/data/saved_models/Ma_fit.rds")

imr_i$rank_M <- imr_i$meta$rank
imr_ixz$rank_M <- imr_ixz$meta$rank
imr_i$time <- imr_i$time_secs / 60
imr_ixz$time <- imr_ixz$time_secs / 60

out <- list()
out[[1]] <- IMR::reconstruct(imr_i, model_data)
out[[2]] <- IMR::reconstruct(imr_ixz, model_data)

fit_si_imr <- structure(list(coefficients = fit_si$fit), class = "imr_fit")
out[[3]] <- IMR::reconstruct(fit_si_imr, model_data)

out[[4]] <- list()
fit_ma$rank_M <- fit_ma$fit$fit[[1]]$rank
ffit_ma <- fit_ma$fit$fit[[1]]
ffit_ma$M <- ffit_ma$L %*% t(ffit_ma$R)
ffit_ma$xbeta <- cbind(1, data$X) %*% t(ffit_ma$B)
ffit_ma$estimates <- ffit_ma$M + ffit_ma$xbeta
ffit_ma$beta <- t(ffit_ma$B[, -1])

out[[4]] <- ffit_ma
fit_ma$fit <- ffit_ma

res <- list()
fits <- list(imr_i, imr_ixz, fit_si, fit_ma)
models <- c("IMR-I", "IMR-IXZ", "SoftImpute", "Ma")

for (i in 1:4) {
  print(models[[i]])
  res[[i]] <- prepare_output_movielens(
    models[[i]],
    time = fits[[i]]$time,
    X = data$X,
    Z = data$Z,
    beta.estim = out[[i]]$beta,
    gamma.estim = out[[i]]$gamma,
    estim.test = out[[i]]$estimates[data$test.idx],
    estim.train = as.Incomplete(out[[i]]$estimates * data$obs_mask)@x,
    obs.test = data$test.truths,
    obs.train = model_data$Y[model_data$Y != 0],
    M.estim = out[[i]]$M,
    test_error = IMR:::error_metrics$rmse,
    rank.M = fits[[i]]$rank_M
  )
}

res[[4]]$rank_beta <- 5
res[[5]] <- list(
  model = "Glocal-K",
  time = 52.36,
  error.test = 0.8516,
  corr.test = 0.6278,
  error.train = 0.7018,
  rank_M = 63,
  rank_beta = NA,
  rank_gamma = NA,
  sparsity_beta = NA,
  sparsity_gamma = NA
)

# 1- results table:
res_df <- do.call(rbind, lapply(res, function(x) if (length(x) == 14) x[c(1:2, 5:12)] else x)) |>
  as.data.frame() |>
  dplyr::mutate(
    dplyr::across(
      c(time, error.test, corr.test, error.train, sparsity_beta, sparsity_gamma, rank_M, rank_beta, rank_gamma),
      as.numeric
    )
  ) |>
  dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 3))) |>
  dplyr::mutate(
    rank_total = rank_M +
      ifelse(is.na(rank_beta), 0, rank_beta) +
      ifelse(is.na(rank_gamma), 0, rank_gamma)
  ) |>
  dplyr::arrange(error.test)

out_all <- list(dat = data, fits = fits, res = res_df, res_list = res, out = out)

#-----------------------------------------------------------------------
# we begin with the table:
dat <- out_all$dat
fits <- out_all$fits
res_df <- out_all$res

best_min_cols <- c("time", "error.test", "error.train", "rank_M")
best_max_cols <- c("corr.test")

best_idx_min <- lapply(best_min_cols, function(v) {
  which.min(replace(res_df[[v]], is.na(res_df[[v]]), Inf))
}) |>
  rlang::set_names(best_min_cols)

best_idx_max <- lapply(best_max_cols, function(v) {
  which.max(replace(res_df[[v]], is.na(res_df[[v]]) | is.nan(res_df[[v]]), -Inf))
}) |>
  rlang::set_names(best_max_cols)

fmt_num <- function(x, digits) {
  ifelse(is.na(x) | is.nan(x), "—", sprintf(paste0("%.", digits, "f"), x))
}

disp <- res_df |>
  dplyr::mutate(
    time = fmt_num(time, 2),
    error.test = fmt_num(error.test, 3),
    corr.test = fmt_num(corr.test, 3),
    error.train = fmt_num(error.train, 3),
    sparsity_beta = fmt_num(sparsity_beta, 3),
    sparsity_gamma = fmt_num(sparsity_gamma, 3),
    rank_M = ifelse(is.na(rank_M), "—", as.character(rank_M)),
    rank_beta = ifelse(is.na(rank_beta), "—", as.character(rank_beta)),
    rank_gamma = ifelse(is.na(rank_gamma), "—", as.character(rank_gamma)),
    rank_total = as.character(rank_total)
  ) |>
  dplyr::select(-rank_total)

for (v in names(best_idx_min)) {
  idx <- best_idx_min[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}

for (v in names(best_idx_max)) {
  idx <- best_idx_max[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}

library(kableExtra)
col_names <- c(
  "Model", "Time (min)",
  "RMSE", "correlation", "RMSE",
  "$M$", "$\\beta$", "$\\Gamma$",
  "$\\beta$", "$\\Gamma$"
)

kbl(
  disp,
  format = "html",
  booktabs = TRUE,
  linesep = "",
  escape = FALSE,
  col.names = col_names,
  caption = paste0(
    "Performance comparison on the MovieLens 1M dataset. ",
    "Best values per column are bolded and IMR models are shaded."
  )
) |>
  add_header_above(c(" " = 2, "Test" = 2, "Train" = 1, " " = 5)) |>
  add_header_above(c(" " = 2, "Performance" = 3, "Rank estimation" = 3, "Sparsity" = 2)) |>
  kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8) |>
  row_spec(which(grepl("^IMR", disp$model)), bold = FALSE, background = "#f7f7f7") |>
  column_spec(1, width = "3.2cm")

#===========================================================================================
# we now create the plot:

