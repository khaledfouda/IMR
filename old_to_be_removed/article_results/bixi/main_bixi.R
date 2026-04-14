library(tidyverse)
library(devtools)
load_all()
library(magrittr)
source("./article_results/bixi/preprocess_data.R")
source("./other_models/SoftImpute_cv.R")
source("./article_results/bixi/fit_bktr.R")
source("./article_results/bixi/helpers.R")

# the following line is needed only once. The data is already created.
# preprocess_bixi_data(.6, "25Sep", 2025)

# DELETE
# dat <- prepare_bixi_data(0.99, "25Sep")
# dat$splits$test@x %>% length()
# bktr.fit <- fit_BKTR_to_Bixi(0.6, "25Sep", 2025)
# BKTR_Bixi_Wrapper(dat, "25Sep", 2025, 0.6)


hpar <- IMR::get_imr_default_hparams()
hpar$beta$n.lambda <-  30 # 30
hpar$beta$lambda_max <- .1
hpar$gamma$lambda_max <- .1
hpar$gamma$n.lambda <-  30 # 30
hpar$M$n.lambda <-  20 # 80
hpar$M$rank.step <- 2
hpar$M$early.stopping <- 1

library(IMR)
IMR::initialize_parallel_workers(9)

# select the models to fit inside the loop below
models <- list(
  SI = 0,
  BKTR = 0,
  Intercept = 0,
  covariate = 0,
  similarity = 1
)

for (miss in c( .55, .65, .75, .85, .99, .6, .7, .8, .9, .95)){
  # run BKTR
  if (models$BKTR) {
    bktr.fit <- fit_BKTR_to_Bixi(miss, "25Sep", 2025)
  }
  kernels <- generate_similarity_bixi(miss, "25Sep")

  rand <- 2025
  # 1. generate the data and select covariates
  odat <- prepare_bixi_data(miss, "25Sep")
  x_vars <- c("x_humidity", "x_mean_temp_c")
  z_vars <- c(
    "z_walkscore", "z_len_minor_road", "z_num_restaurants",
    "z_num_pop", "z_num_bus_stations",
    "z_len_major_road",
    "z_capacity"
  )
  odat$X <- odat$X[, x_vars]
  odat$Z <- odat$Z[, z_vars]
  dat <- IMR::prepare_data(odat$Y.inc, odat$X, odat$Z, 0.2, 2025)
  if (models$SI) {
    bench::bench_time(fit.si <-
      simpute.cv(dat$model$Y,
        y_train = dat$model$y_train,
        y_valid = dat$model$y_valid,
        seed = 2025,
        rank.limit = 30,
        tol = 2,
        lambda_max = 10,
        n.lambda = 30,
        trace = T
      )) -> time.si
    fit.si$time <- round(lubridate::time_length(time.si[2], "minute"), 2)
    saveRDS(fit.si, paste0("./article_results/bixi/data/SI_", round(100 * miss), "_fit.rds"))
  }
  #  IMR + covariates
  if (models$covariate) {
    bench::bench_time(fit.imr <- IMR::imr.cv(dat$model,
      row_intercept = T,
      col_intercept = T,
      hpar = hpar,
      verbose = 1,
      separate_tuning = T,
      lambda_gamma_default = 999,
      fast.cv = F,
      seed = rand
    )) -> time.imr
    fit.imr$time <- round(lubridate::time_length(time.imr[2], "minute"), 2)
    saveRDS(fit.imr, paste0("./article_results/bixi/data/imr_", round(100 * miss), "_fit.rds"))
  }
  symmetrize <- function(x) (x + t(x))/2

  # IMR + covariates + similarity
  if (models$similarity) {
    old_error <- 999
    #dat <- IMR::prepare_data(odat$Y.inc, NULL, NULL, 0.2, 2025)
    dat$model$similarity_rows <- solve(kernels$temporal) %>% symmetrize()
    dat$model$similarity_cols <- solve(kernels$spatial) %>% symmetrize()
    hpar <- get_imr_default_hparams(dat$model$similarity_rows, dat$model$similarity_cols, 0, 0)
    hpar$M$n.lambda <- 20
    hpar$M$rank.step <- 1
    hpar$laplace$step_sizes <- c(5, 1, 0.1)
    hpar$laplace$start_value <- 30
    hpar$laplace$end_value <- 0
    hpar$rank$step_sizes <- c(2, 1)
    hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1

      hpar$laplacian_col <- IMR::decompose_symmetric_matrix(dat$model$similarity_cols)
      hpar$laplacian_row <- IMR::decompose_symmetric_matrix(dat$model$similarity_rows)
      for (rand in 2035:2035) {

        bench::bench_time(fit.imr <- IMR::imr.cv_laplace(dat$model,
          row_intercept = T,
          col_intercept = T,lambda_beta = 0, lambda_gamma=0,
          hpar = hpar,
          trace = T,
          seed = rand
        )) -> time.imr
        fit.imr$time <- round(lubridate::time_length(time.imr[2], "minute"), 2)
        fit.imr$rand <- rand

        new_error <- output_wrapper_bixi(fit.imr, dat, odat)$error.test
        print(paste(new_error))
        if (new_error < old_error) {
          old_error <- new_error
          saveRDS(fit.imr, paste0(
            "./article_results/bixi/data/imr_", round(100 * miss),
            "_similarity_fit.rds"
          ))
        }

    }
  }

  # -- intercept only model
  if (models$Intercept) {
    hpar$laplacian_col <- hpar$laplacian_col <- list()
    dat <- IMR::prepare_data(odat$Y.inc, NULL, NULL, 0.2, 2025)
    bench::bench_time(fit.imr <- IMR::imr.cv(dat$model,
      row_intercept = T,
      col_intercept = T,
      hpar = hpar,
      verbose = 0,
      separate_tuning = T,
      lambda_gamma_default = 999,
      fast.cv = F,
      seed = rand
    )) -> time.imr
    fit.imr$time <- round(lubridate::time_length(time.imr[2], "minute"), 2)
    saveRDS(fit.imr, paste0(
      "./article_results/bixi/data/imr_", round(100 * miss),
      "_intercept_fit.rds"
    ))
  }


  #--- end test
  # print(paste(new_error_imr, "-", rand))
  # if(new_error_imr < current_error_imr){
  # current_error_imr = new_error_imr
  # }
  # }
  #  current_error_imr <- 999

  # for (rand in 2025:2040) {


  #



  # -- intercept only model
  #   dat <- IMR::prepare_data(odat$Y.inc, NULL, NULL, 0.2, rand)
  #   bench::bench_time(fit.imr <- IMR::imr.cv(dat$model,
  #     row_intercept = T,
  #     col_intercept = T,
  #     hpar = hpar,
  #     verbose = 0,
  #     separate_tuning = T,
  #     lambda_gamma_default = 999,
  #     fast.cv = F,
  #     seed = rand
  #   )) -> time.imr
  #   fit.imr$time <- round(lubridate::time_length(time.imr[2], "minute"), 2)
  #
  #   out <- IMR:::reconstruct(fit.imr$fit, dat)
  #   prepare_output_bixi(
  #     NULL,
  #     X = dat$X,
  #     estim.test = out$estimates[as.matrix(odat$splits$test != 0)],
  #     estim.train = out$estimates[as.matrix(odat$Y.inc != 0)],
  #     obs.test = odat$splits$test@x,
  #     obs.train = odat$Y.inc@x,
  #     beta.estim = out$beta,
  #     M.estim = out$M,
  #     time_per_fit = fit.imr$time_per_fit,
  #     total_num_fits = fit.imr$total_num_fits
  #   )$error.test -> new_error_imr
  #   #--- end test
  #   print(paste(new_error_imr, "-", rand))
  #   if (new_error_imr < current_error_imr) {
  #     current_error_imr <- new_error_imr
  #
  #     saveRDS(fit.imr, paste0(
  #       "./article_results/bixi/data/imr_", round(100 * miss),
  #       "_intercept_fit.rds"
  #     ))
  #   }
  # }
}

total_res <- data.frame()
for (miss in c(.6, .7, .8, .9, .95, .55, .65, .75, .85, .99)) {
  if (miss == 0.6) total_res <- data.frame()
  dat <- prepare_bixi_data(miss, "25Sep")
  x_vars <- c("x_humidity", "x_mean_temp_c")
  z_vars <- c(
    "z_walkscore", "z_len_minor_road", "z_num_restaurants",
    "z_num_pop", "z_num_bus_stations",
    "z_len_major_road",
    "z_capacity"
  )
  dat$X <- dat$X[, x_vars]
  dat$Z <- dat$Z[, z_vars]
  mdat <- IMR::prepare_data(dat$Y.inc, dat$X, dat$Z, 0.2, 2025)
  #--- models >>
  # 1. bktr >
  res.bktr <- BKTR_Bixi_Wrapper(dat, "25Sep", 2025, miss)
  # 2. covariates >
  fit.imr <- readRDS(paste0("./article_results/bixi/data/imr_", round(100 * miss), "_fit.rds"))
  out <- IMR:::reconstruct(fit.imr$fit, mdat)
  prepare_output_bixi(
    fit.imr$time * 60,
    X = dat$X,
    estim.test = out$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test = dat$splits$test@x,
    obs.train = dat$Y.inc@x,
    beta.estim = out$beta,
    M.estim = out$M,
    time_per_fit = fit.imr$time_per_fit,
    total_num_fits = fit.imr$total_num_fits
  ) -> res.imr

  # 2.5 similarity >
  fit.imr <- readRDS(paste0(
    "./article_results/bixi/data/imr_", round(100 * miss),
    "_similarity_fit.rds"
  ))
  out <- IMR:::reconstruct(fit.imr$fit, mdat)
  prepare_output_bixi(
    fit.imr$time * 60,
    X = dat$X,
    estim.test = out$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test = dat$splits$test@x,
    obs.train = dat$Y.inc@x,
    beta.estim = out$beta,
    M.estim = out$M,
    time_per_fit = fit.imr$time_per_fit,
    total_num_fits = fit.imr$total_num_fits
  ) -> res.imr3
  # 3. intercept >
  fit.imr <- readRDS(paste0(
    "./article_results/bixi/data/imr_", round(100 * miss),
    "_intercept_fit.rds"
  ))
  out <- IMR:::reconstruct(fit.imr$fit, dat)
  prepare_output_bixi(
    fit.imr$time * 60,
    X = dat$X,
    estim.test = out$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test = dat$splits$test@x,
    obs.train = dat$Y.inc@x,
    beta.estim = out$beta,
    M.estim = out$M,
    time_per_fit = fit.imr$time_per_fit,
    total_num_fits = fit.imr$total_num_fits
  ) -> res.imr2
  #--- Simpute
  fit.si <- readRDS(paste0(
    "./article_results/bixi/data/SI_", round(100 * miss),
    "_fit.rds"
  ))
  out <- IMR:::reconstruct(fit.si$fit, dat)
  prepare_output_bixi(
    fit.si$time * 60,
    X = dat$X,
    estim.test = out$estimates[as.matrix(dat$splits$test != 0)],
    estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
    obs.test = dat$splits$test@x,
    obs.train = dat$Y.inc@x,
    beta.estim = out$beta,
    M.estim = out$M,
    time_per_fit = fit.imr$time_per_fit,
    total_num_fits = fit.imr$total_num_fits
  ) -> res.si
  #---- combine
  total_res <- rbind(
    total_res,
    res.imr[-10] %>%
      rbind(res.imr2[-10]) %>%
      rbind(res.imr3[-10]) %>%
      rbind(res.si[-10]) %>%
      rbind(res.bktr[-c(13, 1:3)]) %>%
      cbind(miss = miss, model = c(
        "IMR", "IMR-Intercept",
        "IMR-Similarity", "SoftImpute", "BKTR"
      ))
  )
}
# ==========================================================================

total_res %>%
  as.data.frame() %>%
  mutate(
    error.test = as.numeric(error.test),
    corr.test = as.numeric(corr.test),
    error.train = as.numeric(error.train),
    miss = as.numeric(miss),
    time = as.numeric(time),
    sparsity = as.numeric(sparsity)
  ) %>%
  # group_by(miss) %>%
  # mutate(diff = error.test[1] - error.test[4]) %>%
  # ungroup() %>%
  as.data.frame() %>%
  dplyr::arrange(miss, error.test) %>%
  dplyr::select(-time_per_fit, -total_num_fits) %>%
  filter(miss != .99) %>%
  mutate(miss = 1 - miss) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) ->
results.df

# plot
require(scales)
require(ggrepel)


cols <- c(
  "IMR" = "#56B4E9", "IMR-Intercept" = "#0072B2",
  "IMR-Similarity" = "#3393D1",
  "SoftImpute" = "#E69F00", "BKTR" = "#009E73"
)
# cols <- c("IMR"="#009E73", "SoftImpute"="#D55E00", "BKTR"="#0072B2")
# lts  <- c("IMR"="solid", "SoftImpute"="dashed", "BKTR"="dotdash")
# shp  <- c("IMR"=16, "SoftImpute"=17, "BKTR"=15)

results.df %>%
  mutate(time = (time / 60)) %>%
  mutate(
    model = factor(model, levels = c(
      "IMR", "IMR-Intercept",
      "IMR-Similarity",
      "SoftImpute", "BKTR"
    ))
  ) %>%
  dplyr::select(model, miss, time, error.test, error.train) %>%
  pivot_longer(
    cols = c(time, error.test, error.train),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("error.test", "error.train", "time"),
      labels = c("Test error", "Train error", "Time (minutes, log scale)")
    )
  ) %>%
  filter(model != "IMR-Intercept") %>%
  ggplot(aes(miss, value, color = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7, stroke = 0.2) +
  facet_wrap(~metric,
    scales = "free_y", ncol = 3,
    labeller = labeller(metric = label_value)
  ) +
  scale_x_continuous("Observation rate",
    breaks = seq(0.05, 0.45, by = 0.1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_color_manual(values = cols, name = "Model") +
  # scale_linetype_manual(values = lts, name = "Model") +
  # scale_shape_manual(values = shp, name = "Model") +
  labs(y = NULL) +
  # theme_minimal(base_size = 11) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.justification = "left",
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8))
    # plot.margin = margin(t = 6, r = 14, b = 6, l = 6)
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Time (minutes, log scale)" ~ scale_y_continuous(
        trans  = "log",
        # breaks = breaks_pretty(n = 12),
        # labels = function(b) round(b,2)
        labels = function(b) paste0("log(", number((b), 0.01), ")")
      )
    )
  ) -> g1;g1


#====================================================================
ggsave("./article_results/bixi/data/plot_bixi.png",
  g1,
  width = 320 / 25.4, height = 150 / 25.4, dpi = 600
)
# =============================================================
# added on Nov 24. Goal: investigate covariate effects with similarity model
missp = .6; timestamp = "25Sep"
results.df %>% filter(as.numeric(miss)==.2)
dat <- prepare_bixi_data(missp, "25Sep")
x_vars <- c("x_humidity", "x_mean_temp_c")
z_vars <- c(
  "z_walkscore", "z_len_minor_road", "z_num_restaurants",
  "z_num_pop", "z_num_bus_stations",
  "z_len_major_road",
  "z_capacity"
)
dat$X <- dat$X[, x_vars]
dat$Z <- dat$Z[, z_vars]
mdat <- IMR::prepare_data(dat$Y.inc, dat$X, dat$Z, 0.2, 2025)
fit.imr <- readRDS(paste0(
  "./article_results/bixi/data/imr_", round(100 * missp),
  "_similarity_fit.rds"
))
out <- IMR:::reconstruct(fit.imr$fit, mdat)
mean(out$beta!=0)
mean(out$gamma!=0)

out$gamma %>% dim()
gamma_matrix <- out$gamma

# Convert to long format for ggplot
gamma_df <- gamma_matrix |>
  as.data.frame() |>
  mutate(day = 1:n()) |>
  pivot_longer(
    cols = -day,
    names_to = "covariate",
    values_to = "coefficient"
  )

# 2. Plot
gamma_df |>
  ggplot(aes(x = day, y = covariate, fill = coefficient)) +
  geom_tile() +
  # Use a diverging scale: Red (negative), White (0), Blue (positive)
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0
  ) +
  labs(
    title = "Temporal Evolution of Covariate Coefficients",
    subtitle = "White areas indicate days where the coefficient was shrunk to exactly zero",
    x = "Day Index",
    y = NULL,
    fill = "Effect Size"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

#==================================
# lab_ln <- function(breaks, show_eq = FALSE, acc = 0.01) {
#   v <- exp(breaks)                                   # back-transform
#   base <- paste0("log(", number(v, accuracy = acc), ")")
#   if (show_eq) paste0(base, " = ", round(breaks, 1)) else base
# }

# ggplot( aes(x = miss, y = value, color = model, group = model)) +
#   geom_line(linewidth = 0.7) +
#   geom_point(size = 1.6) +
#   facet_wrap(~ metric, scales = "free_y", ncol = 3) +
#   scale_x_continuous(
#     name   = "Sparsity (fraction of zeros)",
#     breaks = seq(.05, 0.45, by = 0.05),
#     labels = percent_format(accuracy = 1)
#   ) +
#   scale_color_brewer(palette = "Dark2", name = "Model") +
#   labs(y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(
#     legend.position = "top",
#     legend.title = element_text(size = 10),
#     panel.grid.minor = element_blank(),
#     strip.text = element_text(face = "bold"),
#     axis.title.x = element_text(margin = margin(t = 8)))
# ================================
results.df %>% filter(model %in% c("IMR","IMR-Similarity"))
dat$X[1:5, ]
dat$Z[1:5, ]
results.df %>% filter(sparsity <= 0.8 & model == "IMR")
results.df %>% filter(miss == 0.15)
# plot to present the covariates.
dat <- prepare_bixi_data(0.6, "25Sep")
x_vars <- c("x_humidity", "x_mean_temp_c")
z_vars <- c(
  "z_walkscore", "z_len_minor_road", "z_num_restaurants",
  "z_num_pop", "z_num_bus_stations",
  "z_len_major_road",
  "z_capacity"
)
dat$X <- odat$X[, x_vars]
dat$Z <- odat$Z[, z_vars]
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

# ----------------------- helpers -------------------------------

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0),
      strip.text       = element_text(face = "bold")
    )
}

hist_with_density <- function(x, xlab, binwidth, color_bars = "#B9C6D6", color_lines = "#2F4858") {
  x <- as.numeric(x)
  n <- length(x)
  # y-scale for density overlay = density * n * binwidth
  y_max <- hist(x,
    breaks = seq(floor(min(x) / binwidth) * binwidth,
      ceiling(max(x) / binwidth) * binwidth,
      by = binwidth
    ),
    plot = FALSE
  )$counts |> max(0)
  qs <- quantile(x, probs = c(.25, .5, .75), names = FALSE)

  ggplot(data.frame(x = x), aes(x)) +
    geom_histogram(
      binwidth = binwidth, boundary = 0,
      fill = color_bars, color = "#4B4B4B"
    ) +
    stat_density(aes(y = after_stat(density) * n * binwidth),
      color = color_lines, linewidth = 0.7, kernel = "gaussian"
    ) +
    geom_vline(xintercept = qs[2], color = "#232323") +
    geom_vline(xintercept = qs[c(1, 3)], linetype = "dashed", color = "#232323") +
    annotate("text",
      x = qs[2], y = y_max, label = paste0("Med = ", round(qs[2], 1)),
      vjust = -0.4, size = 3.2
    ) +
    labs(x = xlab, y = "Days") +
    theme_paper()
}

timeseries_loess <- function(dates, y, ylab, color_line = "#2F4858", color_smooth = "#BC412B") {
  # Monthly ticks (align to month starts); keep first and last day visible
  x_min <- min(dates)
  x_max <- max(dates)
  monthly <- seq(as.Date(format(x_min, "%Y-%m-01")), as.Date(format(x_max, "%Y-%m-01")), by = "1 month")
  ggplot(tibble(date = dates, y = as.numeric(y)), aes(date, y)) +
    geom_line(color = color_line, linewidth = 0.4) +
    geom_smooth(
      method = "loess", se = FALSE, span = 0.35,
      color = color_smooth, linewidth = 0.9
    ) +
    scale_x_date(breaks = monthly, date_labels = "%b %d") +
    labs(x = NULL, y = ylab) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

fd_bins_or_fixed <- function(x, fallback_bins = 30) {
  x <- as.numeric(x)
  bw <- 2 * IQR(x) / (length(x)^(1 / 3))
  if (!is.finite(bw) || bw <= 0) {
    # convert desired bin count to an approximate binwidth
    (max(x) - min(x)) / fallback_bins
  } else {
    bw
  }
}

# -------------------- day-level (rows: X) ----------------------

hum_name <- "x_humidity"
temp_name <- "x_mean_temp_c"

stopifnot(is.matrix(dat$X), !is.null(colnames(dat$X)))
stopifnot(hum_name %in% colnames(dat$X), temp_name %in% colnames(dat$X))

n_days <- nrow(dat$X)
# Days are consecutive starting 2019-04-15 over 196 days
date_seq <- as.Date("2019-04-15") + seq_len(n_days) - 1L

hum <- as.numeric(dat$X[, hum_name])
mnt <- as.numeric(dat$X[, temp_name])

# Plots: histograms (counts) + time series (loess)
p_hum_hist <- hist_with_density(hum, xlab = "Relative humidity (%)", binwidth = 5)
p_hum_ts <- timeseries_loess(date_seq, hum, ylab = "Relative humidity ")

# Temperature: choose a sane binwidth if range is narrow
bw_temp <- max(1.5, fd_bins_or_fixed(mnt, fallback_bins = 30))
# Round to a nice increment (0.5°C)
bw_temp <- max(0.5, round(bw_temp / 0.5) * 0.5)

p_tmp_hist <- hist_with_density(mnt, xlab = "Mean temperature (°C)", binwidth = bw_temp)
p_tmp_ts <- timeseries_loess(date_seq, mnt, ylab = "Mean temperature ") +
  geom_hline(yintercept = 0, color = "#7F7F7F", linewidth = 0.3)

# Arrange 2×2 panel
# fig_bixi_daycovs <- (p_hum_hist + p_hum_ts) / (p_tmp_hist + p_tmp_ts) +
fig_bixi_daycovs <- (p_hum_ts) / (p_tmp_ts) +
  plot_annotation(
    title = "BIXI day-level covariates (2019-04-15 to 2019-10-27)",
    subtitle = "daily standardized values with LOESS smoothing.",
    theme = theme_paper(11)
  )
fig_bixi_daycovs

ggsave("./article_results/bixi/data/plot_bixi_cov_x.png",
  fig_bixi_daycovs,
  width = 320 / 25.4, height = 150 / 25.4, dpi = 600
)

#--- do to Z
colnames(dat$Z) <- map(colnames(dat$Z), function(x) str_split(x, "z_")[[1]][2]) %>% unlist()
z_df <- as_tibble(dat$Z, .name_repair = "minimal")
# If there are more than 7 columns, keep the first 7; adjust as needed
if (ncol(z_df) > 7) z_df <- z_df[, 1:7]

# Long form for faceting
z_long <- z_df |>
  mutate(.row = row_number()) |>
  pivot_longer(-.row, names_to = "feature", values_to = "value") |>
  dplyr::select(-.row)

# Per-feature summaries for median/quantiles
z_summ <- z_long |>
  summarise(
    med = median(value),
    q1  = quantile(value, 0.25),
    q3  = quantile(value, 0.75),
    .by = feature
  )

# Faceted mini-histograms (uniform bin count for comparability)
p_station_hist <- ggplot(z_long, aes(value)) +
  geom_histogram(bins = 30, fill = "#B9C6D6", color = "grey60") + # "#4B4B4B") +
  # geom_density(color = "#2F4858", linewidth = 0.7) +
  geom_vline(data = z_summ, aes(xintercept = med), color = "#232323") +
  # geom_vline(data = z_summ, aes(xintercept = q1), linetype = "dashed", color = "#232323") +
  # geom_vline(data = z_summ, aes(xintercept = q3), linetype = "dashed", color = "#232323") +
  facet_wrap(~feature, scales = "free", ncol = 4) +
  scale_x_continuous(NULL, breaks = pretty_breaks(3)) +
  scale_y_continuous(
    "Number of stations",
    breaks = pretty_breaks(3),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    title = "BIXI station-level covariates",
    subtitle = paste0(
      "Per-station histograms; ",
      "dark line = median (values standardized)."
    )
  ) +
  theme_paper() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 8),
    panel.spacing = unit(6, "pt"),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  )
p_station_hist
ggsave("./article_results/bixi/data/plot_bixi_cov_z.png",
  p_station_hist,
  width = 320 / 25.4, height = 150 / 25.4, dpi = 600
)
# ===========================================================
out <- IMR:::reconstruct(fit.imr$fit, dat)

prepare_output_bixi(
  time = NULL,
  X = dat$Z,
  estim.test = out$estimates[as.matrix(dat$splits$test != 0)],
  estim.train = out$estimates[as.matrix(dat$Y.inc != 0)],
  obs.test = dat$splits$test@x,
  obs.train = dat$Y.inc@x,
  beta.estim = out$gamma,
  M.estim = out$M,
  time_per_fit = fit.imr$time_per_fit,
  total_num_fits = fit.imr$total_num_fits
)





fit <- simpute.cv(
  y_train = dat$splits$train,
  y_valid = dat$splits$valid,
  # W_valid   = dat$masks$valid,
  y_full = dat$Y.inc,
  n.lambda = hpar$M$n.lambda,
  trace = FALSE,
  print.best = FALSE,
  tol = 5,
  thresh = 1e-6,
  rank.init = hpar$M$rank.init,
  rank.limit = hpar$M$rank.limit,
  rank.step = hpar$M$rank.step,
  maxit = 600,
  seed = NULL
)
fit <- IMR:::reconstruct(fit$fit, dat)
prepare_output_bixi(
  time = fit$time,
  X = NULL,
  estim.test = fit$estimates[as.matrix(dat$splits$test != 0)],
  estim.train = fit$estimates[as.matrix(dat$Y.inc != 0)],
  obs.test = dat$splits$test@x,
  obs.train = dat$Y.inc@x,
  M.estim = fit$estimates,
  total_num_fits = fit$total_num_fits,
  time_per_fit = fit$time_per_fit
)


# BKTR .0870
# rows .0877
# both .0876
# interc .0878
# none  .088
# Simpute .101


