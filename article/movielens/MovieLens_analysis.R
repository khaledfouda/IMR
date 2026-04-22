# first, we read the data and prepare it for the analysis (table and figure)
require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")


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
    estim.train = as_incomplete(out[[i]]$estimates * data$obs_mask)@x,
    obs.test = data$test.truths,
    obs.train = model_data$Y[model_data$Y != 0],
    M.estim = out[[i]]$M,
    test_error = get_metric("rmse"),
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
  format = "simple",
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

dat <- out_all$dat
fits <- out_all$fits
res_df <- out_all$res

# f$lambda_gamma
# f <- fit.imr3; f$fit <- NULL; f
# --- movie→genre map for the 5 genres ------------------------------------
genres <- c(
  "Documentary", "Musical", "Drama", "Fantasy", "Children's",
  "War", "Action", "Sci-Fi", "Horror", "Animation"
)

z_sel <- dat$Z[, genres, drop = FALSE]
arr <- which(z_sel != 0, arr.ind = TRUE)

map_z <- data.frame(
  genre = colnames(z_sel)[arr[, 2]],
  movie = as.character(arr[, 1]),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(movie = as.numeric(movie))

map_x <- as.data.frame(dat$X) |>
  dplyr::mutate(user = seq_len(nrow(dat$X))) |>
  dplyr::mutate(
    group = dplyr::case_when(
      .G == 1 & .A1 == 1 ~ "Male (25-34)",
      .G == 0 & .A1 == 1 ~ "Female (25-34)",
      .G == 1 & .A2 == 1 ~ "Male (35-49)",
      .G == 0 & .A2 == 1 ~ "Female (35-49)",
      .G == 1 & .A3 == 1 ~ "Male (50+)",
      .G == 0 & .A3 == 1 ~ "Female (50+)",
      .G == 1            ~ "Male (0-24)",
      .G == 0            ~ "Female (0-24)"
    )
  ) |>
  dplyr::select(user, group)

map_z <- map_z |>
  dplyr::group_by(movie) |>
  dplyr::mutate(n = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::filter(n == 1) |>
  dplyr::select(-n)

movies <- unique(arr[, 1])
sub_yh <- (out_all$out[[2]]$xbeta + out_all$out[[2]]$gammaz)[, movies]
sub_yh <- out_all$out[[2]]$estimates[, movies] - sub_yh
sub_yh <- out_all$out[[2]]$estimates[, movies]
colnames(sub_yh) <- movies

ume <- as.data.frame(sub_yh) |>
  tibble::rownames_to_column("row") |>
  tidyr::pivot_longer(-row, names_to = "col", values_to = "value") |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
  dplyr::rename(movie = col, user = row, estimate = value) |>
  dplyr::inner_join(map_z, by = "movie", relationship = "many-to-one") |>
  dplyr::inner_join(map_x, by = "user", relationship = "many-to-one")

ume |>
  # count(movie, genre, group, name="n_cell") |>
  dplyr::group_by(genre, movie, group) |>
  dplyr::mutate(median_estim = median(estimate)) |>
  head()

q7 <- function(x, p) {
  stats::quantile(x, probs = p, na.rm = TRUE, type = 7, names = FALSE)
}

pooled_iqr_tbl <- ume |>
  dplyr::summarise(
    q1 = q7(estimate, 0.25),
    q3 = q7(estimate, 0.75),
    pooled_iqr = q3 - q1,
    .by = c(movie, genre)
  ) |>
  dplyr::select(movie, genre, pooled_iqr)

median_gap_tbl <- ume |>
  dplyr::summarise(
    med = median(estimate, na.rm = TRUE),
    .by = c(movie, genre, group)
  ) |>
  dplyr::summarise(
    s_diff = diff(range(med)),
    .by = c(movie, genre)
  )

# get movie titles:
titles <- data.table::fread(
  file = "article_results/movielens/data/movies_Z.dat",
  sep = NULL,
  encoding = "Latin-1",
  header = FALSE
) |>
  tidyr::separate(
    V1,
    into = c("movie", "title", "genres"),
    sep = "::"
  ) |>
  as.data.frame() |>
  dplyr::mutate(movie = as.numeric(movie)) |>
  dplyr::select(-genres) |>
  dplyr::filter(movie %in% map_z$movie)

selected_movies <- pooled_iqr_tbl |>
  dplyr::inner_join(median_gap_tbl, by = c("movie", "genre")) |>
  dplyr::mutate(
    .by = genre,
    z_iqr = as.numeric(scale(pooled_iqr)),
    z_diff = as.numeric(scale(s_diff)),
    s = 0.5 * z_diff + 0.5 * z_iqr
  ) |>
  dplyr::arrange(genre, dplyr::desc(s), dplyr::desc(s_diff), dplyr::desc(pooled_iqr), movie) |>
  dplyr::group_by(genre) |>
  dplyr::slice_max(order_by = s, n = 3, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::filter(genre %in% c("Children's")) |>
  dplyr::arrange(genre, dplyr::desc(s))

movies_picked <- c(586, 1367, 1592)
# c("Home Alone (1990)","101 Dalmatians (1996)", "Air Bud (1997)")

selected_movies <- pooled_iqr_tbl |>
  dplyr::inner_join(median_gap_tbl, by = c("movie", "genre")) |>
  dplyr::mutate(
    .by = genre,
    z_iqr = as.numeric(scale(pooled_iqr)),
    z_diff = as.numeric(scale(s_diff)),
    s = 0.5 * z_diff + 0.5 * z_iqr
  ) |>
  dplyr::arrange(genre, dplyr::desc(s), dplyr::desc(s_diff), dplyr::desc(pooled_iqr), movie) |>
  dplyr::filter(movie %in% movies_picked) |>
  dplyr::filter(genre %in% c("Children's")) |>
  dplyr::arrange(genre, dplyr::desc(s))

selected_movies <- titles |>
  dplyr::filter(movie %in% selected_movies$movie) |>
  dplyr::inner_join(selected_movies, by = "movie", relationship = "one-to-one") |>
  dplyr::select(-z_iqr, -z_diff)

ume_sel <- ume |>
  dplyr::semi_join(selected_movies, by = c("movie", "genre")) |>
  dplyr::inner_join(titles, by = "movie", relationship = "many-to-one")

gender_lv <- c("Female", "Male")
age_lv <- c("0-24", "25-34", "35-49", "50+")
# group_lv  <- as.vector(outer(gender_lv, age_lv, \(g,a) paste0(g, " (", a, ")")))
group_lv <- c(
  "Female (0-24)",
  "Female (25-34)",
  "Female (35-49)",
  "Female (50+)",
  "Male (0-24)",
  "Male (25-34)",
  "Male (35-49)",
  "Male (50+)"
)

ume_sel <- ume_sel |>
  dplyr::mutate(
    group = factor(group, levels = group_lv),
    age = stringr::str_extract(group, "(?<=\\().+(?=\\))") |> factor(levels = age_lv),
    gender = ifelse(stringr::str_detect(group, "^Female.*"), "Female", "Male") |> factor()
  )

# Order movies within genre by selection score (if available)
if (exists("selected_movies")) {
  movie_order <- selected_movies |>
    dplyr::arrange(genre, dplyr::desc(s)) |>
    dplyr::distinct(genre, movie) |>
    dplyr::group_by(genre) |>
    dplyr::mutate(order = dplyr::row_number()) |>
    dplyr::ungroup()

  ume_sel <- ume_sel |>
    dplyr::left_join(movie_order, by = c("movie", "genre"))
}

ume_sel$movie <- factor(
  ume_sel$movie,
  levels = unique(ume_sel$movie[order(ume_sel$genre, ume_sel$order)])
)

library(ggh4x)
library(grid)

g <- ggplot(ume_sel, aes(x = age, y = estimate, fill = age)) +
  geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 0.7, alpha = 0.9) +
  ggh4x::facet_nested(
    cols = vars(title, gender),
    scales = "fixed",
    # switch = "x",
    strip = ggh4x::strip_nested(
      text_x = ggh4x::elem_list_text(face = c("bold"))
    )
  ) +
  scale_fill_brewer(palette = "Set2", name = "Age") +
  scale_y_continuous(
    "Estimated rating",
    limits = c(0.5, 5),
    breaks = seq(0.5, 5, 0.5),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_discrete(name = "Age group") +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background.x = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    # axis.text.x = element_text(angle = 0, hjust = 0),
    legend.position = "none",
    ggh4x.facet.nestline = element_line(colour = "grey70")
  ) +
  ggtitle(
    "Fitted Movie Ratings (Full Model)",
    subtitle = "Selected children's movies by gender and age group"
  )

print(g)

# ggsave("./article_results/movielens/data/plot_intercept_model.png",
#        g, width = 320/25.4, height = 150/25.4, dpi = 600)

ggsave(
  filename = "./article/movielens/data/plot_full_model.png",
  plot = g,
  width = 9,
  height = 4,
  scale = 1.2,
  dpi = 300
)
