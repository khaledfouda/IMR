require(tidyverse)
devtools::load_all()
require(magrittr)

#===== part 1: loading and preparing the data ===============
#===== out:   X, Z, Y, query, test.idx, test.truths, obs_mask
#=============================================================

load("article/movielens/data/Movie_X.Rdata") #X
load("article/movielens/data/Movie_Y.Rdata",verbose = T)
X <- X[,1:4] # keep only main-effects
input_tag = "_c_0_"
#Y <-     readRDS(paste0("article/movielens/data/Movie_Y",input_tag,".Rdata"))
query <- readRDS(paste0("article/movielens/data/Movie_Q",input_tag,".Rdata"))
source("other_models/Ma25.R")
source("other_models/SoftImpute_cv.R")
source("article/movielens/preprocess.R")
source("article/movielens/Ma25_fit.R")
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
hpar_int <- IMR::get_imr_default_hparams()
hpar_int$beta$value = hpar_int$beta$max = 0; hpar_int$beta$length = 1
hpar_int$gamma$value = hpar_int$gamma$max = 0; hpar_int$gamma$length = 1
hpar_int$laplace$max = 10; hpar_int$laplace$min = 10; hpar$laplace$step_sizes = c(1)
hpar_int$laplace$step_sizes = c(1, .1)


hparam <- IMR::get_imr_default_hparams()
hparam$beta$length <- 20
hparam$gamma$length <- 20
hparam$laplace$max <- 20
hparam$laplace$min <- 10
hparam$laplace$step_sizes <- c(1, .1)
hparam$beta$max = 1.8 #4.2 #.706
hparam$gamma$max = .6 #12 # 1.55

hparam$rank$default = 12
hparam$gamma$value = .5

#IMR:::initialize_parallel_workers(9)
seed=2025
run_training = FALSE
if(run_training){

# fit all model variations: [imr+x, imr+xz, imr+intercept, ma, si]
#----------------------
# 1. intercept only >
dat <- IMR::prepare_data(Y, X, Z,val_prop =  0.2, seed = 2025)

# get maximum lambda_beta
hpar$beta$max <- IMR::get_lambda_lasso_max(
  y_train = data$y_train,
  X = data$Xq,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  maxit = 50,
  thresh = 1e-3,
  init_maxit = 100,
  shared_information = TRUE,
  init_thresh = 1e-4,
  r = 5,
  verbose = 100
)
hpar$gamma$max <- IMR::get_lambda_lasso_max(
  y_train = data$y_train,
  Z = data$Zq,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  maxit = 50,
  thresh = 1e-3,
  shared_information = TRUE,
  init_maxit = 100,
  init_thresh = 1e-4,
  r = 5,
  verbose = 100
)

# models to test:
# 1. intercept
# 2. X
# 3. XZ
# 4. X shared
# 5. XZ shared [ need to find an upper-bound on these]

bench::bench_time(fit.imr1 <- IMR::imr.cv_laplace(dat, intercept_row = TRUE, intercept_col = TRUE,
                             hpar = hpar_int, thresh = 1e-3, maxit = 300,
                             trace = 1, ls_initial = TRUE, shared_information = FALSE,
                             seed = seed, num_cores = 9,final_fit = FALSE,
                             final_thresh = 1e-6, final_maxit = 1000)) -> time.imr

bench::bench_time(fit.imr2 <- IMR:::imr.cv_21(dat, intercept_row = TRUE, intercept_col = TRUE,
                                                  hpar = hparam, thresh = 1e-3, maxit = 600,
                                                  trace = 1, ls_initial = TRUE,
                                              shared_information = FALSE,
                                                  seed = seed, num_cores = 9,
                                                  final_thresh = 1e-6, separate = TRUE,
                                              final_maxit = 1000)) -> time.imr


fit.imr2$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)
saveRDS(fit.imr2, paste0("article/movielens/data/saved_models/",
                         "IMR_fit_Feb2.rds"))






fit.imr2$fit$params
fit.imr2$best_fit$lambda_laplace

#-- tmp start
time.imr
i=1
out = list(IMR::reconstruct(fit.imr2$fit, dat,shared_information = F))
fits = list(fit.imr2)
prepare_output_movielens(
  "IMR",
  time = fits[[i]]$time,
  X = dat$X,
  Z = dat$Z,
  beta.estim  = out[[i]]$beta,
  gamma.estim = out[[i]]$gamma,
  estim.test = out[[i]]$estimates[test.idx],
  estim.train = as.Incomplete(out[[i]]$estimates * dat$obs_mask)@x,
  obs.test = test.truths,
  obs.train = dat$Y[dat$Y!=0],
  M.estim = out[[i]]$M,test_error = IMR:::error_metric$rmse,
  rank.M = fits[[i]]$fit$params$rank
)
fit.imr1$best_fit$lambda_laplace
fit.imr1$best_fit$r
#-- tmp end
bench::bench_time(fit.imr1 <- IMR:::imr.cv_3(
  dat,
  intercept_row = T,
  intercept_col = T,
  hpar = hpar_int,
  trace = 1,
  num_cores = 7,
  seed = 2025
)) -> time.imr

fit.imr1$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr1, paste0("article/movielens/data/saved_models/",
                        "IMR_fit_intercept.rds"))

# 2. row covariates >
dat <- IMR::prepare_data(Y, X, NULL, 0.2, seed = 2025)
bench::bench_time(fit.imr2 <- IMR:::imr.cv(
  dat$model,
  intercept_row = T,
  intercept_col = T,
  fast.cv = F,
  hpar = hpar,
  verbose = 1,
  seed = 2025
)) -> time.imr

fit.imr2$time <- round(lubridate::time_length(time.imr[2],  "minute"),2)

saveRDS(fit.imr2, paste0("article/movielens/data/saved_models/",
                        "IMR_fit_rows_nonsp.rds"))


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

saveRDS(fit.imr3, paste0("article/movielens/data/saved_models/",
                        "IMR_fit_rows_cols.rds"))

#--- soft-Impute
bench::bench_time(fit.si <-
                    simpute.cv(y_full = dat$model$Y,
                               y_train = dat$model$y_train,
                               y_valid = dat$model$y_valid,
                               seed = 2025,
                               lambda_max = 100,
                               n.lambda = hpar$M$n.lambda,
                               trace = T)) -> time.si
fit.si$time <- round(lubridate::time_length(time.si[2],  "minute"),2)
saveRDS(fit.si, "article/movielens/data/saved_models/SI_fit.rds")
#------------ Ma
M = fit_MA25_movielens("", 2025)
}
#--------------------------------------------------------------
#==============================================================
#================ part 3 - load the models and show the results
#==============================================================
dat <- IMR::prepare_data(Y, X, Z, 0.2, seed = 2025)
fit.imr1 <- readRDS("article/movielens/data/saved_models/IMR_fit_intercept.rds")
fit.imr2 <- readRDS("article/movielens/data/saved_models/IMR_fit_rows_nonsp.rds")
fit.imr3 <- readRDS("article/movielens/data/saved_models/IMR_fit_rows_cols.rds")
fit.si   <- readRDS("article/movielens/data/saved_models/SI_fit.rds")
fit.ma   <- readRDS("article/movielens/data/saved_models/Ma_fit.rds")



out <- list()
out[[1]] <- IMR::reconstruct(fit.imr1$fit, dat)
out[[2]] <- IMR::reconstruct(fit.imr2$fit, dat)
out[[3]] <- IMR::reconstruct(fit.imr3$fit, dat)
out[[4]]   <- IMR::reconstruct(fit.si$fit, dat)
# out[[5]] <- IMR::reconstruct(fit.ma$fit, dat)
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
    M.estim = out[[i]]$M,test_error = IMR:::error_metric$rmse,
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
#===============================================================
require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")
#===== part 1: loading and preparing the data ===============

out <- prepare_results_for_analysis()
dat = out$dat; fits=out$fits; res.df=out$res



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
    format   = "html",
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
# data tables:
require(scales)
require(tidyr)

fmt_cell <- function(n, denom) sprintf("%s (%.1f%%)", comma(n), 100 * n / denom)

# X covariates:
X <- dat$X
n_users <- nrow(X)
gender <- factor(X[,1], levels = c(1,0), labels = c("Male", "Female"))


age <- factor(apply(X[,2:4], 1, function(x) {
  if(x[1] == 1){ "25-34"
  }else if(x[2] == 1){"35-49"
  }else if(x[3] == 1){"50+"
  }else "0-24"}))

users <- tibble(gender, age)
tab <- users |>
  count(gender, age, .drop = FALSE)

# Wide with formatted cells
tab_fmt <- tab |>
  mutate(cell = fmt_cell(n, n_users)) |>
  dplyr::select(-n) |>
  pivot_wider(
    names_from  = age,
    values_from = cell,
    values_fill = "0 (0.0%)"
  )

tot_by_gender <- users |>
  count(gender, name = "n") |>
  mutate(Total = fmt_cell(n, n_users)) |>
  dplyr::select(gender, Total)

tab_fmt <- left_join(tab_fmt, tot_by_gender, by = "gender")

# Totals row (by age) + grand total
tot_by_age <- users |>
  count(age, name = "n") |>
  mutate(val = fmt_cell(n, n_users)) |>
  dplyr::select(age, val) |>
  pivot_wider(names_from = age, values_from = val)

grand_total <- fmt_cell(n_users, n_users)
totals_row <- tot_by_age |>
  mutate(gender = "Total", Total = grand_total) |>
  relocate(gender)

# Bind totals row; order columns (age levels + Total)
col_order <- c(levels(age), "Total")
tab_final <- bind_rows(tab_fmt, totals_row) |>
  dplyr::select(gender, all_of(col_order))# %>%
# rename(`gender\\age group` = gender)


kbl(
  tab_final,
  format   = "latex",
  escape   = FALSE,
  align    = c("l", rep("r", 5)),
  caption  = paste0("User demographics in MovieLens-1M (Gender x Age group)")
) |>
  add_header_above(c(" " = 1, "Age group" = 5)) |>
  kable_styling(
    full_width        = FALSE,
    bootstrap_options = c("striped", "hover", "condensed"),
    font_size         = 12
  ) |>
  column_spec(1, width = "12em", bold = TRUE) |>
  column_spec(2:ncol(tab_final), width = "8em") |>
  row_spec(nrow(tab_final), bold = TRUE, background = "#f7f7f7") |>
  footnote(
    general = "Percentages are of all users.",
    threeparttable = TRUE
  )

# we do the same for Z
G <- dat$Z
m_titles <- nrow(G)
counts <- colSums(G > 0, na.rm = TRUE)

tibble(
  genre = colnames(G),
  n     = as.integer(counts)
) |>
  mutate(pct = 100 * n / m_titles) |>
  arrange(desc(n)) |>
  mutate(genre = factor(genre, levels = genre)) %>%
  arrange(desc(n)) %>%
  mutate(
    genre = fct_reorder(genre, n),                          # sort by count
    lbl   = paste0(comma(n), " (", round(pct, 1), "%)")
  ) %>%
  ggplot( aes(n, genre)) +
  geom_segment(aes(x = 0, xend = n, yend = genre), linewidth = 1, colour = "#B8C2CC") +
  geom_point(size = 3, colour = "#355C7D") +
  geom_text(aes(label = lbl), nudge_x = 0.02 * 1605, hjust = 0, size = 3.2) +
  scale_x_continuous("Number of movies", labels = comma, expand = expansion(mult = c(0, .1))) +
  labs(y = NULL,
       title = "Distribution of movie genres in MovieLens-1M",
       subtitle = paste("Per-title counts with percentages relative to 3,952 movies,",
                        "titles may carry multiple genres, so counts sum to more than 3,952.")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) ->pg;pg
ggsave("./article/movielens/data/plot_genres.png",
       pg, width = 320/25.4, height = 150/25.4, dpi = 600)




# grid size:

# Hyperparameter grids (MovieLens–1M) → LaTeX table with kableExtra
library(tibble)
library(dplyr)
library(kableExtra)

# Grid sizes (per your setup)
sz_imr_xz  <- 20*20 + 20*20
sz_imr_x   <- 20*20
sz_imr_int <- 20
sz_soft    <- 20

tab_hyper <- tribble(
  ~Model,           ~`Tuned hyperparameters`,                         ~`Grid structure`,                                                                 ~`Grid size`,                                      ~Notes,
  "IMR–Intercept",  "$\\lambda_{M}$",                                  "$\\lambda_{M}$ over 20 values (log–spaced)",                                     "$20$",                                           "$r$ chosen on validation; no $\\beta,\\gamma$ coefficients.",
  "IMR–X",          "$\\lambda_{\\beta},\\ \\lambda_{M}$",             "Cartesian grid $20\\times 20$",                                                   "$20\\times 20 = 400$",                           "$r$ chosen on validation.",
  "IMR–XZ",         "$\\lambda_{\\beta},\\ \\lambda_{\\gamma},\\ \\lambda_{M}\\ (+\\ r)$",
  "Blockwise grids: $(\\lambda_{\\beta},\\lambda_{M})$ and $(\\lambda_{\\gamma},\\lambda_{M})$, each $20\\times 20$",
  "$20\\times 20 + 20\\times 20 = 800$",           "$r$ chosen on validation; both $\\mathbf{X}$ and $\\mathbf{Z}$ active.",
  "Soft–Impute",    "$\\lambda,\\ r$",                                 "$\\lambda$ over 20 values (log–spaced)",                                          "$20$",                                           "Effective rank adapts with $\\lambda$; no side covariates.",
  "GLocal–K",       "---",                                             "---",                                                                              "---",                                            "Trained for 100 epochs (authors’ settings); no grid search."
)

# Build kable (LaTeX). Use escape = FALSE to render math.
kbl(
  tab_hyper,
  format   = "latex",
  booktabs = TRUE,
  escape   = FALSE,
  caption  = "Hyperparameter grids and search budget on MovieLens–1M. Counts report the number of $\\lambda$ configurations evaluated; ranks $r$ are selected on validation unless noted.",
  label    = "tab:movielens:hyperparameters",
  align    = c("l","l","l","r","l")
) |>
  kable_styling(latex_options = c("hold_position"))
# To save: save_kable(last_kable(), "tab-movielens-hyperparameters.tex")


#=========================================
#================== part 4 > insights ? =====================================
#############################################################################
fit <- out[[2]]
G <- fit$beta[1,]
which.max(G)
which.min(G)

#==============================================================================
#==== examine non-sparse coefficients in IMR =====================
out[[2]]$beta -> beta
Z <- t(Z)
`%||%` <- function(a, b) if (is.null(a)) b else a
#
# summarize_beta_by_genre(out[[2]]$beta, t(dat$Z)) %>%
#   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
#   arrange(covariate, mean) %>% view
#




genres_keep <- c("Action", "Musical", "Romance", "Sci-Fi", "Comedy")

# --- defensive names ------------------------------------------------------
covar_labs <- c("Sex (male)", "Age (25-34)", "Age (35-54)",
                "Age (55+)")
covar_names <- rownames(beta)
movie_ids   <-  as.character(seq_len(ncol(beta)))

if (!all(genres_keep %in% rownames(Z))) {
  stop("Some requested genres are missing from rownames(Z).")
}

# --- long beta ------------------------------------------------------------
df_beta <- data.frame(
  covariate = rep(covar_labs, times = ncol(beta)),
  movie     = rep(movie_ids,   each = nrow(beta)),
  value     = as.vector(beta),
  check.names = FALSE
)

# --- movie→genre map for the 5 genres ------------------------------------
Z_sel <- Z#Z[genres_keep, , drop = FALSE]
arr   <- which(Z_sel != 0, arr.ind = TRUE)

map <- data.frame(
  genre = rownames(Z_sel)[arr[, 1]],
  movie =  as.character(arr[, 2]),
  stringsAsFactors = FALSE
)

# --- join + clean ---------------------------------------------------------
plot_df <- df_beta %>%
  inner_join(map, by = "movie", relationship = "many-to-many") %>%
  # keep only movies belonging to the 5 genres
  filter(!is.na(value), value!=0) %>%
  mutate(
    # genre     = factor(genre, levels = genres_keep),
    genre     = factor(genre),
    covariate = factor(covariate, levels = covar_labs)  # keep your covariate order
  )

#-- >
X %>%
  as.data.frame() %>%
  mutate(group = case_when(
    .G == 1 & .A1 == 1 ~ "Male (25-34)",
    .G == 0 & .A1 == 1 ~ "Female (25-34)",
    .G == 1 & .A2 == 1 ~ "Male (35-54)",
    .G == 0 & .A2 == 1 ~ "Female (35-54)",
    .G == 1 & .A3 == 1 ~ "Male (55+)",
    .G == 0 & .A3 == 1 ~ "Female (55+)",
    .G == 1  ~ "Male (0-24)",
    .G == 0  ~ "Female (0-24)",
  )) -> Xg
out[[2]]$xbeta -> xbetahat
#groups <- unique(Xg$group)
groups <- c(
  #"Female (0-24)",
  "Female (25-34)",
  "Female (35-54)",
  "Female (55+)",
  #"Male (0-24)",
  "Male (25-34)",
  "Male (35-54)",
  "Male (55+)"
)


j <- 2
df.xbeta <- data.frame()
for(j in 1:length(groups)){

  uid <- which(Xg$group == groups[j])
  xbetau <- xbetahat[uid,]
  xbetau %>% dim()
  df.xbeta %<>% rbind( data.frame(
    user = rep(uid, times = ncol(xbetau)),
    movie = rep(movie_ids, each = nrow(xbetau)),
    value = as.vector(xbetau)
  ) %>%
    filter(value != 0) %>%
    mutate(group = groups[j]) %>%
    inner_join(map, by = "movie", relationship = "many-to-many")
  )
}

df.xbeta %<>%
  mutate(genre = factor(genre),
         group = factor(group, levels=groups)) %>%
  rename(covariate = group)
# df.xbeta[1:5,]
# dim(df.xbeta)
# plot_df[1:5,]
# map[1:5,]
# df_beta[1:5,]
# df
plot_df <- df.xbeta #%>%
#filter(genre %in% genres_keep)


#-- >>

dim(plot_df)
# --- boxplot --------------------------------------------------------------
# plot_df %>%
#   filter(genre %in% genres_keep) %>%
# ggplot(aes(genre, value, fill = genre)) +
#   geom_boxplot(width = 0.7, outlier.alpha = 0.2) +
#   facet_wrap(~ covariate, ncol = 2, scales = "free") +
#   # symmetric, log-like transform around 0
#   scale_y_continuous(
#   #  trans  = scales::asinh_trans(),     # or: scales::modulus_trans(p = 0.5)
#     breaks = scales::pretty_breaks(3)
#   ) +
#   scale_fill_brewer(palette = "Set2", guide = "none") +
#   labs(x = "Genre", y = expression(beta~"coefficient (asinh scale)")) +
#   theme_classic(base_size = 11) +
#   theme(strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 0, hjust = 1))


library(tidytext)
plot_df_ord <- plot_df %>%
  mutate(
    # order genres within each covariate facet by mean beta (lowest first)
    genre_ord = tidytext::reorder_within(genre,- value, within = covariate,
                                         fun = median)
  )

# ggplot(plot_df_ord, aes(x = genre_ord, y = value, fill = genre)) +
#   geom_boxplot(width = 0.7, outlier.alpha = 0.2) +
#   facet_wrap(~ covariate, ncol = 2, scales = "free") +
#   # symmetric transform to reveal structure while keeping signs
#   scale_y_continuous(
#     # trans = scales::asinh_trans(),
#                      breaks = scales::pretty_breaks(3),
#                      name = expression(beta~"coefficient (asinh scale)")) +
#   scale_x_reordered(name = "Genre (ordered by mean β within facet)") +
#   scale_fill_brewer(palette = "Set2", guide = "none") +
#   theme_classic(base_size = 11) +
#   theme(strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 0, hjust = 1))

plot_df %>%
  mutate(value = as.numeric(value)) %>%
  filter(covariate =="Age (25-34)") %>%
  dplyr::group_by(covariate, genre) %>%
  summarise(m = round(mean(value),3), mm = round(median(value),3)) %>%
  ungroup() %>%
  arrange(desc(mm))



stat_fun <- function(z) median(z, na.rm = TRUE)

by_cov_genre <- plot_df %>%
  filter(str_detect(covariate, "Female")) %>%
  filter(!(genre %in% c("Western", "Film-Noir"))) %>%
  group_by(covariate, genre) %>%
  summarise(stat = stat_fun(value), .groups = "drop")

extremes <- bind_rows(
  by_cov_genre %>% group_by(covariate) %>% slice_max(stat, n = 3, with_ties = FALSE),
  by_cov_genre %>% group_by(covariate) %>% slice_min(stat, n = 3, with_ties = FALSE)
) %>%
  distinct(covariate, genre, .keep_all = TRUE)   # avoid overlap if <6 genres


extreme_genres <- as.character(unique(extremes$genre))

require(scales)
# 2) Keep only those genres, and order x within each facet (DESC by mean β)
plot_df_ext <- plot_df %>%
  inner_join(extremes, by = c("covariate", "genre")) %>%   # adds 'stat'
  mutate(
    genre_ord = reorder_within(genre, -stat, within = covariate) # descending
  )

# 3) Plot (asinh y to tame outliers)
# ggplot(plot_df_ext, aes(x = genre_ord, y = value, fill = genre)) +
#   geom_boxplot(width = 0.7, outlier.alpha = 0.2) +
#   facet_wrap(~ covariate, ncol = 2, scales = "free") +   # different x-order per facet
#   scale_y_continuous(
#     trans  = asinh_trans(),
#     breaks = pretty_breaks(5),
#     name   = expression(beta~"coefficient (asinh scale)")
#   ) +
#   scale_x_reordered(name = "Genre (top 3 and bottom 3 by mean β, descending)") +
#   ggthemes::scale_color_tableau("Tableau 20") +
#   #scale_fill_brewer(palette = "Set2", guide = "none") +
#   theme_classic(base_size = 11) +
#   theme(
#     legend.position = "none",
#     strip.text = element_text(face = "bold"),
#     axis.text.x = element_text(angle = 25, hjust = 1)
#   )

###
library(dplyr)
library(ggplot2)
library(tidytext)  # reorder_within, scale_y_reordered (we'll flip axes)

# Assume you have one row per movie with columns:
#   group  (e.g., "Female (25-34)"),
#   genre  (e.g., "Comedy"),
#   value  (coefficient or predicted rating for that movie)
# If your column is named 'pred' replace value -> pred below.

# emoji_map_options <- list(
#   "Documentary" = c("📚", "🎥"),
#   "Musical"     = c("🎵", "🎼"),
#   "Film-Noir"   = c("🕵️", "🎞️"),
#   "Children's"  = c("🧒", "🧸"),
#   "War"         = c("🪖", "⚔️"),
#   "Western"     = c("🤠", "🌵"),
#   "Horror"      = c("👻", "🧛"),
#   "Action"      = c("💥", "🏃‍♂️"),
#   "Sci-Fi"      = c("🚀", "👽"),
#   "Animation"   = c("🎬", "🎨"),
#   "Fantasy"     = c("🧙‍♂️", "")
# )

emoji_map <- c(
  "Documentary"="📚","Musical"="🎵","Film-Noir"="🎞️","Children's"="🧸",
  "War"="⚔️","Western"="🤠","Horror"="👻","Action"="💥",
  "Sci-Fi"="🚀","Animation"="🎨","Fantasy"="🐉",
  "Romance" = "❤️", "Adventure" = "⛰️", "Drama" = "🎭",
  "Crime"  = "👮"
)
#
# emoji_map <- c(
#   "Documentary" = "books",           #
#   "Musical"     = "musical_note",    #
#   "Film-Noir"   = "film_frames",     #
#   "Children's"  = "teddy_bear",      #
#   "War"         = "crossed_swords",  #
#   "Western"     = "cowboy_hat_face", #
#   "Horror"      = "ghost",           #
#   "Action"      = "boom",            #
#   "Sci-Fi"      = "rocket",          #
#   "Animation"   = "art",             #
#   "Fantasy"     = "dragon"           #
# )


library(sysfonts)
library(ragg)
# library(showtext)

# library(emojifont)
# sysfonts::font_add("emoji", regular = "/System/Library/Fonts/Apple Color Emoji.ttc")

showtext::showtext_auto(FALSE)

plot_df %>%
  filter(str_detect(covariate, "Female")) %>%
  mutate(genre = as.character(genre)) %>%
  filter(genre %in% as.character(extreme_genres)) %>%
  mutate(emoji = emoji_map[as.character(genre)]) ->
  fplot_df


for(group in groups){

  gg = as.character(filter(extremes, covariate == group)$genre)
  fplot_df %<>%
    mutate(emoji = ifelse((!genre %in% gg) & covariate == group, "", emoji))%>%
    mutate(genre = ifelse((!genre %in% gg) & covariate == group, "", genre))
}

fplot_df %>%
  rename(group = covariate) %>%
  #filter(genre %in% genres_keep) %>%
  group_by(group, genre) %>%
  summarise(
    emoji = emoji[1],
    n     = n(),
    mean  = median(value, na.rm = TRUE),
    sd    = sd(value,   na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, pmax(n - 1, 1)),
    low   = mean - tcrit * se,
    high  = mean + tcrit * se,
    .groups = "drop"
  ) %>%
  # order genres within each facet by mean (descending)
  mutate(genre_ord = tidytext::reorder_within(genre, mean, within = group)) %>%


  ggplot( aes(x = mean, y = genre_ord)) +
  # geom_vline(xintercept = 0.5, colour = "grey75", linewidth = 0.4) +
  #geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.6, colour = "grey30") +
  #geom_point(size = 2.2, colour = "black") +
  #ggimage::geom_emoji(aes(emoji=emoji), size= 6) +
  geom_text(aes(label = genre), size=2) +
  # geom_text(aes(label = genre), size=3, family = "Apple Color Emoji") +
  facet_wrap(~ group, scales = "free_y") +
  scale_y_reordered(name = NULL)+
  labs(x = "",
       title = "Mean contribution to the predictor for each group, averaged over movies in selected genres",
       subtitle = "Genre selection was based on the top and bottom 3 genres for each group.") +
  theme_minimal() +
  #ggthemes::theme_few() +
  theme(strip.text = element_text(size = 9, face = "bold"),
        title = element_text(size=9),
        text = element_text(size = 8), panel.grid = element_blank(),
        axis.text.y = element_blank()) -> p; p

ragg::agg_png("article/fig.png", width = 2200, height = 1500, res = 300);print(p); dev.off()


# ggplot( aes(x = group, y = genre, fill = mean)) +
#   geom_tile() +
#    #geom_text(aes(label = paste0("n=", n)), colour = "grey15", size = 3) +
#    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
#                                                 midpoint = 0.5, name = "Mean\ncoefficient") +
#    labs(x = NULL, y = NULL, title = "Mean genre effect by group") +
#    theme_minimal(base_size = 11) +
#    theme(axis.text.x = element_text(angle = 0, hjust = 1),
#                panel.grid = element_blank())
#

#-----------------------------
#=== putting them in a table:
sumry <- plot_df %>%
  rename(group = covariate) %>%
  group_by(group, genre) %>%
  summarise(
    n     = n(),
    mean  = mean(value, na.rm = TRUE),
    sd    = sd(value,   na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, pmax(n - 1, 1)),
    low   = mean - tcrit * se,
    high  = mean + tcrit * se,
    .groups = "drop"
  )
# --- helpers -----------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "—", sprintf(paste0("%.", d, "f"), x))

# builds a comparison table Female vs Male for each age band
build_topbot_comp <- function(sumry, top_n = 3, bottom_n = 3) {
  df <- sumry %>%
    tidyr::extract(group, c("sex", "age"), "^(Female|Male)\\s*\\(([^)]+)\\)$",
                   remove = FALSE)

  # pick extremes per sex×age (top & bottom)
  extremes <- df %>%
    group_by(sex, age) %>%
    slice_max(mean, n = top_n, with_ties = FALSE) %>%
    bind_rows(slice_min(.,mean, n = bottom_n, with_ties = FALSE)) %>%
    ungroup()

  # focus genres per age = union of extremes across sexes
  focus <- extremes %>% group_by(age, genre) %>% summarise(.groups = "drop")

  # keep needed stats; compute se if only CI is available
  df2 <- df %>%
    semi_join(focus, by = c("age", "genre")) %>%
    mutate(
      se = if (!"se" %in% names(.)) (high - low) / (2 * qt(0.975, pmax(n - 1, 1))) else se
    ) %>%
    dplyr::select(age, sex, genre, n, mean, se)

  # ranks within each sex×age
  ranks <- df %>%
    group_by(age, sex) %>% arrange(desc(mean), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>% ungroup() %>%
    dplyr::select(age, sex, genre, rank)

  # female vs male side-by-side
  wide <- df2 %>%
    tidyr::pivot_wider(names_from = sex, values_from = c(n, mean, se), names_sep = ".") %>%
    # difference, CI for difference (conservative t on min(n)-1)
    mutate(
      delta    = mean.Female - mean.Male,
      se_delta = sqrt(se.Female^2 + se.Male^2),
      df_min   = pmax(pmin(n.Female, n.Male) - 1, 1),
      tcrit    = qt(0.975, df = df_min),
      lo_delta = delta - tcrit * se_delta,
      hi_delta = delta + tcrit * se_delta,
      sig      = ifelse(!is.na(lo_delta) & (lo_delta > 0 | hi_delta < 0), "*", "")
    ) %>%
    left_join(ranks %>% filter(sex == "Female") %>% rename(rank.Female = rank),
              by = c("age", "genre")) %>%
    left_join(ranks %>% filter(sex == "Male") %>% rename(rank.Male = rank),
              by = c("age", "genre")) %>%
    group_by(age) %>% arrange(desc(abs(delta)), .by_group = TRUE) %>% ungroup()
  wide
}

render_topbot_table <- function(tbl) {
  # color scale for Δ
  pal <- scales::col_numeric(c("#b2182b", "#f7f7f7", "#2166ac"),
                             domain = range(tbl$delta, na.rm = TRUE))
  disp <- tbl %>%
    transmute(
      Age = age,
      Genre = genre,
      n_F = n.Female, Mean_F = fmt(mean.Female),
      n_M = n.Male,   Mean_M = fmt(mean.Male),
      Delta = kableExtra::cell_spec(
        sprintf("%s%s", fmt(delta), sig),
        "latex", escape = FALSE,
        background = pal(delta)
      ),
      `95% CI (Δ)` = ifelse(is.na(lo_delta), "—",
                            sprintf("[%s, %s]", fmt(lo_delta), fmt(hi_delta))),
      Rank_F = rank.Female,
      Rank_M = rank.Male
    )

  kbl(disp, format = "html", booktabs = TRUE, linesep = "",
      escape = FALSE,
      col.names = c("Age", "Genre", "$n_F$", "$\\bar x_F$",
                    "$n_M$", "$\\bar x_M$", "$\\Delta = \\bar x_F-\\bar x_M$",
                    "95\\% CI", "Rank$_F$", "Rank$_M$"),
      caption = "Top & bottom genres per age band. Values are mean Xβ across movies with that genre. Δ>0 favors Female. Asterisk indicates 95\\% CI for Δ excludes 0."
  ) %>%
    kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 9) %>%
    collapse_rows(columns = 1, valign = "top")
}

# --- Build + render -----------------------------------------------------------
tab_comp <- build_topbot_comp(sumry, top_n = 3, bottom_n = 3)
render_topbot_table(tab_comp)






