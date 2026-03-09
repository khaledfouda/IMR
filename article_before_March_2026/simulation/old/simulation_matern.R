library(devtools)
# clean_dll()
# Rcpp::compileAttributes()
# document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article_results/simulation/generate_simu_dat.R")
# ===============================================================================
# you need the first function in test_cross_validation.R
kernel_type <- "AR1"
filename <- paste0("./article_results/saved/simulation_", kernel_type, "_res")

n <- m <- 800
r <- 5
p <- 0.8

seed <- 2025
B <- 15
run_all <- F
if (run_all) {
  all_res <- data.frame()
  for (r in c(5, 15)) {
    res <- data.frame()
    for (B in 1:15) {
      seed <- seed + 1
      set.seed(seed)
      dat <-
        generate_simulated_data(n, m, r, 0, 0, p,
          sparsity_beta = .5, sparsity_gamma = 0.0,
          structured_error_A = T,
          structured_error_B = T,
          prepare_for_fitting = T, mv_coeffs = T, seed = seed,
          kernel_type = kernel_type
        )

      dat$similarity_rows %<>% solve()
      dat$similarity_cols %<>% solve()

      inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
      data <- inp.dat$model
      hpar <- get_imr_default_hparams(dat$similarity_rows, data$similarity_cols, 0, 0)
      hpar$M$n.lambda <- 20
      hpar$M$rank.step <- 1
      hpar$laplace$step_sizes <- c(5, 1, 0.1)
      hpar$laplace$start_value <- 30
      hpar$laplace$end_value <- 0
      hpar$rank$step_sizes <- c(2, 1)
      hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1


      timesim <- system.time(
        fitsim <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar,
          trace = F, maxit = 600,
          seed = seed
        )
      )
      timesi <- system.time(fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
        trace = T,
        tol = 2, seed = seed,
        n.lambda = 20,
        test_error = IMR:::error_metric$rmse,
        rank.init = hpar$M$rank.min, maxit = 600,
        rank.step = 1, rank.limit = hpar$M$rank.max
        # lambda0_fun = IMR::get_lambda_m_max
      ))


      timer <- system.time(fit.imrS <- IMR::imr.cv(inp.dat$model,
        intercept_row = F,
        seed = seed, ls_initial = FALSE, hpar = hpar,
        maxit = 600,
        intercept_col = F, verbose = 0
      ))

      res <- bind_rows(res, rbind(
        sim_res(dat, fitsi$fit, "SoftImpute") %>%
          dplyr::mutate(valerr = fitsi$error, time = timesi[[3]]),
        sim_res(dat, fit.imrS$fit, "IMR") %>%
          dplyr::mutate(valerr = fit.imrS$error, time = timer[[3]]),
        sim_res(dat, fitsim$fit, "IMR-Similarity") %>%
          dplyr::mutate(valerr = fitsim$cols$best_error, time = timesim[[3]])
      ))
      print(res)
    }
    res %<>% mutate(n = n, m = m, true_rank = r, pmiss = p)
    all_res %<>% rbind(res)
    saveRDS(all_res, paste0(filename, "_1.rds"))
  }
  # =========
  # in the following, we provide the true rank to the models


  all_res <- data.frame()
  for (r in c(5, 15)) {
    res <- data.frame()
    for (B in 1:15) {
      seed <- seed + 1
      set.seed(seed)
      dat <-
        generate_simulated_data(n, m, r, 0, 0, p,
          sparsity_beta = .5, sparsity_gamma = 0.0,
          structured_error_A = T,
          structured_error_B = T,
          prepare_for_fitting = T, mv_coeffs = T, seed = seed
        )

      dat$similarity_rows %<>% solve()
      dat$similarity_cols %<>% solve()

      inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
      data <- inp.dat$model
      hpar <- get_imr_default_hparams(dat$similarity_rows, data$similarity_cols, 0, 0)
      hpar$M$n.lambda <- 20
      hpar$M$rank.step <- 1
      hpar$rank$rank.min <- hpar$rank$rank.max <- hpar$M$rank.min <-
        hpar$M$rank.max <- hpar$M$rank.init <- r
      hpar$laplace$step_sizes <- c(5, 1, 0.1)
      hpar$laplace$start_value <- 30
      hpar$laplace$end_value <- 0
      hpar$rank$step_sizes <- c(2, 1)
      hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1


      timesim <- system.time(
        fitsim <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar,
          trace = F, maxit = 600,
          seed = seed
        )
      )
      timesi <- system.time(fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
        trace = T,
        tol = 2, seed = seed,
        n.lambda = 20,
        test_error = IMR:::error_metric$rmse,
        rank.init = hpar$M$rank.min, maxit = 600,
        rank.step = 1, rank.limit = hpar$M$rank.max
        # lambda0_fun = IMR::get_lambda_m_max
      ))


      timer <- system.time(fit.imrS <- IMR::imr.cv(inp.dat$model,
        intercept_row = F,
        seed = seed, ls_initial = FALSE, hpar = hpar,
        maxit = 600,
        intercept_col = F, verbose = 0
      ))

      res <- bind_rows(res, rbind(
        sim_res(dat, fitsi$fit, "SoftImpute") %>%
          dplyr::mutate(valerr = fitsi$error, time = timesi[[3]]),
        sim_res(dat, fit.imrS$fit, "IMR") %>%
          dplyr::mutate(valerr = fit.imrS$error, time = timer[[3]]),
        sim_res(dat, fitsim$fit, "IMR-Similarity") %>%
          dplyr::mutate(valerr = fitsim$cols$best_error, time = timesim[[3]])
      ))
      print(res)
    }
    res %<>% mutate(n = n, m = m, true_rank = r, pmiss = p)
    all_res %<>% rbind(res)
    saveRDS(all_res, paste0(filename, "_2.rds"))
  }
}
# =========== Done training
# --- analyzing results

res1 <- readRDS(paste0(filename, "_1.rds"))
res2 <- readRDS(paste0(filename, "_2.rds"))

res1 %>%
  group_by(model, n, m, true_rank, pmiss) %>%
  summarise(
    theta_rmse_mean = mean(theta_rmse),
    theta_rmse_sd = sd(theta_rmse),
    test_rmse_mean = mean(test_rmse),
    test_rmse_sd = sd(test_rmse),
    rank_mean = mean(rank),
    rank_sd = sd(rank),
    eigen_rmse_mean = mean(eigen_rmse),
    eigen_rmse_sd = sd(eigen_rmse),
    .groups = "drop"
  ) %>%
  arrange(true_rank, test_rmse_mean) %>%
  transmute(
    model,
    true_rank = sprintf("(%d,%d) with rank = %d", n, m, true_rank),
    `M RMSE` = sprintf("%.2f (%.2f)", theta_rmse_mean, theta_rmse_sd),
    `test RMSE` = sprintf("%.2f (%.2f)", test_rmse_mean, test_rmse_sd),
    `rank estimate` = sprintf("%.1f (%.1f)", rank_mean, rank_sd),
    `eigen RMSE` = sprintf("%.2f (%.2f)", eigen_rmse_mean, eigen_rmse_sd)
  ) %>%
  gt::gt(rowname_col = "model", groupname_col = "true_rank") %>%
  gt::tab_header(
    title = "Simulation M = AB^T where both A and B are generated from GP with Matern kernel",
    subtitle = "Mean (SD) over 15 replications. 80% of the data is missing (test set)"
  )


res2 %>%
  group_by(model, n, m, true_rank, pmiss) %>%
  summarise(
    theta_rmse_mean = mean(theta_rmse),
    theta_rmse_sd = sd(theta_rmse),
    test_rmse_mean = mean(test_rmse),
    test_rmse_sd = sd(test_rmse),
    rank_mean = mean(rank),
    rank_sd = sd(rank),
    eigen_rmse_mean = mean(eigen_rmse),
    eigen_rmse_sd = sd(eigen_rmse),
    .groups = "drop"
  ) %>%
  arrange(true_rank, test_rmse_mean) %>%
  transmute(
    model,
    true_rank = sprintf("(%d,%d) with rank = %d", n, m, true_rank),
    `M RMSE` = sprintf("%.2f (%.2f)", theta_rmse_mean, theta_rmse_sd),
    `test RMSE` = sprintf("%.2f (%.2f)", test_rmse_mean, test_rmse_sd),
    `rank estimate` = sprintf("%.1f (%.1f)", rank_mean, rank_sd),
    `eigen RMSE` = sprintf("%.2f (%.2f)", eigen_rmse_mean, eigen_rmse_sd)
  ) %>%
  gt::gt(rowname_col = "model", groupname_col = "true_rank") %>%
  gt::tab_header(
    title = "Simulation M = AB^T where both A and B are generated from GP with Matern kernel",
    subtitle = "Mean (SD) over 15 replications. 80% of the data is missing (test set)"
  )


# ==============================================================================
# Results 2 ------  one run - analyze results - show estimated eigen-values
n <- m <- 800
r <- 5
p <- 0.8
seed <- 2025
if (run_all) {
  shuffle_mat <- function(x) {
    x <- matrix(sample(as.vector(x)), nrow(x), ncol(x))
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
    diag(x) <- 1
    x
  }
  res <- data.frame()
  set.seed(seed)
  dat <-
    generate_simulated_data(n, m, r, 0, 0, p,
      sparsity_beta = .5, sparsity_gamma = 0.0,
      structured_error_A = T,
      structured_error_B = T, kernel_type = kernel_type,
      prepare_for_fitting = T, mv_coeffs = T, seed = seed
    )

  dat$similarity_rows %<>% solve()
  dat$similarity_cols %<>% solve()

  inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
  data <- inp.dat$model
  hpar <- get_imr_default_hparams(dat$similarity_rows, data$similarity_cols, 0, 0)
  hpar$M$n.lambda <- 20
  hpar$M$rank.step <- 1
  hpar$laplace$step_sizes <- c(5, 1, 0.1)
  hpar$laplace$start_value <- 30
  hpar$laplace$end_value <- 0
  hpar$rank$step_sizes <- c(2, 1)
  hpar$rank$n_streaks <- hpar$laplace$n_streaks <- 1


  timesim <- system.time(
    fitsim <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar,
      trace = F, maxit = 600,
      seed = seed
    )
  )

  # same model but with shuffled precision matrix
  inp.dat$model$similarity_rows <- shuffle_mat(inp.dat$model$similarity_rows)
  inp.dat$model$similarity_cols <- shuffle_mat(inp.dat$model$similarity_cols)
  hpar$laplacian_row <- IMR::decompose_symmetric_matrix(inp.dat$model$similarity_rows)
  hpar$laplacian_col <- IMR::decompose_symmetric_matrix(inp.dat$model$similarity_cols)

  timesimF <- system.time(
    fitsimF <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar,
      trace = F, maxit = 600,
      seed = seed
    )
  )


  # same model but with identity precision matrix
  inp.dat$model$similarity_rows <- diag(1, n, n)
  inp.dat$model$similarity_cols <- diag(1, m, m)
  hpar$laplacian_row <- IMR::decompose_symmetric_matrix(inp.dat$model$similarity_rows)
  hpar$laplacian_col <- IMR::decompose_symmetric_matrix(inp.dat$model$similarity_cols)

  timesimI <- system.time(
    fitsimI <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, F, F, hpar,
      trace = F, maxit = 600,
      seed = seed
    )
  )


  # soft-impute
  timesi <- system.time(fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
    trace = T,
    tol = 2, seed = seed,
    n.lambda = 20,
    test_error = IMR:::error_metric$rmse,
    rank.init = hpar$M$rank.min, maxit = 600,
    rank.step = 1, rank.limit = hpar$M$rank.max
    # lambda0_fun = IMR::get_lambda_m_max
  ))


  # original IMR
  timer <- system.time(fit.imrS <- IMR::imr.cv(inp.dat$model,
    intercept_row = F,
    seed = seed, ls_initial = FALSE, hpar = hpar,
    maxit = 600,
    intercept_col = F, verbose = 0
  ))

  res <- bind_rows(res, rbind(
    sim_res(dat, fitsi$fit, "SoftImpute") %>%
      dplyr::mutate(valerr = fitsi$error, time = timesi[[3]]),
    sim_res(dat, fit.imrS$fit, "IMR") %>%
      dplyr::mutate(valerr = fit.imrS$error, time = timer[[3]]),
    sim_res(dat, fitsim$fit, "IMR-Similarity") %>%
      dplyr::mutate(valerr = fitsim$cols$best_error, time = timesim[[3]]),
    sim_res(dat, fitsimF$fit, "IMR-Similarity-Shuffled") %>%
      dplyr::mutate(valerr = fitsimF$cols$best_error, time = timesimF[[3]]),
    sim_res(dat, fitsimI$fit, "IMR-Similarity-Identity") %>%
      dplyr::mutate(valerr = fitsimI$cols$best_error, time = timesimI[[3]])
  ))

  all_res <- list(
    IMR = fit.imrS, IMRS = fitsim,
    IMRSF = fitsimF, IMRSI = fitsimI,
    true_eig = IMR::svd_opt(dat$theta, tol = 1e-6)$d,
    dat = dat,
    fitSI = fitsi, res = res
  )
  print(res)
  saveRDS(all_res, paste0(filename, "_3.rds"))
}
#----------------------------
# we analyze here
match_length <- function(x, n) {
  if (length(x) > n) {
    x <- x[1:n]
  } else if (length(x) < n) {
    x <- c(x, rep(0, n - length(x)))
  }
  x
}
# put it in a table


all_res <- readRDS(paste0(filename, "_3.rds"))
fit.imrS <- all_res$IMR
fitsim <- all_res$IMRS
fitsimF <- all_res$IMRSF
fitsimI <- all_res$IMRSI
fitsi <- all_res$fitSI
res <- all_res$res
true_eig <- all_res$true_eig

data.frame(
  IMR = match_length(fit.imrS$fit$d, r) - true_eig,
  IMR_matern = match_length(fitsim$fit$d, r) - true_eig,
  IMR_shuffled = match_length(fitsimF$fit$d, r) - true_eig,
  IMR_identity = match_length(fitsimI$fit$d, r) - true_eig,
  SoftImpute = match_length(fitsi$fit$d, r) - true_eig,
  true = true_eig
) |>
  abs() |>
  round(1) |>
  mutate(eigenvalue = paste0("λ", row_number())) |>
  relocate(eigenvalue, true, .before = IMR) %>%
  gt::gt() |>
  gt::tab_header(
    title = "Absolute error in the estimated singular values",
    subtitle = "Only the first 5 values were considered"
  ) |>
  gt::cols_label(
    eigenvalue   = "Eigenvalue",
    true         = "true",
    IMR_matern   = "(Matern)",
    IMR_identity = "(identity)",
    IMR          = "(original)",
    SoftImpute   = "SoftImpute",
    IMR_shuffled = "(shuffled)",
  ) |>
  gt::fmt_number(
    columns = c(true, IMR_matern, IMR_identity, IMR, SoftImpute, IMR_shuffled),
    decimals = 1,
    use_seps = TRUE
  ) |>
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#f0f4ff"),
      gt::cell_text(weight = "bold")
    ),
    locations = gt::cells_column_labels()
  )
####
data.frame(
  IMR = c(fit.imrS$lambda_m, fit.imrS$lambda_m),
  IMR_matern = c(fitsim$rows$best_parameter, fitsim$cols$best_parameter),
  IMR_shuffled = c(fitsimF$rows$best_parameter, fitsimF$cols$best_parameter),
  IMR_identity = c(fitsimI$rows$best_parameter, fitsimI$cols$best_parameter),
  SoftImpute = c(fitsi$lambda, fitsi$lambda),
  row.names = c("λ_row", "λ_column")
) |>
  round(2) |>
  tibble::rownames_to_column("parameter") %>%
  dplyr::select(parameter, IMR_matern, IMR_identity, IMR, SoftImpute, IMR_shuffled) %>%
  gt::gt() |>
  gt::tab_header(
    title = "Frobeinus regularization parameter value",
    subtitle = "Both parameters are equal for the original model and SoftImpute"
  ) |>
  gt::cols_label(
    parameter    = "Parameter",
    IMR_matern   = "(Matern)",
    IMR_identity = "(identity)",
    IMR          = "(original)",
    SoftImpute   = "SoftImpute",
    IMR_shuffled = "(shuffled)"
  ) |>
  gt::fmt_number(
    columns = c(IMR_matern, IMR_identity, IMR, SoftImpute, IMR_shuffled),
    decimals = 1,
    use_seps = T
  ) |>
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#f0f4ff"),
      gt::cell_text(weight = "bold")
    ),
    locations = gt::cells_column_labels()
  )

res %>% gt::gt()
