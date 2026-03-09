
library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article_results/simulation/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")



#--------------------------------------------------------------------------
sim_res <- function(dat, fit, name = "",
                    shared_information = TRUE,
                    error_metric = IMR:::error_metric$rmse,
                    digits = 5) {
  # dat should include: Y, theta, beta, X, Xr, Z, Zr
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  refit <- IMR::reconstruct(fit, dat,shared_information = shared_information)
  test.obs <- dat$Y == 0

  out <- data.frame(model = name)

  out$Theta <- error_metric(refit$estimates, dat$theta)
  out$Test <- error_metric(refit$estimates[test.obs], dat$theta[test.obs])
  out$Beta <- error_metric(refit$beta, dat$beta)
  out$M <- error_metric(refit$M, dat$M)

  out$rank <- length(fit$d)#qr(refit$estimates)$rank
  # true_singular <- IMR::svd_opt(dat$theta, tol = 1e-6)$d
  # num_singular <- length(true_singular)
  # if (length(fit$d) > num_singular) {
  #   estim_singular <- fit$d[1:num_singular]
  # } else if (length(fit$d) < num_singular) {
  #   estim_singular <- c(fit$d, rep(0, num_singular - length(fit$d)))
  # } else {
  #   estim_singular <- fit$d
  # }
  out$lambda <- if(is.null(fit$lambda_laplace)) NA else fit$lambda_laplace
  out <- out %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
  return(as_tibble(out))
}
#----------------------------------------------------------------------
get_metrics <- function(results, models, names, shared_information=TRUE, DATA=FALSE){
  stopifnot(length(models) == length(names))
  map_dfr(1:length(models), function(j) {
    mod <- models[j]
    reps <- results[[mod]]

    map_dfr(seq_along(reps), function(i) {
      res <- reps[[i]]
      tibble(
       # model = mod,
        rep   = i,
        `time in seconds`  = res$time #/ 60
      ) |>
        bind_cols(sim_res(if(DATA) results$DATA$dat else results$data[[i]],
                          if(mod=="LS") res else res$best_fit, names[j], shared_information))
    })
  })
}

require(gt)
metrics_to_table <- function(metrics,
                             title = "Model performance across replications",
                             order_metric = "Test",
                             decreasing = FALSE) {
  metrics_long <-
    metrics |>
    tidyr::pivot_longer(
      cols = -c(model, rep),
      names_to = "metric",
      values_to = "value"
    )

  metrics_sum <-
    metrics_long |>
    dplyr::group_by(model, metric) |>
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      .groups = "drop"
    )

  # Order models by the mean of the chosen metric (e.g., "Test")
  order_levels <-
    metrics_sum |>
    dplyr::filter(metric == order_metric) |>
    dplyr::arrange(if (decreasing) dplyr::desc(mean) else mean) |>
    dplyr::pull(model)

  # Put models missing `order_metric` at the end
  order_levels <- c(order_levels, setdiff(unique(metrics_sum$model), order_levels))

  metrics_sum |>
    dplyr::mutate(mean_sd = sprintf("%.3f (%.4f)", mean, sd)) |>
    dplyr::select(model, metric, mean_sd) |>
    tidyr::pivot_wider(names_from = metric, values_from = mean_sd) |>
    dplyr::mutate(model = factor(model, levels = order_levels)) |>
    dplyr::arrange(model) |>
    dplyr::mutate(model = as.character(model)) |>
    gt::gt(rowname_col = "model") |>
    gt::tab_header(title = title) |>
    gt::opt_table_outline() |>
    gt::tab_options(table.font.size = "small")
}
#----------------------------------------------------------------------
results0 <- readRDS("article_results/simulation/data/december_2025_results0.rds")
results1 <- readRDS("article_results/simulation/data/december_2025_results1.rds")
results2 <- readRDS("article_results/simulation/data/december_2025_results2.rds")
results3 <- readRDS("article_results/simulation/data/december_2025_results3.rds")
results4 <- readRDS("article_results/simulation/data/december_2025_results4.rds")

#-------------------------------------------------------------------
models <- c("IMRXZR", "IMRXZLS", "LS")
names <- c("[Covariate] with 0 initials", "[Covariate] with least-squares initial",
           "Least-squares initial only")
metrics0 <- get_metrics(results0, models, names, F, T)
metrics1 <- get_metrics(results1, models, names, T, T)

models <- c("IMR", "IMRXZ")
names <- c("[Intercept] ", "[Covariate][Intercept]")
metrics2 <- get_metrics(results2, models, names)
metrics21 <- get_metrics(results2_1, models, names)
metrics2 <- rbind(metrics2, metrics21)

models <- c("IMR", "IMRL", "IMRXZ", "IMRXZL")
names <- c("[None]", "[Laplace]", "[Covariate]", "[Covariate]  [Laplace]")
metrics3 <- get_metrics(results3, models, names)
metrics4 <- get_metrics(results4, models, names)
#-----------------------------------------------------------------
metrics_to_table(metrics0, title = html(paste0("mean(sd); ",
        "&beta; &isin; &Ropf;<sup>2&times800</sup>, ",
                                               "&Gamma; &isin; &Ropf;<sup>800&times2</sup>")))
metrics_to_table(metrics1, title = html(paste0("mean(sd); ",
                                             "&beta;, &Gamma; &isin; &Ropf;<sup>2</sup>")))
metrics_to_table(metrics2, title = "[covariates] [no correlation] [shared information]")
metrics_to_table(metrics3, title = "[covariates] [Matern on rows] [identity on columns] [shared informaiton]")
metrics_to_table(metrics4, title = "[covariates] [Matern on rows] [AR1 on columns] [shared information]")
#---------------------------------------------------------------------------------------

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)


# ---- True values ----
beta_true  <- as.numeric(results1$DATA$dat$fit_data$Rbeta[, 1])
gamma_true <- as.numeric(results1$DATA$dat$fit_data$gammaRt[1, ])

true_df <- tibble(
  term  = c("x1", "x2", "z1", "z2"),
  truth = c(beta_true, gamma_true)
)
get_hat <- function(res, what = c("beta", "gamma")) {
  what <- match.arg(what)

  # Preferred location
  if (!is.null(res$best_fit) && !is.null(res$best_fit[[what]])) {
    return(as.numeric(res$best_fit[[what]]))
  }

  # Common fallback for simple model objects (e.g., LS)
  if (!is.null(res[[what]])) {
    return(as.numeric(res[[what]]))
  }

  stop("Could not find ", what, " in result object.")
}

# ---- Build long data: 4 coefficients per replication ----
coef_long <- map_dfr(models, function(mod) {
  reps <- results1[[mod]]

  map_dfr(seq_along(reps), function(i) {
    res <- reps[[i]]

    beta_hat  <- get_hat(res, "beta")
    gamma_hat <- get_hat(res, "gamma")

    tibble(
      model    = mod,
      rep      = i,
      term     = c("x1", "x2", "g1", "g2"),
      estimate = c(beta_hat, gamma_hat)
    )
  })
})
models <- c(
  "IMRXZR"  = "fit with 0 initials",
  "IMRXZLS" = "fit with least-squares initials",
  "LS"      = "Least-squares initials only"
)

coef_long <- coef_long |>
  dplyr::mutate(model = factor(model, levels = names(models)))

ggplot(coef_long, aes(x = model, y = estimate, fill = model, color = model)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.20) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.6) +
  geom_hline(
    data = true_df,
    aes(yintercept = truth),
    color = "black",
    linewidth = 0.6,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ term, scales = "free_y", nrow = 1, ncol = 4) +
  scale_fill_discrete(labels = unname(models)) +
  scale_color_discrete(labels = unname(models)) +
  guides(
    color = "none",
    fill  = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_x_discrete(labels = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  )
