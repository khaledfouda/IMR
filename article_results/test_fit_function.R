library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document();
load_all()
#devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

n = 800; m = 900;
dat <-
generate_simulated_data(n, m, 3, 5, 5, 0.9,sparsity_beta = 0.5, sparsity_gamma = 0.5,
                        prepare_for_fitting = F,mv_coeffs = T,seed = 2025)
colnames(dat$X) <- paste0("R",  1:ncol(dat$X))
colnames(dat$Z) <- paste0("C",  1:ncol(dat$Z))
#======
# d1 <-  matrix(rbinom(n*n,n,.2), n, n);d1 <- (d1 + t(d1)) / 2
d2 <-  matrix(rbinom(m*m,m,.2), m, m);d2 <- (d2 + t(d2)) / 2
# S <- generate_similarity(d1, invert = T, jitter = 1);S
S2 <- imr_similarity(d2, invert = T, jitter = 1);S2
data <- IMR:::imr_data(dat$Y, dat$X, dat$Z,val_prop = 0.2,
                          similarity_rows = NULL, similarity_cols = S2);data
data <- update(data, col_similarity = FALSE, row_similarity = FALSE,
                shared_beta = FALSE, intercept_row = T, intercept_col = T);data

grid <- IMR::imr_tune_grid(laplace = c(0,NA,40,3),
                           rank = c(2,15, 1));grid
convergence <- IMR::imr_convergence(600, 1e-5, FALSE,ls_initial = T); convergence

# grid$beta$max <- 1.695
# grid$gamma$max <- 1.815
# grid$laplace$max <- 92
grid <- imr_set_grid_limits(data, grid, convergence=convergence, verbose=2); grid

cv_out <- imr_tune(data, grid, warm_start = NULL, n_cores = 9,
                           final_fit = TRUE, convergence=convergence,
                     verbose=1, seed=2025)
fit <- cv_out$fit
fit
summary(fit)
data

fit <- imr_fit(data=data, rank = 3, lambda_m =.002,
                    lambda_beta = 0, lambda_gamma = 0,
                     warm_start = NULL,
                    convergence = convergence);

data$model$row_covariates <- FALSE

fit <- softImpute::softImpute(data$Y, cv_out$fit$meta$rank, cv_out$fit$meta$lambdas["M"])
fit <- structure(list(coefficients=fit), class = "imr_fit")

coefs <- coef(fit)
rec <- IMR::reconstruct(fit, data)
dat$test <- (dat$theta * (1 - dat$mask)) %>% IMR::as.Incomplete()
recp <- IMR::reconstruct_partial(fit, data, dat$test, TRUE)
IMR::evaluate(recp@x, dat$test@x, "all") %>% as_tibble()

beta.estim <- IMR:::inv(data$Xr) %*% fit$coefficients$beta
gamma.estim <- tcrossprod(fit$coefficients$gamma, IMR:::inv(data$Zr))
IMR::evaluate(beta.estim, dat$beta)
IMR::evaluate(gamma.estim, dat$gamma)
IMR::evaluate(coefs$u %*% (t(coefs$v) * (coefs$d)), dat$M )

mean(dat$beta == 0 & beta.estim == 0)

beta.estim[dat$beta==0]

cv_out$history -> history
history %>% head(20)
#--------------------------------------------------------------------------

# parameters for lambda_lasso_max
target = "beta"
rank = 2
lambda_m = lambda_beta = lambda_gamma= 0.1
intercept_row = FALSE
intercept_col = FALSE
shared_beta = shared_gamma = shared_effects = FALSE
bisection_iter = 15
verbose = 3
n_cores = 9
seed = 2021
iter = 1
fixed_other_lasso = 0
warm_start = NULL
error_function = IMR::error_metrics$rmse
default_lambda_beta = default_lambda_gamma = 0

rank = 2
lambda_m = .2
target="beta"
verbose = 2
bisection_iter = 5
i=1

#====
fit
summary(fit)

x <- coef(fit)

fit <- imr_fit(data,
               rank = 0,
               lambda_m = 2,
               lambda_beta = lambda_beta,
               lambda_gamma = lambda_gamma,
               convergence = convergence
)

imr_get_lambda_lasso_max(data, "gamma",
                         rank = 2,
                         lambda_m = .1,
                         convergence = convergence,
                         bisection_iter = bisection_iter,
                         verbose = verbose
)

baseline_fit <- imr_fit(
  data = data,
  rank = rank,
  lambda_m = lambda_m,
  lambda_beta = 0,
  lambda_gamma = 0,
  convergence = convergence
)


plot_tuning_trace(history, "path", "last", 9)
plot_tuning_trace(history, "convergence", "all", 13)
#' Plot IMR tuning diagnostics
#' @export
#' Plot IMR tuning diagnostics
#' @export
plot_tuning_trace <- function(history,
                              type = c("path", "convergence"),
                              show_iters = "last",
                              top_n_paths = 1) { # <--- NEW ARGUMENT

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to plot.")
  }

  type <- match.arg(type)

  # =========================================================================
  # TYPE 1: CONVERGENCE (The Macro View)
  # =========================================================================
  if (type == "convergence") {

    # Find the best Validation Error for every iteration overall
    best_list <- lapply(split(history, history$iter), function(df) {
      df[which.min(df$verror), ]
    })
    best_points <- do.call(rbind, best_list)
    best_points <- best_points[order(best_points$iter), ]

    p <- ggplot2::ggplot(best_points, ggplot2::aes(x = iter, y = verror)) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(color = "firebrick", size = 3) +
      ggplot2::scale_x_continuous(breaks = unique(best_points$iter)) +
      ggplot2::labs(
        title = "Alternating Optimization Convergence",
        subtitle = "Tracking the minimum validation error across iterations",
        x = "Tuning Iteration",
        y = "Best Validation Error"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())

    return(p)
  }

  # =========================================================================
  # TYPE 2: THE LAPLACE PATH (The Micro View)
  # =========================================================================
  if (type == "path") {

    max_iter <- max(history$iter)
    min_iter <- min(history$iter)

    # 1. Filter by requested iterations
    if (is.numeric(show_iters)) {
      history <- history[history$iter %in% show_iters, ]
    } else if (show_iters == "last") {
      history <- history[history$iter == max_iter, ]
    } else if (show_iters == "first_last") {
      history <- history[history$iter %in% c(min_iter, max_iter), ]
    } else if (show_iters == "all") {
      # Proceed without filtering
    } else {
      stop("show_iters must be 'last', 'first_last', 'all', or a numeric vector.")
    }

    if (nrow(history) == 0) stop("No data found for the specified iterations.")

    # 2. Filter down to the Top N paths per iteration to avoid 20+ facets
    # Create a unique ID for each [iter, beta, gamma] combination
    history$path_id <- paste(history$iter, history$lambda_beta, history$lambda_gamma, sep = "_")

    # Find the minimum error for each unique path
    path_summary <- aggregate(verror ~ path_id + iter, data = history, FUN = min)
    path_summary <- path_summary[order(path_summary$iter, path_summary$verror), ]

    # Keep only the top_n_paths per iteration
    top_paths <- do.call(rbind, lapply(split(path_summary, path_summary$iter), function(df) {
      head(df, top_n_paths)
    }))

    # Subset the history to only those top paths
    history <- history[history$path_id %in% top_paths$path_id, ]

    # 3. Format a clean label for the facets
    history$facet_label <- sprintf(
      "Iter %.1f\n\u03B2: %.4f | \u03B3: %.4f",
      history$iter, history$lambda_beta, history$lambda_gamma
    )

    # Force the facet levels to strictly follow the iteration order, then by performance
    history$facet_label <- factor(history$facet_label, levels = unique(history$facet_label))

    # 4. Find best points for the dashed line
    best_list <- lapply(split(history, history$facet_label), function(df) df[which.min(df$verror), ])
    best_points <- do.call(rbind, best_list)

    p <- ggplot2::ggplot(history, ggplot2::aes(x = lambda_laplace, y = verror)) +
      ggplot2::geom_line(color = "gray60", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(fill = rank_out), shape = 21, color = "white", size = 3) +
      ggplot2::geom_vline(data = best_points, ggplot2::aes(xintercept = lambda_laplace),
                          linetype = "dashed", color = "firebrick", alpha = 0.7) +
      ggplot2::geom_point(data = best_points, shape = 1, size = 5, color = "firebrick", stroke = 1.2) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_fill_viridis_c(name = "Active Rank (M)", option = "plasma") +
      ggplot2::facet_wrap(~ facet_label, scales = "free_y", ncol = min(3, length(unique(history$facet_label)))) +
      ggplot2::labs(
        title = "Laplace Tuning Path",
        subtitle = sprintf("Showing the top %d best-performing path(s) per iteration", top_n_paths),
        x = "Laplace Penalty (\u03BB_M) [Log Scale]",
        y = "Validation Error"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "gray95"),
                     legend.position = "bottom")

    return(p)
  }
}

cv_out$history
ploti(cv_out, "trace")
ploti(cv_out, "profile")

type = "profile"
x <- cv_out

ploti <- function(x, type = c("profile", "trace"), ...) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required to plot the tuning trace.")
    }

    type <- match.arg(type)

    # =========================================================================
    # COMMON SETUP
    # =========================================================================
    h <- x$history
    prm <- x$params

    if (is.null(h) || nrow(h) == 0 || is.null(prm) || nrow(prm) == 0) {
      stop("Missing tuning history or params in the model object.")
    }

    # Derive step and target for the History (used for Profile plot)
    h$tuning_target <- ifelse(h$iter %% 1 == 0, "Beta", "Gamma")
    h$step <- h$iter * 2 - 1

    # Derive step and target for the Params (used for Trace plot)
    prm$tuning_target <- ifelse(prm$iter %% 1 == 0, "Beta", "Gamma")
    prm$step <- prm$iter * 2 - 1

    # =========================================================================
    # TYPE 1: THE MICRO LOSS-SURFACE (Profile)
    # =========================================================================
    if (type == "profile") {

      # Create a unique ID for every single sweep
      h$sweep_id <- paste(h$step, h$tuning_target, h$lambda_beta, h$lambda_gamma, sep = "_")
      best_global <- h[which.min(h$verror), ]

      # Highlight only the BEST path from each step
      best_paths <- do.call(rbind, lapply(split(h, h$step), function(df) {
        best_sweep_id <- df$sweep_id[which.min(df$verror)]
        df[df$sweep_id == best_sweep_id, ]
      }))

      p <- ggplot2::ggplot(best_paths, ggplot2::aes(x = lambda_laplace, y = verror,
                                                    group = sweep_id,
                                                    color = step,
                                                    linetype = tuning_target)) +
        ggplot2::geom_line(linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_viridis_c(name = "Chronological Step", option = "plasma",
                                       breaks = unique(best_paths$step)) +
        ggplot2::scale_linetype_manual(name = "Tuning Phase", values = c("Beta" = "solid", "Gamma" = "dashed")) +
        ggplot2::geom_vline(xintercept = best_global$lambda_laplace,
                            color = "red", linetype = "dotted", alpha = 0.6) +
        ggplot2::geom_point(data = best_global, color = "red", size = 4, shape = 18) +
        ggplot2::labs(
          title = "Laplace Validation Profiles",
          subtitle = "Showing the optimal Laplace path for each alternating step",
          x = "Laplace Penalty (\u03BB_M) [Log Scale]",
          y = "Validation Error"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "right")

      return(p)
    }

    # =========================================================================
    # TYPE 2: THE MACRO CONVERGENCE (Trace)
    # =========================================================================
    if (type == "trace") {

      # No need to search the history anymore! Just use the pre-calculated params.
      best_df <- prm[order(prm$step), ]

      # Reshape into long format for facet_grid
      n_steps <- nrow(best_df)

      long_df <- data.frame(
        step = rep(best_df$step, 4),
        iter = rep(best_df$iter, 4),
        target = rep(best_df$tuning_target, 4),
        value = c(best_df$verror, best_df$lambda_beta, best_df$lambda_gamma, best_df$rank_out),
        metric = factor(rep(c("1. Val Error", "2. \u03BB_\u03B2", "3. \u03BB_\u03B3", "4. Latent Rank"),
                            each = n_steps))
      )

      # Custom x-axis labels (e.g., "1.0\n(Beta)", "1.5\n(Gamma)")
      x_labels <- sprintf("%.1f\n(%s)", best_df$iter, substr(best_df$tuning_target, 1, 1))

      p <- ggplot2::ggplot(long_df, ggplot2::aes(x = step, y = value)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 1) +
        ggplot2::geom_point(ggplot2::aes(color = target), size = 3) +
        ggplot2::scale_color_manual(name = "Active Tune", values = c("Beta" = "firebrick", "Gamma" = "darkorange")) +
        ggplot2::facet_grid(metric ~ ., scales = "free_y", switch = "y") +
        ggplot2::scale_x_continuous(breaks = best_df$step, labels = x_labels) +
        ggplot2::labs(
          title = "Alternating Optimization Convergence Trace",
          subtitle = "Tracking optimal hyperparameters across block coordinate descent",
          x = "Iteration & Tuning Target",
          y = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.placement = "outside",
          strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0),
          panel.grid.minor.x = ggplot2::element_blank(),
          legend.position = "bottom"
        )

      return(p)
    }
  }
