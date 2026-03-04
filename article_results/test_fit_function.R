library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document();
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

n = 800; m = 900;
dat <-
generate_simulated_data(n, m, 3, 5, 5, 0.9,sparsity_beta = 0, sparsity_gamma = 0,
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
                shared_beta = FALSE, intercept_row = F, intercept_col = F);data

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

fit <- IMR::imr_fit(data=data, rank = 3, lambda_m =.002,
                    lambda_beta = 0, lambda_gamma = 0,
                     warm_start = NULL,
                    convergence = convergence);

data$model$row_covariates <- FALSE

fit <- softImpute::softImpute(data$Y, cv_out$fit$meta$rank, cv_out$fit$meta$lambdas["M"])
fit <- structure(list(coefficients=fit), class = "imr_fit")


rec <- IMR::reconstruct(fit, data)
dat$test <- (dat$theta * (1 - dat$mask)) %>% IMR::as.Incomplete()
recp <- IMR::reconstruct_partial(fit, data, dat$test, TRUE)
IMR::evaluate(recp@x, dat$test@x, "all") %>% as_tibble()

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
