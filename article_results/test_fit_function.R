library(devtools)
clean_dll(); Rcpp::compileAttributes(); document();
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

n = 800; m = 900;
dat <-
generate_simulated_data(n, m, 3, 5, 5, 0.7,sparsity_beta = .50, sparsity_gamma = .50,
                        prepare_for_fitting = F,mv_coeffs = T,seed = 2025)

#======
# d1 <-  matrix(rbinom(n*n,n,.2), n, n);d1 <- (d1 + t(d1)) / 2
d2 <-  matrix(rbinom(m*m,m,.2), m, m);d2 <- (d2 + t(d2)) / 2
# S <- generate_similarity(d1, invert = T, jitter = 1);S
S2 <- imr_similarity(d2, invert = T, jitter = 1);S2
data <- IMR:::imr_data(dat$Y, dat$X, dat$Z,val_prop = 0.2,
                          similarity_rows = NULL, similarity_cols = S2);data
data <- update(data,with_col_similarity = FALSE, with_row_similarity = FALSE,
                shared_beta = TRUE, intercept_row = T);data

grid <- IMR::imr_tune_grid(laplace = c(0,NA,40,2));grid
convergence <- IMR::imr_convergence(600, 1e-5, FALSE,ls_initial = T); convergence

grid$beta$max <- 1.695
grid$gamma$max <- 1.815
grid$laplace$max <- 92
grid <- imr_set_grid_limits(data, grid, convergence=convergence, verbose=1); grid

cv_out_g <- imr_tune(data, grid, warm_start = NULL, n_cores = 9,
                           final_fit = TRUE, convergence=convergence,
                     verbose=1, seed=2025)


fit <- IMR::imr_fit(data, rank = 5, lambda_m = 5,
                    lambda_beta = 0.2, lambda_gamma = 0.2,
                     warm_start = NULL,
                    convergence = convergence);fit
fit <- cv_out_g$fit
rec <- IMR::reconstruct(fit, data)
dat$test <- (dat$theta * (1 - dat$mask)) %>% IMR::as.Incomplete()
recp <- IMR::reconstruct_partial(fit, data, dat$test, TRUE)
IMR::evaluate(recp@x, dat$test@x, "all") %>% as_tibble()
#--------------------------------------------------------------------------

# parameters for lambda_lasso_max
target = "beta"
rank = 2
lambda_m = lambda_beta = lambda_gamma= 0.1
intercept_row = FALSE
intercept_col = FALSE
shared_beta = shared_gamma = shared_effects = FALSE
bisection_iter = 15
verbose = 2
n_cores = 9
seed = 2021
iter = 1
fixed_other_lasso = 0
warm_start = NULL
error_function = IMR::error_metrics$rmse
default_lambda_beta = default_lambda_gamma = 0
#====
fit
