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
colnames(dat$X) <- paste0("R",  1:ncol(dat$X))
colnames(dat$Z) <- paste0("C",  1:ncol(dat$Z))
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
                    lambda_beta = 0, lambda_gamma = 0,
                     warm_start = NULL,
                    convergence = convergence);
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
colnames(data$Zq)


#' @export
coef.imr <- function(object, ...) {
  # Wrap the coefficients, active model, and variable names in a custom class
  structure(list(
    coefs = object$coefficients,
    active = object$data$model_model,
    x_names = object$data$meta$x_names, # Assuming you store col names
    z_names = object$data$meta$z_names
  ), class = "coef.imr")
}

#' @export
#'
coeffsprint(fit)

x = fit

colnames(data$Xq)
object=fit
summarise(fit)
coef.imr(fit)

coef.imr <- function(object, ...) {
  cat("\n=====================================================")
  cat("\n=== Summary of Incomplete Matrix Regression (IMR) ===")
  cat("\n=====================================================\n")

  #-- first part: training object metrics:

  sst <- object$meta$sum_squares$sst
  sse <- object$meta$sum_squares$sse
  ssrc <- object$meta$sum_squares$ssrc
  sscc <- object$meta$sum_squares$sscc
  ssri <- object$meta$sum_squares$ssri
  ssci <- object$meta$sum_squares$ssci
  ssm <- object$meta$sum_squares$ssm
  sss <- ssrc + sscc + ssri + ssci + ssm



  exp_var <- 1-(sse/sst)
  n <- length(object$residuals@x)
  mse <- sse / n
  rmse <- sqrt(mse)
  prepare_fmt <- function(x) ifelse(round(x,1) >= 0.1, sprintf("%2.1f%%",x), "< 0.1%")
  prepare_d <- function(x) ifelse(round(x,4) >= 1e-4, sprintf("%.4f",x), "< 1e-4")
  cat("\nGoodness of Fit (on observed entries)")
  cat("\n-------------------------------------------------------------\n")
  cat(sprintf("RMSE                                   : %s\n", prepare_d(rmse)))
  cat(sprintf("MSE                                    : %s\n", prepare_d(mse)))
  cat(sprintf("Total Explained Variance (1 - SSE/SST) : %s\n", prepare_fmt(100*exp_var)))
  cat(sprintf("Unexplained Variance (SSE/SST)         : %s\n", prepare_fmt(100*(sse/sst))))


  cat("\nRelative Contribution to Explained Variance")
  cat("\n-------------------------------------------------------------\n")
    cat(sprintf("Additive Components        : %2.1f%% of explained variance\n",100* sss/(sst-sse)))
  if(object$model$low_rank_component)
    cat(sprintf("    |-- Latent Matrix      :    %s\n", prepare_fmt(100*ssm/sss)))
  if(object$model$row_covariates)
    cat(sprintf("    |-- Row Covariates     :    %s\n", prepare_fmt(100*ssrc/sss)))
  if(object$model$col_covariates)
    cat(sprintf("    |-- Column Covariates  :    %s\n", prepare_fmt(100*sscc/sss)))
  if(object$model$intercept_row)
    cat(sprintf("    |-- Row Intercepts     :    %s\n", prepare_fmt(100*ssri/sss)))
  if(object$model$intercept_col)
    cat(sprintf("    |-- Column Intercepts  :    %s\n", prepare_fmt(100*ssci/sss)))

    cat(sprintf("\nRegularization & Overlap   : %s of explained variance\n", prepare_fmt(100*(1-(sss/(sst-sse))))))
  cat("(Measures the deviation from orthogonal additivity due to shrinkage penalties and missing data.)\n")
  cat("-------------------------------------------------------------\n")

  cat("\n Model Estimates")
  cat("\n-------------------------------------------------------------\n")

  summarize_covariates <- function(coefs, names, rows=FALSE, len=NA){

    shared <- nrow(coefs) == 1

    cat(sprintf("\n-- %s Covariates --\nMode: %s (%s)\n",
            ifelse(rows, "Row", "Column"),
            ifelse(shared, paste0("Shared across ",ifelse(rows, "columns", "rows")),
                   paste0(ifelse(rows, "Column", "Row"),"-specific")),
            ifelse(shared, sprintf("%s = %d", ifelse(rows,"p","q"), ncol(coefs)),
                   sprintf("%d x %d matrix", nrow(coefs), ncol(coefs)))))
    if(shared){
      vec <- data.frame(Estimate = t(coefs))
      rownames(vec) <- names
      vec <- round(vec, 4)
      print(format(vec, justify = "right"))
    }else{
      cat(sprintf("Summary of effects across %d %s\n", len,
                  ifelse(rows, "columns", "rows")))
      summary <- data.frame(
        Mean = colMeans(coefs),
        SD  = apply(coefs, 2, sd),
        Min = apply(coefs, 2, min),
        Max = apply(coefs, 2, max),
        Sparsity = sprintf("%.2f%%",colMeans(round(coefs,6)==0) * 100),
        `L2-Norm` = apply(coefs, 2, function(x)sqrt(sum(x^2))),
        row.names = names
      )
      summary <- summary[order(summary$L2.Norm, decreasing = TRUE),]
      summary[,-5] <- round(summary[,-5], 4)
      print(format(summary, justify = "right"))

    }


  }

  if(object$model$row_covariates)
    summarize_covariates(t(object$coefficients$beta), object$meta_data$names_X_vars,
                         TRUE,
                         object$meta_data$dimensions[2])
  if(object$model$col_covariates)
  summarize_covariates(object$coefficients$gamma, object$meta_data$names_Z_vars,
                       FALSE,
                       object$meta_data$dimensions[1])

  if(object$model$intercept_row || object$model$intercept_col){
    cat("\n-- Intercepts --")
    if(object$model$intercept_row){
      coefs <- object$coefficients$beta0
      cat(sprintf("\nRow Intercepts    (n=%d) | ", length(coefs)))
      cat(sprintf("Mean: %.4f | SD: %.4f | Min: %.4f | Max: %.4f\n",
                  mean(coefs), sd(coefs), min(coefs), max(coefs)))

    }
    if(object$model$intercept_col){
      coefs <- object$coefficients$gamma0
      cat(sprintf("Column Intercepts (m=%d) | ", length(coefs)))
      cat(sprintf("Mean: %.4f | SD: %.4f | Min: %.4f | Max: %.4f\n",
                  mean(coefs), sd(coefs), min(coefs), max(coefs)))

    }
  }
  if(object$model$low_rank_component){
    cat("\n-- Latent Component --")
    cat(sprintf("\nRank (r): %d\nSingular Values:\n", length(object$coefficients$d) ))
    print(object$coefficients$d)
  }

  invisible()
}


x <- fit
print.imr(fit)
