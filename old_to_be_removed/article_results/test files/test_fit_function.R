library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

# goal -> simulate fit function
n = 10
m = 8
r = 2
set.seed(2025)
dat <-
  generate_simulated_data(n, m, 2, 2, 5, 0.8,
                          sparsity_beta = .5, sparsity_gamma = 0.0,
                          structured_error_A = T,
                          structured_error_B = T,
                          prepare_for_fitting = T,mv_coeffs = T,seed = 2025)
inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)

hpar <- IMR::get_imr_default_hparams()
#hpar$laplacian_col <- IMR::decompose_symmetric_matrix(dat$similarity_cols,1e-6,grid[i])

Y <- inp.dat$Y

yn <- IMR:::naive_MC(Y) %>% IMR::svd_opt(k=r)

U <- yn$u
Dsq <- yn$d
V <- yn$v
lambda_m <- 1.2
xsvd <- base::eigen(dat$similarity_cols, symmetric=TRUE)
dc <- xsvd$values
Uc <- xsvd$vectors

xsvd <- base::eigen(dat$similarity_rows, symmetric=TRUE)
dr <- xsvd$values
Ur <- xsvd$vectors

#dr
#Ur
#-------------------------------------------------------------------------------
# test 1
VDsq <- t(t(V)*Dsq)
all.equal(VDsq, V %*% diag(Dsq))
# test 2
BD <- IMR:::update_B_cpp(Y, U, V, Dsq, lambda_m)
BD2 <- as.matrix(crossprod(Y, U) + VDsq)
BD3 <- BD2 %*% diag(1/(1+lambda_m/Dsq))
BD2 <- t(t(BD2)*(1/(1+lambda_m/Dsq)))
all.equal(BD, BD2)
all.equal(BD2, BD3)
# test 3 with lambda_c = 0
BD4 <- matrix(NA, m, r)
W <- as.matrix(crossprod(Y, U) + VDsq)
W <- crossprod(Uc, W)
for(i in 1:r){
  BD4[,i] <- ( Uc) %*% diag(Dsq[i]/(Dsq[i]+lambda_m+dc)) %*% W[,i]
}
all.equal(BD, BD4)
BD5 <- IMR:::update_B_sim_cpp(Y, U, V, Dsq, lambda_m, Uc, dc)
all.equal(BD, BD5)
all.equal(BD4,BD5)
BD[1:5,1:4]
BD4[1:5,1:4]
#----------------------------
# we now test on A
# test 1
UDsq <- t(t(U)*Dsq)
all.equal(UDsq, U %*% diag(Dsq))
# test 2
AD <- IMR:::update_A_cpp(Y, U, V, Dsq, lambda_m)
AD2 <- as.matrix(Y %*% V + UDsq)
AD3 <- AD2 %*% diag(1/(1+lambda_m/Dsq))
AD2 <- t(t(AD2)* (1/(1+lambda_m/Dsq)))
all.equal(AD, AD2)
all.equal(AD2, AD3)
# test 3 with lambda_r = 0
AD4 <- matrix(NA, n, r)
W <- as.matrix(Y %*% V + UDsq)
W <- crossprod(Ur, W)
for(i in 1:r){
  AD4[,i] <- Ur %*% diag(Dsq[i]/(Dsq[i]+lambda_m+dr)) %*% W[,i]
}
all.equal(AD, AD4)
AD5 <- IMR:::update_A_sim_cpp(Y, U, V, Dsq, lambda_m, Ur, dr)
all.equal(AD5, AD4)
all.equal(AD, AD5)
#================================================
# test the svd function
M <- matrix(rnorm(1000*10), 1000, 10)
microbenchmark::microbenchmark(
  svd1 = IMR:::svd_small_nc_cpp(M),
  svd2 = base::svd(M),
  svd3 = IMR::svd_opt(M, 4),
  times = 1000

)
svd1 = IMR:::svd_small_nc_cpp(M)
svd2 = base::svd(M)
svd3 = IMR::svd_opt(M, 4)
undo <- function(x) x$u %*% t(x$v * x$d)
all.equal( undo(svd1), undo(svd2))
all.equal(svd1$d, svd2$d)
all.equal(svd1$u, svd2$u)
norm(crossprod(svd1$v) - diag(ncol(svd1$v)), "F")
#===================================================================
#== this parts test the cv function with laplace
dat <-
  generate_simulated_data(600, 700, 5, 0, 0, 0.8,
                          sparsity_beta = .5, sparsity_gamma = 0.0,
                          structured_error_A = T,
                          structured_error_B = T,
                          prepare_for_fitting = T,mv_coeffs = T,seed = seed)
inp.dat <- IMR::prepare_data(dat$Y, dat$X, dat$Z, dat$similarity_rows, dat$similarity_cols)
hpar = get_imr_default_hparams(data$similarity_rows, data$similarity_cols, 0, 0)
out <- IMR:::imr.cv_laplace(inp.dat$model, 0, 0, T, T, hpar)
data <- inp.dat$model
lambda_beta = 0
lambda_gamma = 0
intercept_row = T
intercept_col = T
error_function = error_metric$rmse
thresh = 1e-6
maxit = 300
trace = TRUE
old_fit = NULL
ls_initial = FALSE



library(ggplot2)

ggplot(results_rows$history, aes(x = parameter, y = error, color = factor(step_size))) +
  geom_line() +
  geom_point() +
  labs(
    x = expression(lambda[r]),
    y = "Validation error",
    color = "Step size"
  ) +
  theme_minimal()
#==============================================
Y = inp.dat$Y
X = inp.dat$Xq
Z = inp.dat$Zq
intercept_row = FALSE
intercept_col = FALSE
r = 2
lambda_m = 0.1
lambda_beta = 0.1
lambda_gamma = 0.2
Ur = NULL;dr = NULL;Uc = NULL;dc = NULL
maxit = 300;thresh = 1e-5;trace = TRUE
warm_start = NULL
ls_initial = TRUE





