library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")

# goal -> simulate fit function
n = 400
m = 300
r = 4
set.seed(2025)
dat <-
  generate_simulated_data(n, m, 10, 20, 0, 0.8,
                          sparsity_beta = .5, sparsity_gamma = 0.0,
                          structured_error_A = T,
                          structured_error_B = T,
                          prepare_for_fitting = T,mv_coeffs = T,seed = 2025)
inp.dat <- list(
  Y = dat$fit_data$Y_full,
  y_train = dat$fit_data$train,
  y_valid = dat$fit_data$valid,
  Xq = dat$fit_data$X$Q,
  Zq = dat$fit_data$Z$Q
)
hpar <- IMR::get_imr_default_hparams()
#hpar$laplacian_col <- IMR::decompose_symmetric_matrix(dat$similarity_cols,1e-6,grid[i])

Y <- inp.dat$Y

yn <- IMR:::naive_MC(Y) %>% IMR::opt_svd(k=r)

U <- yn$u
Dsq <- yn$d
V <- yn$v
lambda_M <- 1.2
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
BD <- IMR:::update_B_cpp(Y, U, V, Dsq, lambda_M)
BD2 <- as.matrix(crossprod(Y, U) + VDsq)
BD3 <- BD2 %*% diag(1/(1+lambda_M/Dsq))
BD2 <- t(t(BD2)*(1/(1+lambda_M/Dsq)))
all.equal(BD, BD2)
all.equal(BD2, BD3)
# test 3 with lambda_c = 0
BD4 <- matrix(NA, m, r)
W <- as.matrix(crossprod(Y, U) + VDsq)
W <- crossprod(Uc, W)
for(i in 1:r){
  BD4[,i] <- ( Uc) %*% diag(Dsq[i]/(Dsq[i]+lambda_M+dc)) %*% W[,i]
}
all.equal(BD, BD4)
BD5 <- IMR:::update_B_sim_cpp(Y, U, V, Dsq, lambda_M, Uc, dc)
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
AD <- IMR:::update_A_cpp(Y, U, V, Dsq, lambda_M)
AD2 <- as.matrix(Y %*% V + UDsq)
AD3 <- AD2 %*% diag(1/(1+lambda_M/Dsq))
AD2 <- t(t(AD2)* (1/(1+lambda_M/Dsq)))
all.equal(AD, AD2)
all.equal(AD2, AD3)
# test 3 with lambda_r = 0
AD4 <- matrix(NA, n, r)
W <- as.matrix(Y %*% V + UDsq)
W <- crossprod(Ur, W)
for(i in 1:r){
  AD4[,i] <- Ur %*% diag(Dsq[i]/(Dsq[i]+lambda_M+dr)) %*% W[,i]
}
all.equal(AD, AD4)
AD5 <- IMR:::update_A_sim_cpp(Y, U, V, Dsq, lambda_M, Ur, dr)
all.equal(AD5, AD4)
all.equal(AD, AD5)

