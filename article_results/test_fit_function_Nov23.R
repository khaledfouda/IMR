library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
source("./article_results/simulation/generate_simu_dat.R")


dat <-
generate_simulated_data(587, 196, 3, 5, 5, 0.7,sparsity_beta = .50, sparsity_gamma = .50,
                        prepare_for_fitting = T,mv_coeffs = T,seed = 2025)

dim(dat$Y)

# generate laplace:
# make this a function: BKTR laplacian matrices
library(BKTR)
bkdat <- BixiData$new()
spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 3)
temporal_kernel = BKTR::KernelSE$new()
spatial_kernel$set_positions(bkdat$spatial_positions_df)
temporal_kernel$set_positions(bkdat$temporal_positions_df)
spatial_kernel$kernel_gen()  %>% as.matrix() -> spatial_kernel
temporal_kernel$kernel_gen() %>% as.matrix() -> temporal_kernel
#====================================================================
hpar <- IMR::get_imr_default_hparams(spatial_kernel, temporal_kernel, 1, 1)
#====================================================================
fit <- IMR::imr.fit(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,
                    #Uc = hpar$laplacian_col$U, dc = hpar$laplacian_col$d,
                    #Ur = hpar$laplacian_row$U, dr = hpar$laplacian_row$d,
                    r=6, lambda_M = 3.23, lambda_beta=.000, lambda_gamma=0,
                    trace=F, ls_initial = T,intercept_row = T, intercept_col = T)
quick_camc_simu_res(dat, fit)
#========================================================
fit2 <- IMR::imr.cv_M(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$X$Q,
                     # hpar=hpar,
                      dat$fit_data$Z$Q, dat$fit_data$Y_full, 0, 0, T, T,ls_initial = T)
quick_camc_simu_res(dat, fit2$fit)
#===========================================================

## ---------------------------------------------------------
## 1. Original (loop + diag) implementation
## ---------------------------------------------------------

method_loop <- function(Y, V, U, Dsq, Ur, dr, lambda_M) {
  # Y:  m x n
  # V:  n x r
  # U:  m x r
  # Dsq: length r
  # Ur: m x m
  # dr: length m
  # lambda_M: scalar
  #
  # Returns: A_mat (m x r), using your original structure

  m <- nrow(Y)
  r <- length(Dsq)

  # Your original code:
  partial <- Y %*% V + U %*% diag(Dsq)
  partial <- crossprod(Ur, partial %*% diag(Dsq))

  A_mat <- matrix(0, nrow = nrow(Ur), ncol = r)

  for (j in seq_len(r)) {
    A_mat[, j] <- Ur %*% diag(1 / (dr + Dsq[j] + lambda_M)) %*% partial[, j]
  }

  A_mat
}


## ---------------------------------------------------------
## 2. Optimized vectorized implementation (no diag, no loop)
## ---------------------------------------------------------

method_vec <- function(Y, V, U, Dsq, Ur, dr, lambda_M) {
  # Same arguments as method_loop, but vectorized.

  m <- nrow(Y)
  r <- length(Dsq)

  # 1) Build partial, avoiding diag()
  # Y %*% V  -> m x r
  partial <- Y %*% V

  # U %*% diag(Dsq) -> scale columns of U by Dsq
  partial <- partial + sweep(U, 2L, Dsq, `*`)   # m x r

  # partial %*% diag(Dsq) -> scale columns of partial by Dsq again
  partial <- crossprod(Ur, sweep(partial, 2L, Dsq, `*`))  # m x r

  # 2) Build scaling factors for each (row i, column j)
  # denom[i, j] = dr[i] + Dsq[j] + lambda_M
  denom <- outer(dr + lambda_M, Dsq, `+`)  # m x r
  coef  <- 1 / denom                       # m x r

  # 3) Apply scaling and final multiply
  tmp   <- partial * coef                  # m x r
  A_mat <- Ur %*% tmp                      # m x r

  A_mat
}


## ---------------------------------------------------------
## 3. Simulation / comparison function
## ---------------------------------------------------------

run_sim <- function(m = 500, n = 400, r = 10, lambda_M = 0.1, seed = 1L) {
  set.seed(seed)

  # Generate random matrices with compatible dimensions
  Y   <- matrix(rnorm(m * n), nrow = m, ncol = n)
  V   <- matrix(rnorm(n * r), nrow = n, ncol = r)
  U   <- matrix(rnorm(m * r), nrow = m, ncol = r)
  Dsq <- abs(rnorm(r)) + 0.1                 # keep positive
  dr  <- abs(rnorm(m)) + 0.1                 # positive diagonal-like

  # Make Ur approximately orthonormal (so it's well-conditioned)
  tmp <- matrix(rnorm(m * m), nrow = m, ncol = m)
  Ur  <- qr.Q(qr(tmp))                       # m x m, orthonormal columns

  # 1) Check numerical agreement
  res_loop <- method_loop(Y, V, U, Dsq, Ur, dr, lambda_M)
  res_vec  <- method_vec (Y, V, U, Dsq, Ur, dr, lambda_M)

  max_diff <- max(abs(res_loop - res_vec))
  all.equal(res_loop, res_vec)
  cat("Dimensions: m =", m, " n =", n, " r =", r, "\n")
  cat("Max |loop - vec| =", max_diff, "\n\n")

  # 2) Benchmark performance (if microbenchmark is available)
  if (requireNamespace("microbenchmark", quietly = TRUE)) {
    library(microbenchmark)

    bench <- microbenchmark(
      loop = method_loop(Y, V, U, Dsq, Ur, dr, lambda_M),
      vec  = method_vec (Y, V, U, Dsq, Ur, dr, lambda_M),
      times = 50L
    )

    print(bench)
  } else {
    cat("Package 'microbenchmark' not installed.\n",
        "Install it with install.packages('microbenchmark') to see timing comparison.\n")
  }
}


## ---------------------------------------------------------
## 4. Run a few scenarios
## ---------------------------------------------------------

# Medium size
run_sim(m = 300, n = 300, r = 10, lambda_M = 0.1, seed = 1L)

# Larger example (adjust to your machine’s RAM/CPU)
run_sim(m = 800, n = 600, r = 20, lambda_M = 0.1, seed = 2L)










#========================================================================
m = 10; r=5;
W1 <- crossprod(matrix(rnorm(m*m),m)) + diag(0.1, m)
eig <- IMR::decompose_symmetric_matrix(W1,1e-10)
W2 <- rexp(r) + 0.1
W3 <- matrix(rnorm(m*r), m, r)
d <- runif(r, 0.5, 1.5)
eig$U <- svd(W1)$u
eig$d <- svd(W1)$d

compute_B_naive <- function(W1, W2, W3) {
  m <- nrow(W1)
  r <- ncol(W3)
  B <- matrix(NA_real_, nrow = m, ncol = r)
  I_m <- diag(1, m, m)

  for (j in seq_len(r)) {
    W_j <- W1 + W2[j] * I_m
    # Solve W_j %*% B[, j] = W3[, j]
    B[, j] <- solve(W_j) %*% W3[, j]
  }

  B
}
compute_B_fast <- function(W3, Uc, dc, W2) {
  # Rotate W3 into eigenbasis of W1: t(Uc) %*% W3
  tmp <- crossprod(Uc, W3)  # m x r
  # denom[i, j] = dc[i] + w2_diag[j]
  denom <- outer(dc, W2, "+")
  tmp <- tmp / denom
  # Rotate back: B = Uc %*% tmp
  B <- Uc %*% tmp
  #-------
  B <- matrix(NA, m, r)
  for(j in seq_len(r)){
    parta <- 1/(dc+W2[j])
    B[,j] <- Uc %*% diag(parta) %*% t(Uc) %*% W3[,j]
  }
  B
}

B1 <- compute_B_naive(W1, W2, W3)
B2 <- compute_B_fast(W3, eig$U, eig$d, W2)
B1[1:5,1:5]
verify_sol(W1, W2, W3, B1)
verify_sol(W1, W2, W3, B2)
B = B2

verify_sol <- function(W1, W2, W3, B){
  all.equal(W3, B %*% diag(W2) + W1 %*% B, tolerance=1e-4)
}

#----------
b = 0.4; 40
W1 <- crossprod(matrix(rnorm(m*m),m)) + diag(0.1, m)
W1 <- temporal_kernel; m = nrow(W1)

I <- diag(1, nrow(W1))
eig1 <- decompose_symmetric_matrix(W1, 1e-6)
eig2 <- list(U = svd(W1)$u, d=svd(W1)$d)
length(eig2$d)
length(eig1$d)
all.equal(W1, eig1$U %*% diag(eig1$d) %*% t(eig1$U), tolerance=1e-3 )
eig1$d
m1 <- solve(W1 + b * I)
m2 <- eig2$U %*% diag(1/(eig2$d+b) ) %*% t(eig2$U)
all.equal(m1, m2)
m3 <- (eig1$U %*% diag(1/(eig1$d)) %*% t(eig1$U)  )
all.equal(round(m1,3), round(m3,3))
m4 <- (eig1$U %*% diag(1/(eig1$d+b) ) %*% t(eig1$U) )
all.equal(m1, m4)

# test the fit function
fit <- IMR::imr.fit(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,

                r=6, lambda_M = 3.23, lambda_beta=.000, lambda_gamma=0,
                trace=T, ls_initial = T,intercept_row = T, intercept_col = T)
quick_camc_simu_res(dat, fit)

fitm <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q, dat$fit_data$Z$Q,
                                 .0001, .0001, T, T, trace=T)

fit2 <- IMR::imr.cv_M(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$X$Q,
                      dat$fit_data$Z$Q, dat$fit_data$Y_full, 0, 0, T, T,ls_initial = T)
quick_camc_simu_res(dat, fit2$fit)


mfit <- IMR::imr.cv_M(
  y_train     = dat$fit_data$train,
  y_valid     = dat$fit_data$valid,
  X           = NULL,
  Z           = dat$fit_data$Z$Q,
  Y_full      = dat$fit_data$Y_full,
  lambda_beta = 0,
  lambda_gamma = 0,
  intercept_row = T,
  intercept_col = T,
  trace       = T
)


fit <- IMR::imr.fit(
  Y = dat$fit_data$train,
  X = NULL,
  Z = dat$fit_data$Z$Q,
  r = 7,
  lambda_M = 62.34,
  lambda_beta = 0,
  lambda_gamma = 0,
  intercept_row = T,
  intercept_col = T,
  #warm_start = old_fit,
  trace = T,
  thresh = 1e-6,
  maxit = 300
)



IMR::get_lambda_lasso_max(dat$fit_data$train, NULL,
                          dat$fit_data$Z$Q, dat$fit_data$valid,
                          dat$fit_data$Y_full, T, T, verbose = T
                          )


future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()
hpar$M$lambda_factor <- 1
hpar$M$n.lambda <- 20
hpar$M$rank.init <- 2
hpar$M$rank.step <- 2
hpar$M$rank.max  <- 15
hpar$M$lambda_max <- 80
hpar$beta$n.lambda <- 10
hpar$gamma$n.lambda <- 10
fit32 <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                    Z = dat$fit_data$Z$Q,intercept_row = F,
                    hpar = hpar, seed = 2025,
                    intercept_col = F, verbose=2)
quick_camc_simu_res(dat, fit32$fit)



fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                    trace=T,tol = 5, n.lambda=80, rank.init = 2,
                    rank.step = 1
                    #lambda0_fun = IMR::get_lambda_M_max
                    )
quick_camc_simu_res(dat, fitsi$fit)




fitmmci <- MCCI.cv(dat$Y, dat$X, dat$mask, numCores = 9)
quick_camc_simu_res(dat, fitmmci$fit, T)

#-------------------------------------------------------
#---- testing IMC no bad local minima 2022
library(reticulate)
reticulate::py_run_string("import sys; sys.path.insert(0, r'~/Research/IMR/other_models/')")
np <- reticulate::import("numpy")
sp    <- import("scipy.sparse")
reticulate::source_python("./other_models/GNIMC.py")
GNIMC <- reticulate::py$GNIMC

pyopt <- list(
  X = dat$Y,
  omega = dat$mask,#sp$csr_matrix(dat$mask, dtype = np$float64),
  A = dat$X,
  B = dat$Z,
  verbose = TRUE,
  init_option = py$INIT_WITH_SVD
)
out <- do.call(GNIMC, c(list(rank=10L), pyopt))
estimate <- out[[1]]

xnorm <- np$linalg$norm(dat$theta[dat$Y==0])
np$linalg$norm(estimate[dat$Y==0] - dat$theta[dat$Y==0])/xnorm

IMR::error_metric$rmse(estimate[dat$Y==0], dat$theta[dat$Y==0])

#---------------------------------------------------------
f <- IMR::imr.cv_M

worker <- function(args) {
  on.exit(if (!is.null(p)) p(), add = TRUE)
  do.call(f, c(args, list(y_train = y_train,
                          y_valid = y_valid,
                          X = X,
                          Z = Z,
                          Y_full = Y,
                          intercept_row = intercept_row,
                          intercept_col = intercept_col,
                          hpar = hpar,
                          error_function = error_function,
                          thresh = thresh,
                          maxit = maxit,
                          trace = inner_trace,
                          seed = seed)))
}

hpar$M$lambda_max = 30

IMR::imr.cv_M(
  lambda_beta = 0,
  lambda_gamma = 0,
  y_train = y_train,
  y_valid = y_valid,
  X = X,
  Z = Z,
  Y_full = Y,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  hpar = hpar,
  error_function = error_function,
  thresh = thresh,
  maxit = maxit,
  trace = T,
  seed = seed
)


IMR::imr.fit(
  lambda_beta = 0,
  lambda_gamma = 0,
  Y = y_train,
  r = 2,
  X = X,
  Z = Z,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  thresh = thresh,
  maxit = 300,
  trace = T,
  ls_initial = F
) -> o

dims <- dim(y_train)
nr <- dims[1]
nc <- dims[2]
init <- opt_svd(IMR:::naive_MC(y_train), 2, nr, nc, FALSE, FALSE)

mat <- as.matrix(y_train)


microbenchmark::microbenchmark(
a = partial_crossprod(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
b =  partial_crossprod_cpp(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
times = 500
)
