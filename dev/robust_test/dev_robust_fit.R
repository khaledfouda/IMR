require(magrittr)
require(tidyverse)
pkgbuild::clean_dll(); Rcpp::compileAttributes(); devtools::document()
devtools::load_all(recompile = TRUE)
devtools::test() # check

set.seed(1)
Y <- matrix(rnorm(60 * 30), 60, 30)
outlier_idx <- sample(length(Y), 60)
Y[outlier_idx] <- rcauchy(60, scale = 30)
A <- Y
test_idx <- sample(length(Y), 300)
Y[test_idx] <- NA
d <- update(imr_data(Y = Y, val_prop = 0), row_intercept = TRUE)

f0 <- imr_fit(d, rank = 3, convergence = imr_convergence(maxit = 30, trace=TRUE))
f1 <- imr_fit(d, rank = 3, huber_shift = 5,
              convergence = imr_convergence(maxit = 600, thresh=1e-6, trace = TRUE,ls_initial = FALSE))

print(f1); summary(f1)
print(f0); summary(f0)

o0 = reconstruct(f0, d)$estimates
o1 = reconstruct(f1, d)$estimates
#------------------------------------------------------------------------------------------------
idx1 <- intersect(test_idx, outlier_idx)
idx2 <- setdiff(test_idx, outlier_idx)
list(
  test      = test_idx,
  outliers  = idx1,
  nonoutliers = idx2
) |>
  imap(\(idx, subset_name) {
    bind_rows(
      evaluate(o0[idx], A[idx]) |> mutate(method = "Frob"),
      evaluate(o1[idx], A[idx]) |> mutate(method = "Huber")
    ) |>
      mutate(subset = subset_name)
  }) |>
  bind_rows()
#=========================================================================

set.seed(2026)
n <- 200; q <- 100; r_true <- 3

U0 <- matrix(rnorm(n * r_true), n, r_true)
V0 <- matrix(rnorm(q * r_true), q, r_true)
M  <- U0 %*% t(V0) / sqrt(r_true)
Y_clean <- M + matrix(rnorm(n * q), n, q)
n_out   <- round(0.03 * n * q)
out_idx <- sample(n * q, n_out)
Y_obs   <- Y_clean
Y_obs[out_idx] <- Y_obs[out_idx] +
  sample(c(-1, 1), n_out, TRUE) * runif(n_out, 10, 30)

test_idx <- sample(n * q, round(0.2 * n * q))
Y_in <- Y_obs; Y_in[test_idx] <- NA

dat <- imr_data(Y = Y_in, val_prop = 0)
cv  <- imr_convergence(maxit = 600, thresh = 1e-5)

f_ls <- imr_fit(dat, rank = 3, lambda_m = 0, convergence = cv)
f_hb <- imr_fit(dat, rank = 3, lambda_m = 0, huber_shift = 1.345, convergence = cv)

o_ls <- reconstruct(f_ls, dat, trace = FALSE)$estimates
o_hb <- reconstruct(f_hb, dat, trace = FALSE)$estimates

clean <- setdiff(test_idx, out_idx)

cmp <- function(idx, label) rbind(
  data.frame(evaluate(rep(mean(Y_in, na.rm = TRUE), length(idx)), M[idx]),
             method = "Constant"),
  data.frame(evaluate(o_ls[idx], M[idx]), method = "Frob"),
  data.frame(evaluate(o_hb[idx], M[idx]), method = "Huber")
) |> transform(subset = label)

rbind(cmp(clean, "clean"), cmp(intersect(test_idx, out_idx), "contaminated"))
cat(sprintf("\nHuber: c = %.4f, clipped = %.1f%%, iters = %d (LS: %d)\n",
            f_hb$meta$huber$huber_c, 100 * f_hb$meta$huber$prop_clipped,
            f_hb$meta$n_iter, f_ls$meta$n_iter))


