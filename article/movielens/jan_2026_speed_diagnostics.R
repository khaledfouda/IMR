# after runing main.R >>

# function paramter (cv_lambda)
data = fdata = dat;
intercept_row = TRUE
intercept_col = TRUE
hpar = hparam
error_function = IMR:::error_metric$rmse
thresh = 1e-3
maxit = 300
trace = 2
warm_start = NULL
ls_initial = TRUE
shared_information = FALSE
num_cores = 9
final_fit = TRUE
seed = 2025
final_thresh = 1e-6
final_maxit = 1000
init_thresh = 1e-4
init_maxit = 500
# extras
lambda_laplace = 15
r = 12
lambda_beta = 1.5
lambda_gamma = .3
#--------------------------
require(bench)
# testing effect of warm_start and which to use to speed convergence / improve performance
warm_start = fit_init
bench_time(
  {
    results = data.frame()
    for(lambda_laplace in lamseq[5:1]){

      res <- single_fit(
        lambda_laplace = lambda_laplace,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        data = data,
        intercept_row = intercept_row,
        intercept_col = intercept_col,
        hpar = hpar,
        shared_information = shared_information,
        error_function = error_function,
        thresh = thresh,
        trace = trace,
        maxit = maxit,
        ls_initial = ls_initial,
        seed = seed,
        return_fit = TRUE,
        warm_start = NULL
      )
      warm_start = res$fit
      results <- rbind(results, res$res)
    }
  }
)

results1
results5

results6 <- results
results1 <- results


results6 %>%
  arrange(error) %>% head(1)

# 1.  rank fit function speed (single fit)
bench::bench_time(rank_fit <- rank_fit_function(r, data, hpar, shared_information, lambda_laplace,
                  intercept_row, intercept_col, trace, thresh, maxit, TRUE)) ->
  single_fit_time; single_fit_time



# 2. one fit
require(bench)
bench_time(fit_init <- IMR::imr.fit(
  Y = fdata$y_train,
  X = fdata$Xq,
  Z = fdata$Zq,
  r = hpar$rank$default,
  lambda_M = lambda_laplace,
  lambda_beta = hpar$beta$value,
  lambda_gamma = hpar$gamma$value,
  intercept_row = intercept_row,
  intercept_col = intercept_col,
  shared_information = shared_information,
  Ur = hpar$laplacian_row$U,
  dr = hpar$laplacian_row$d,
  Uc = hpar$laplacian_col$U,
  dc = hpar$laplacian_col$d,
  warm_start = NULL,
  trace = T,
  thresh = 1e-3,
  maxit = maxit,
  ls_initial = F
))



bench_time(
  bmax <- IMR::get_lambda_lasso_max(
    y_train = data$y_train,
    X = data$Xq,
    # y_valid = data$y_valid,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    maxit = 50,
    thresh = 1e-3,
    init_maxit = 100,
    init_thresh = 1e-4,
    r = 2,
    verbose = 100
  )
);

bmax


for(lambda_laplace in seq(10,0,-1)){


bench_time(
  onfitf <- single_fit(
    lambda_laplace = lambda_laplace,
    lambda_beta = 0,
    lambda_gamma = 0,
    data = data,
    intercept_row = intercept_row,
    intercept_col = intercept_col,
    hpar = hpar,
    shared_information = shared_information,
    error_function = error_function,
    thresh = 1e-3,
    trace = trace,
    maxit = maxit,
    ls_initial = ls_initial,
    seed = seed,
    warm_start = warm_start
  )
) %>% print()
  warm_start = onfitf$fit

}
