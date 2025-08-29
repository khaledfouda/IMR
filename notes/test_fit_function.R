library(devtools)
clean_dll(); Rcpp::compileAttributes(); document(); load_all()
require(tidyverse)


f <- function() {
  dlls <- Filter(function(d) d[["name"]] == "IMR", getLoadedDLLs())
  vapply(dlls, `[[`, "", "path")  %>% print()
  getDLLRegisteredRoutines(dlls[[1]])$`.Call` %>% print()
  getDLLRegisteredRoutines(dlls[[1]])$Fortran     %>% print()
}
f()

source("./notes/generate_simu_dat.R")
dat <-
generate_simulated_data(300, 400, 3, 5, 0, 0.7,
                        prepare_for_fitting = T,mv_coeffs = T,seed = 2025)
x <- matrix(rbinom(150, 1, .4), 15,10); x[1,5] <- x[3:4,7] <- NA; as(x, "Incomplete")

fit <- IMR::fit(dat$fit_data$train, dat$X, J=3, lambda_M = .001, lambda_beta=.0001,
                trace=T, ls_initial = FALSE)


col_means_cpp(dat$fit_data$Y_full, 400)



microbenchmark::microbenchmark(
a = partial_crossprod(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
b =  partial_crossprod_cpp(dat$X, (dat$beta), dat$fit_data$train@i, dat$fit_data$train@p,F),
times = 500
)

pcr
all(
==
)
