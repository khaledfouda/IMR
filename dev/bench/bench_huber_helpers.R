# ---------------------------------------------------------------------------
#  Benchmarks for the Huber helpers.  Scratch script: dev/ is .Rbuildignore'd,
#  so `bench` is not a package dependency.
#
#  RUN AGAINST AN INSTALLED BUILD, NOT devtools::load_all().
#      devtools::install()          # or:  R CMD INSTALL --preclean .
#  then, in a FRESH session, from the package root:
#      library(IMR); source("dev/bench/bench-huber-helpers.R")
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("bench", quietly = TRUE))

# --- provenance guard: load_all() may compile with -O0 -g -----------------
dll <- getLoadedDLLs()[["IMR"]]
if (is.null(dll)) stop("IMR is not loaded. Run library(IMR) first.")
if (!any(startsWith(normalizePath(dll[["path"]]), normalizePath(.libPaths())))) {
  warning("The IMR DLL is not from an installed library:\n  ", dll[["path"]],
    "\nThis looks like a devtools::load_all() build, which may be ",
    "compiled at -O0. Timings would be meaningless. Install first.",
    call. = FALSE, immediate. = TRUE
  )
}
if (!isTRUE(capabilities("profmem"))) {
  message("R built without memory profiling: mem_alloc will be NA.")
}

# --- the kernels are internal, so bind them locally -----------------------
# Not exported in NAMESPACE, and `IMR:::` inside a benchmarked expression
# would charge a namespace lookup to every iteration.
ns <- asNamespace("IMR")
huber_clip_cpp <- get("huber_clip_cpp", ns)
huber_clip_into_cpp <- get("huber_clip_into_cpp", ns)
huber_clip_inplace_cpp <- get("huber_clip_inplace_cpp", ns)
update_huber_c_cpp <- get("update_huber_c_cpp", ns)

# --- data ------------------------------------------------------------------
set.seed(20260810)
m <- 1e7L
yx <- c(rnorm(9.9e6), rcauchy(1e5, scale = 30)) # 1% contamination
m <- length(yx)
cc <- 1.345
buf <- numeric(m) # numeric() zero-fills, so pages are already resident;
# first-touch faults are not charged to iteration 1
cat(sprintf(
  "m = %s observed entries (%.0f MB per vector)\n\n",
  format(m, big.mark = ","), m * 8 / 1024^2
))


# ===========================================================================
#  1. Clipping: where the allocations go
# ===========================================================================
b_mem <- bench::mark(
  R_pmin_pmax = pmin(pmax(yx, -cc), cc),
  cpp_alloc = huber_clip_cpp(yx, cc),
  cpp_into = {
    huber_clip_into_cpp(yx, cc, buf)
    buf
  },
  check = TRUE, # all three must agree elementwise
  memory = TRUE, min_time = 2,
  filter_gc = FALSE, # GC *is* the cost being measured; do not drop it
  min_iterations = 10
)
cat("--- 1. clip: allocation ---\n")
print(b_mem[c("expression", "median", "mem_alloc", "n_gc", "n_itr")])


# ===========================================================================
#  2. Clipping: timing, with profmem off
#     Rprofmem perturbs timings, so the allocation and timing runs are split.
# ===========================================================================
b_time <- bench::mark(
  R_pmin_pmax = pmin(pmax(yx, -cc), cc),
  cpp_alloc = huber_clip_cpp(yx, cc),
  cpp_into = {
    huber_clip_into_cpp(yx, cc, buf)
    buf
  }, min_time = 2,
  check = TRUE, memory = FALSE, filter_gc = FALSE, min_iterations = 20
)
cat("\n--- 2. clip: timing ---\n")
print(b_time[c("expression", "min", "median", "itr/sec")])
cat("\nrelative to fastest:\n")
print(summary(b_time, relative = TRUE)[c("expression", "median")])


# ===========================================================================
#  3. Clipping: scaling in m
#     The allocation penalty only bites once 8m leaves L3 and GC engages.
# ===========================================================================
b_scale <- bench::press(
  n = c(1e4, 1e5, 1e6, 1e7),
  {
    y <- rnorm(n)
    b <- numeric(n)
    bench::mark(
      R_pmin_pmax = pmin(pmax(y, -cc), cc),
      cpp_alloc = huber_clip_cpp(y, cc),
      cpp_into = {
        huber_clip_into_cpp(y, cc, b)
        b
      }, min_time = 2,
      check = TRUE, memory = FALSE, filter_gc = FALSE, min_iterations = 10
    )
  }
)
b_scale$ns_per_elt <- as.numeric(b_scale$median) / b_scale$n * 1e9
cat("\n--- 3. clip: ns per element ---\n")
print(b_scale[c("expression", "n", "ns_per_elt")], n = Inf)


# ===========================================================================
#  4. Tuning constant: R vs exact C++ vs subsampled C++
#     The R arm is the bare expression, not update_huber_c_R(), so R is not
#     charged for validation it need not do.  The sampled arm is approximate
#     by design, so equality is checked to 5% rather than to tolerance.
# ===========================================================================
tol_check <- function(a, b) isTRUE(all.equal(a, b, tolerance = 0.05))

b_c <- bench::mark(
  R_IQR = min(1.345 * stats::IQR(yx) / 1.349, Inf),
  cpp_exact = update_huber_c_cpp(yx, 1.345, Inf, max_sample = 0L),
  cpp_sampled = update_huber_c_cpp(yx, 1.345, Inf),min_time = 2,
  check = tol_check, memory = TRUE, filter_gc = FALSE, min_iterations = 10
)
cat("\n--- 4. tuning constant ---\n")
print(b_c[c("expression", "median", "mem_alloc", "n_gc")])

cat(sprintf("\nsubsampling relative error: %.4f%%\n", 100 * abs(
  update_huber_c_cpp(yx, 1.345, Inf) /
    update_huber_c_cpp(yx, 1.345, Inf, max_sample = 0L) - 1
)))


# ===========================================================================
#  5. OPTIONAL: the in-place variant.
#     It destroys its input, so it cannot share a bench::mark() call: after
#     iteration 1 the data is already clipped and the branches become
#     perfectly predicted.  Timing a fresh copy alone and copy+clip together
#     brackets the cost.  Difference of medians is an ESTIMATE, not a
#     measurement -- do not report it as one.
# ===========================================================================
b_ip <- bench::mark(
  copy_only = {
    z <- yx + 0
    z
  },
  copy_inplace = {
    z <- yx + 0
    huber_clip_inplace_cpp(z, cc)
    z
  }, min_time = 2,
  check = FALSE, memory = FALSE, filter_gc = FALSE, min_iterations = 20
)
cat("\n--- 5. in-place (bracketed) ---\n")
print(b_ip[c("expression", "median")])
cat(sprintf(
  "implied in-place cost ~ %s\n",
  format(diff(as.numeric(b_ip$median)) |> bench::as_bench_time())
))

print(b_scale)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(b_scale)


#=============
# temp: testthat for huber
testthat::test_that("update_huber_c_cpp reproduces IQR(yx)/1.349 and the min() rule", {

  set.seed(2026)

  # ---------------------------------------------------------------------
  # Every residue of n mod 4 exercises a different type-7 interpolation
  # branch (g == 0 vs g > 0, at each quartile).
  # ---------------------------------------------------------------------
  for (n in c( 4L, 5L, 6L, 7L, 8L, 9L, 101L, 1000L, 4001L)) {
    yx <- rnorm(n)
    testthat::expect_equal(
      update_huber_c_cpp(yx, 1.345, Inf),
      update_huber_c_R(yx, 1.345, Inf),
      info = sprintf("type-7 IQR mismatch at n = %d.", n)
    )
    testthat::expect_equal(
      update_huber_c_cpp(yx, 1.345, Inf),
      stats::IQR(yx) / 1.349 * 1.345,
      info = sprintf("does not agree with stats::IQR() at n = %d.", n)
    )
  }

  # Heavy-tailed residuals: the whole point of the robust scale.
  yx <- c(rnorm(5000), rcauchy(200, scale = 50))
  testthat::expect_equal(update_huber_c_cpp(yx, 2, 10), update_huber_c_R(yx, 2, 10))

  # ---------------------------------------------------------------------
  # The min() clamp: monotone non-increasing across iterations.
  # ---------------------------------------------------------------------
  yx <- rnorm(1000)
  d  <- stats::IQR(yx) / 1.349
  testthat::expect_equal(update_huber_c_cpp(yx, 1.345, 1e-6), 1e-6,
                         info = "c_old must dominate when it is the smaller value.")
  testthat::expect_equal(update_huber_c_cpp(yx, 1.345, 1e6), 1.345 * d,
                         info = "the scale must dominate when c_old is large.")
  c_seq <- Reduce(function(cc, i) update_huber_c_cpp(yx, 1.345, cc), 1:5, Inf,
                  accumulate = TRUE)[-1]
  testthat::expect_true(all(diff(c_seq) <= 0),
                        info = "the tuning-constant sequence is not non-increasing.")

  # ---------------------------------------------------------------------
  # Unhappy paths
  # ---------------------------------------------------------------------
  # Fully-missing response: Y@x is numeric(0).
  testthat::expect_identical(update_huber_c_cpp(numeric(0), 1.345, 3), 0)

  # Single observed entry: IQR is 0, hence c collapses to 0.
  testthat::expect_identical(update_huber_c_cpp(5, 1.345, 3), 0)
  testthat::expect_equal(update_huber_c_cpp(5, 1.345, 3), update_huber_c_R(5, 1.345, 3))

  # Degenerate (rank-deficient / exactly-fitted) residuals: IQR == 0.
  testthat::expect_identical(update_huber_c_cpp(rep(2.5, 1000L), 1.345, 3), 0)

  # NA / NaN in the residual vector.
  #testthat::expect_error(update_huber_c_cpp(c(1, NA, 3), 1.345, 3),   "missing values")
  #testthat::expect_error(update_huber_c_cpp(c(1, NaN, 3), 1.345, 3),  "missing values")

  # NA scalars: base::min() would silently return NA, so both must stop.
  testthat::expect_error(update_huber_c_cpp(rnorm(10), NA_real_, 3), "NA/NaN")
  testthat::expect_error(update_huber_c_cpp(rnorm(10), 1.345, NA_real_), "NA/NaN")


  # The input must survive: nth_element operates on a copy.
  yx  <- rnorm(500)
  ref <- yx + 0
  invisible(update_huber_c_cpp(yx, 1.345, Inf))
  testthat::expect_identical(yx, ref, info = "update_huber_c_cpp mutated its input.")
})
testthat::test_that("huber_clip_* agree with the R reference and with each other", {

  set.seed(2027)
  yx <- c(rnorm(10000), rcauchy(500, scale = 30))

  for (cc in c(0, 1e-8, 0.5, 1.345, 1e6, Inf)) {
    expected <- huber_clip_R(yx, cc)

    testthat::expect_equal(huber_clip_cpp(yx, cc), expected,
                           info = sprintf("huber_clip_cpp mismatch at c = %g.", cc))

    into <- numeric(length(yx))
    huber_clip_into_cpp(yx, cc, into)
    testthat::expect_equal(into, expected,
                           info = sprintf("huber_clip_into_cpp mismatch at c = %g.", cc))

  }

  # psi_c is idempotent and bounded.
  cl <- huber_clip_cpp(yx, 1.345)
  testthat::expect_identical(huber_clip_cpp(cl, 1.345), cl)
  testthat::expect_true(max(abs(cl)) <= 1.345)

  # psi_c(u) = min(1, c/|u|) * u, i.e. clipping is the IRLS reweighting.
  w <- pmin(1, 1.345 / abs(yx))
  testthat::expect_equal(cl, w * yx)

  # Non-destructive variants leave the source untouched.
  ref <- yx + 0
  invisible(huber_clip_cpp(yx, 1.345))
  into <- numeric(length(yx)); huber_clip_into_cpp(yx, 1.345, into)
  testthat::expect_identical(yx, ref, info = "the source residual vector was mutated.")

  # ---------------------------------------------------------------------
  # Unhappy paths
  # ---------------------------------------------------------------------
  # Fully-missing response.
  testthat::expect_identical(huber_clip_cpp(numeric(0), 1.345), numeric(0))
  into0 <- numeric(0)
  testthat::expect_silent(huber_clip_into_cpp(numeric(0), 1.345, into0))

  # NA / NaN propagate rather than becoming +/- c.  This verifies the
  # comparison semantics of the generated code on *this* toolchain.
  probe <- c(1, NA_real_, -50, NaN, 50)
  testthat::expect_equal(huber_clip_cpp(probe, 2), huber_clip_R(probe, 2))
  testthat::expect_true(is.na(huber_clip_cpp(probe, 2)[2]))
  testthat::expect_true(is.nan(huber_clip_cpp(probe, 2)[4]))

  # +/- Inf are clipped like any other extreme value.
  testthat::expect_equal(huber_clip_cpp(c(-Inf, 0, Inf), 2), c(-2, 0, 2))

  # Length mismatch on the buffer variant.
  testthat::expect_error(huber_clip_into_cpp(rnorm(10), 1.345, numeric(9)),
                         "same length")
  testthat::expect_error(huber_clip_into_cpp(rnorm(10), 1.345, numeric(11)),
                         "same length")

  # Invalid c.
  for (f in list(huber_clip_cpp)) {
    testthat::expect_error(f(rnorm(10), NA_real_), "NA/NaN")
    testthat::expect_error(f(rnorm(10), -1),       "non-negative")
  }
  testthat::expect_error(huber_clip_into_cpp(rnorm(10), -1, numeric(10)), "non-negative")
})
testthat::test_that("the Huber helpers compose on an actual Incomplete object", {

  # Column 2 and row 3 are entirely unobserved.
  M <- matrix(c(1, 0, -20,
                0, 0,   0,
                4, 0,   0), nrow = 3, byrow = TRUE)
  mat <- as_incomplete(M)

  cc <- update_huber_c_cpp(mat@x, 1.345, Inf)
  testthat::expect_equal(cc, update_huber_c_R(mat@x, 1.345, Inf))
  testthat::expect_false(is.na(cc))

  into <- numeric(length(mat@x))
  huber_clip_into_cpp(mat@x, cc, into)
  testthat::expect_equal(into, huber_clip_R(mat@x, cc))
  testthat::expect_true(max(abs(into)) <= cc)

  # A completely unobserved matrix yields an empty @x slot.
  empty <- as_incomplete(matrix(0, 4, 4))
  testthat::expect_length(empty@x, 0L)
  testthat::expect_identical(update_huber_c_cpp(empty@x, 1.345, 3), 0)
})
testthat::test_that("update_huber_c_cpp subsamples accurately and deterministically", {

  set.seed(4242)
  n  <- 2e6
  yx <- rnorm(n)                                   # true IQR / 1.349 = 1

  exact   <- update_huber_c_cpp(yx, 1, Inf, max_sample = 0L)
  sampled <- update_huber_c_cpp(yx, 1, Inf)        # default max_sample = 1e5

  testthat::expect_equal(exact, stats::IQR(yx) / 1.349,
                         info = "max_sample = 0 must reproduce stats::IQR() exactly.")

  # 0.37% nominal sd at m' = 1e5; 2% is a ~5 sd envelope.
  testthat::expect_lt(abs(sampled - exact) / exact, 0.02)

  # No RNG is involved at all: identical across calls and across seeds.
  testthat::expect_identical(sampled, update_huber_c_cpp(yx, 1, Inf))
  set.seed(1); a <- update_huber_c_cpp(yx, 1, Inf)
  set.seed(2); b <- update_huber_c_cpp(yx, 1, Inf)
  testthat::expect_identical(a, b)

  # R's RNG stream must be untouched.
  set.seed(99); before <- .Random.seed
  invisible(update_huber_c_cpp(yx, 1, Inf))
  testthat::expect_identical(.Random.seed, before)

  # Sampling is effectively without replacement, so error shrinks with m'.
  err <- vapply(c(1e3, 1e4, 1e5), function(m)
    abs(update_huber_c_cpp(yx, 1, Inf, max_sample = as.integer(m)) - exact) / exact,
    numeric(1))
  testthat::expect_lt(err[3], err[1])

  # Heavy-tailed residuals: the robust scale must stay near the Gaussian core.
  set.seed(7)
  yh <- c(rnorm(n * 0.99), rcauchy(n * 0.01, scale = 50))
  ex_h <- update_huber_c_cpp(yh, 1, Inf, max_sample = 0L)
  testthat::expect_lt(abs(update_huber_c_cpp(yh, 1, Inf) - ex_h) / ex_h, 0.04)

  # Small inputs bypass sampling entirely, so exact agreement is preserved.
  small <- rnorm(500)
  testthat::expect_equal(update_huber_c_cpp(small, 1.345, Inf),
                         update_huber_c_R(small, 1.345, Inf))
})

