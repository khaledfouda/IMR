sample_matrix <- function(n, p, dist, mvnorm_means = NULL, mvnorm_vars = NULL) {
  stopifnot(.imr_check_param(dist, "character", choices = c("uniform", "normal", "mvnorm")),
            .imr_check_param(c(n,p), "numeric_list", 0, integer = TRUE, min_inclusive = FALSE))
  if (p == 0) {
    return(NULL)
  }
  if (dist == "uniform") {
    return(matrix(stats::runif(n * p), n, p))
  }else if(dist == "normal"){
  return(matrix(stats::rnorm(n * p), n, p))
  }else if(dist == "mvnorm"){
    stopifnot(.imr_check_param(mvnorm_means, "numeric_list") &&
              .imr_check_param(mvnorm_vars, "numeric_list", 0, min_inclusive = FALSE,
                               max_inclusive = FALSE) &&
              length(mvnorm_means) == p &&
              length(mvnorm_vars) == p)
    MASS::mvrnorm(
      n = n,
      mu = mvnorm_means,
      Sigma = diag(mvnorm_vars, nrow = p, ncol = p)
    )
  }
}

orthonormal_matrix <- function(n, r, null_of=NULL){
  A <- sample_matrix(n, r, "normal")

  # if orthogonality required, set A <- A - A * (projection of null_of)
  if(! is.null(null_of))
    A <- qr.resid(qr(null_of), A)

  # finally, we want A^TA=I with rank r
  qr_A <- qr(A)
  if(qr_A$rank < r) {
    stop("Rank of A is smaller than r")
  }
  return(qr.Q(qr_A))
}


zero_out <- function(A, frac, rowwise = FALSE) {
  if (frac <= 0) {
    return(A)
  }
  if (rowwise) {
    rows <- sample.int(nrow(A), size = round(frac * nrow(A)))
    A[rows, ] <- 0
  } else {
    A[sample.int(length(A), size = round(frac * length(A)))] <- 0
  }
  A
}

#' Solve for the logit intercept giving a target mean propensity
#' @noRd
calibrate_logit_intercept <- function(lin_pred, target_mean, max_eval = 1e6) {
  lp <- if (length(lin_pred) > max_eval) sample_vec(lin_pred, max_eval) else lin_pred
  stats::uniroot(
    f = function(a) mean(stats::plogis(a + lp)) - target_mean,
    interval = c(-40, 40),
    extendInt = "yes",
    tol = .Machine$double.eps^0.5
  )$root
}


make_propensity <- function(mechanism, missing_rate, strength,
                            n, m, x_mat, z_mat, theta, y_full) {
  target <- 1 - missing_rate
  if (mechanism == "mcar" || strength == 0) {
    return(matrix(target, n, m))
  }
  lin <- switch(
    mechanism,
    mar = {
      if (is.null(x_mat) && is.null(z_mat)) {
        stop("`missing_mechanism = \"mar\"` requires p > 0 or q > 0.")
      }
      row_part <- if (is.null(x_mat)) rep(0, n) else as.vector(x_mat %*% stats::rnorm(ncol(x_mat)))
      col_part <- if (is.null(z_mat)) rep(0, m) else as.vector(z_mat %*% stats::rnorm(ncol(z_mat)))
      outer(row_part, col_part, "+")
    },
    mnar_theta = theta,
    mnar_y     = y_full
  )
  lin <- strength * (lin - mean(lin)) / stats::sd(as.vector(lin))
  matrix(stats::plogis(calibrate_logit_intercept(as.vector(lin), target) + as.vector(lin)), n, m)
}


# the following for testing only, remove later.
n = 300; m = 400; r = 10; p = 6; q = 6; shared_beta = FALSE; shared_gamma=TRUE; sparsity_beta = .2; sparsity_gamma = 0;
coefs_rowwise_sparsity = TRUE; orthogonalise = TRUE; signal_share = c(M = 1/3, beta = 1/3, gamma = 1/3);
signal_scale = 1; sparsity = 0.8; snr = 1.5; outlier_prop = 0; outlier_mag = 10;
outlier_structure = c("cellwise", "rowwise", "colwise"); outlier_law = c("shift", "scale", "cauchy");
outlier_sign = c("symmetric", "positive", "negative"); outlier_within = 1
missing_mechanism  = "mar"; missing_rate = .8

generate_simulated_data <- function(
    # Dimensions
    n = 300,
    m = 400,
    r = 10,
    p = 6,
    q = 6,
    # settings for beta and gamma
    shared_beta = FALSE,
    shared_gamma = FALSE,
    sparsity_beta = 0,
    sparsity_gamma = 0,
    coefs_rowwise_sparsity = FALSE,
    # settings for orthogonality between M, Xbeta, and Xgamma
    orthogonalise = TRUE,
    signal_share = c(M = 1/3, beta = 1/3, gamma = 1/3),
    signal_scale = 1,
    # Sparsity + noise for Y
    sparsity = 0.8,
    snr = 1.5, # Reduced from 0.6 / 0.4
    missing_mechanism  = "mcar",
    missing_rate = .8,
    # settings for outliers
    outlier_prop = 0,
    outlier_mag = 10,
    outlier_structure = c("cellwise", "rowwise", "colwise"),
    outlier_law = c("shift", "scale", "cauchy"),
    outlier_sign = c("symmetric", "positive", "negative"),
    outlier_within = 1,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  #  checks ---------------------------------------------------
  if (sparsity <= 0 || sparsity > 1) stop("`sparsity` must be in (0, 1].")
  if (sparsity_beta < 0 || sparsity_beta > 1) stop("`sparsity_beta` must be in [0, 1].")
  if (sparsity_gamma < 0 || sparsity_gamma > 1) stop("`sparsity_gamma` must be in [0, 1].")

  if (outlier_prop < 0 || outlier_prop > 1) stop("`outlier_prop` must be in [0, 1].")
  if (outlier_mag < 0) stop("`outlier_mag` must be non-negative.")
  if (outlier_within <= 0 || outlier_within > 1) stop("`outlier_within` must be in (0, 1].")
  outlier_structure <- match.arg(outlier_structure)
  outlier_law <- match.arg(outlier_law)
  outlier_sign <- match.arg(outlier_sign)

  # Covariates X and Z ---------------------------------------------
  x_mat <- if(p > 0) sample_matrix(n, p, "uniform") else NULL
  z_mat <- if(q > 0) sample_matrix(m, q, "uniform") else NULL

  # Low-rank structure M  ------------------------
  u_mat <- orthonormal_matrix(n, r, if(orthogonalise) x_mat else NULL)
  v_mat <- orthonormal_matrix(m, r, if(orthogonalise) z_mat else NULL)
  m_mat <- u_mat %*% t(v_mat)

  # Coefficient matrices beta and Gamma and ------------------------

  if(p > 0){
    beta_means <- runif(p, 0.1, 1) * sample(c(-1, 1), p, replace = TRUE)
    beta_vars <-  runif(p, 0.5, 1)^2
    beta_mat <- if(shared_beta) matrix(beta_means, p, 1) else t(sample_matrix(m, p, "mvnorm", beta_means, beta_vars))
  }else
    beta_mat <- NULL


  if(q > 0){
    gamma_means <- runif(q, 0.1, 1) * sample(c(-1, 1), q, replace = TRUE)
    gamma_vars <-  runif(p, 0.5, 1)^2
    gamma_mat <- if(shared_gamma) matrix(gamma_means, 1, q) else sample_matrix(n, q, "mvnorm", gamma_means, gamma_vars)
  }else
    gamma_mat <- NULL


  # Sparsity in coefficients matrices beta and Gamma --------------------
  if(p > 0)
    beta_mat <- t(zero_out(t(beta_mat), sparsity_beta, rowwise = coefs_rowwise_sparsity))
  if(q > 0)
    gamma_mat <- zero_out(gamma_mat, sparsity_gamma, rowwise = coefs_rowwise_sparsity)


  # Signal distribution between M, Xbeta, GammaZ and construct theta -----------
  parts <- list(M = m_mat)
  if(p > 0)
    parts[["beta"]] <-  x_mat %*% beta_mat
  if(q > 0)
    parts[["gamma"]] <-  gamma_mat %*% t(z_mat)

  if(p > 0 && shared_beta) parts[["beta"]] <- matrix(as.vector(parts[["beta"]]), n, m, byrow = FALSE)
  if(q > 0 && shared_gamma) parts[["gamma"]] <- matrix(as.vector(t(parts[["gamma"]])), n, m, byrow = TRUE)

  parts_norm <- vapply(parts, norm, 1, type="F")

  signal_share <- signal_share[names(parts)]
  signal_share[parts_norm < 0] <- 0
  signal_share <- signal_share / sum(signal_share)
  scale_fac <- sqrt(signal_share * n * m) / parts_norm
  if(r > 0) d_vec <- rep(scale_fac[["M"]], r)
  if(p > 0) beta_mat <- beta_mat * scale_fac[["beta"]]
  if(q > 0) gamma_mat <- gamma_mat * scale_fac[["gamma"]]
  parts <- Map(`*`, parts, as.list(scale_fac))
  m_mat <- parts$M
  theta <- parts$M
  if(p > 0) theta <- theta + parts$beta
  if(q > 0) theta <- theta + parts$gamma

  # Noise and construct Y = theta + noise --------------------
  if(snr > 0){
    noise_sd <- sqrt((sum((theta - mean(theta))^2) / (n * m - 1)) / (snr^2))
    e_mat <- matrix(rnorm(n * m, mean = 0, sd = noise_sd), nrow = n, ncol = m)
  }else {
    noise_sd <- 0
    e_mat <- matrix(0, n, m)
  }

  y_full <- theta + e_mat


  # Missingness mask (1 = observed, 0 = missing) and mask Y ----------------------------
  prob_obs <- make_propensity(mechanism = missing_mechanism,
                              missing_rate = missing_rate,
                              strength = 1,
                              n = n, m = m,
                              x_mat = x_mat, z_mat = z_mat,
                              theta = theta, y_full = y_full)

    mask <- matrix(stats::rbinom(n * m, size = 1, prob = as.vector(prob_obs)), nrow = n, ncol = m)
    y_mat <- y_full * mask

    # Outliers, Applied only to Y and only on observed cells
    y_clean <- y_mat

    outlier_mask <- NULL
    if (outlier_prop > 0 && outlier_mag > 0) {
      scale_ref <- if (noise_sd > 0) noise_sd else theta_sd
      outlier_idx <- select_outlier_cells(
        mask, outlier_prop, outlier_structure, outlier_within,
        target = outlier_target,
        lev = if (outlier_structure == "colwise") leverage(z_mat) else leverage(x_mat)
      )
      outlier_shift <- draw_outlier_shifts(
        length(outlier_idx), scale_ref, outlier_mag, outlier_law, outlier_sign
      )
      y_mat[outlier_idx] <- y_mat[outlier_idx] + outlier_shift
      outlier_mask <- matrix(0L, n, m)
      outlier_mask[outlier_idx] <- 1L
    }else{
      outlier_mask <- NULL
      outlier_shift <- 0
    }


  # Combine components and generate noise ---------------------------------------
  theta <- m_mat
  if (p > 0 && !shared) theta <- theta + (x_mat %*% beta_mat)
  if (q > 0 && !shared) theta <- theta + (gamma_mat %*% t(z_mat))
  if (p > 0 && shared) theta <- sweep(theta, 1, as.vector(x_mat %*% beta_mat), "+")
  if (q > 0 && shared) theta <- sweep(theta, 2, as.vector(gamma_mat %*% t(z_mat)), "+")

  if(snr > 0){
    noise_sd <- sqrt((sum((theta - mean(theta))^2) / (n * m - 1)) / (snr^2))
    e_mat <- matrix(rnorm(n * m, mean = 0, sd = noise_sd), nrow = n, ncol = m)
    y_mat <- (theta + e_mat) * mask
  }else {
    noise_sd <- 0
    y_mat <- theta * mask
  }
  rank_theta <- imr_compute_rank(theta, nv = 100)

  #  Done :) ------------------------------------------------------------------
  out <- list(
    theta = theta,
    mask  = mask,
    Y     = y_mat,
    M     = m_mat,
    rank  = rank_theta
  )

  if (p > 0) {
    out$X <- x_mat
    out$beta <- beta_mat
  }
  if (q > 0) {
    out$Z <- z_mat
    out$gamma <- gamma_mat
  }

  out
}

