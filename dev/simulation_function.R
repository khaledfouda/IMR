check_param <- function(...) .imr_check_param(..., raise_error = TRUE)

sample_matrix <- function(n, p, dist, mvnorm_means = NULL, mvnorm_vars = NULL) {
  check_param(dist, "character", choices = c("uniform", "normal", "mvnorm"))
  check_param(c(n,p), "numeric_list", 0, integer = TRUE, min_inclusive = FALSE)
  if (p == 0) {
    return(NULL)
  }
  if (dist == "uniform") {
    return(matrix(stats::runif(n * p), n, p))
  }else if(dist == "normal"){
  return(matrix(stats::rnorm(n * p), n, p))
  }else if(dist == "mvnorm"){
    check_param(mvnorm_means, "numeric_list")
    check_param(mvnorm_vars, "numeric_list", 0, min_inclusive = FALSE, max_inclusive = FALSE)
              stopifnot(length(mvnorm_means) == p &&
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

select_outlier_cells <- function(mask, prop, structure = c("cellwise", "rowwise", "colwise")) {
  n <- nrow(mask)
  obs_idx <- which(mask == 1L)

  if (structure == "cellwise")
    return(obs_idx[sample.int(length(obs_idx),size =  round(prop * length(obs_idx)))])

  n_lines <- if (structure == "rowwise") nrow(mask) else ncol(mask)
  n_out <- round(prop * n_lines)
  lines_out <-  sample.int(n_lines, size = n_out)

  line_of_cell <- if (structure == "rowwise") {
    (obs_idx - 1L) %% n + 1L
  } else {
    (obs_idx - 1L) %/% n + 1L
  }
  return(obs_idx[line_of_cell %in% lines_out])
}

draw_outlier_shifts <- function(n_out, scale_ref, mag, sign) {
  magnitude <- stats::runif(n_out, min = 0, max = mag * scale_ref)
  signs <- switch(
    sign,
    symmetric = sample(c(-1, 1), size = n_out, replace = TRUE),
    positive  = rep(1, n_out),
    negative  = rep(-1, n_out)
  )
  signs * magnitude
}

#
#
# # the following for testing only, remove later.
# n = 300; m = 400; r = 10; p = 6; q = 6; shared_beta = FALSE; shared_gamma=TRUE; sparsity_beta = .2; sparsity_gamma = 0;
# sparse_by_variable = TRUE; orthogonalise = TRUE; signal_share = c(M = 1/3, beta = 1/3, gamma = 1/3);
# signal_scale = 1; sparsity = 0.8; snr = 1.5; outlier_prop = 0; outlier_mag = 10;
# outlier_structure = c("cellwise", "rowwise", "colwise"); outlier_law = c("shift", "scale", "cauchy");
# outlier_sign = c("symmetric", "positive", "negative"); outlier_within = 1
# missing_mechanism  = "mar"; missing_rate = .8

simulate.d <- function(
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
    sparse_by_variable = FALSE,
    # settings for orthogonality between M, Xbeta, and Xgamma
    orthogonalise = TRUE,
    # signal proportion between the components
    signal_share = c(M = 1/3, beta = 1/3, gamma = 1/3),
    # noise for Y
    snr = 1.5, # Reduced from 0.6 / 0.4
    # missing values in Y
    missing_mechanism  = c("mcar", "mar", "mnar_theta", "mnar_y"),
    missing_rate = .8,
    # settings for outliers
    outlier_prop = 0,
    outlier_mag = 10,
    outlier_structure = c("cellwise", "rowwise", "colwise"),
    outlier_sign = c("symmetric", "positive", "negative"),
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  #  checks ---------------------------------------------------
  check_param(n, "numeric", 0, integer = TRUE, min_inclusive = FALSE, max_inclusive = FALSE)
  check_param(m, "numeric", 0, integer = TRUE, min_inclusive = FALSE, max_inclusive = FALSE)
  check_param(p, "numeric", 0, min(n,m), integer = TRUE,  max_inclusive = FALSE)
  check_param(q, "numeric", 0, min(n,m), integer = TRUE,  max_inclusive = FALSE)
  check_param(r, "numeric", 0, min(n,m), integer = TRUE,  max_inclusive = FALSE)
  stopifnot(p>0 || q>0 || r>0)

  check_param(shared_beta, "bool")
  check_param(shared_gamma, "bool")
  check_param(sparsity_beta,  "numeric", 0, 1)
  check_param(sparsity_gamma,  "numeric", 0, 1)
  check_param(sparse_by_variable, "bool")

  check_param(orthogonalise, "bool")
  if(orthogonalise) stopifnot(r <= min(n-p, m-q))

  if(! is.null(signal_share)){
    comp <- c()
    if(r > 0)  comp <- c(comp, "M")
    if(p > 0) comp <- c(comp, "beta")
    if(q > 0) comp <- c(comp, "gamma")
    stopifnot(all(names(signal_share) %in% comp) && all(comp %in% names(signal_share)))
  }
  check_param(snr,  "numeric", 0, max_inclusive = FALSE, min_inclusive = FALSE)


  check_param(missing_rate,  "numeric", 0, 1)
  check_param(missing_mechanism, "character", choices = c("mcar", "mar", "mnar_theta", "mnar_y"))

  check_param(outlier_structure, "character", choices =  c("cellwise", "rowwise", "colwise"))
  check_param(outlier_sign, "character", choices = c("symmetric", "positive", "negative"))
  check_param(outlier_prop,  "numeric", 0, 1)
  check_param(outlier_mag,  "numeric", 0, min_inclusive = FALSE, max_inclusive = FALSE)


  # Covariates X and Z ---------------------------------------------
  x_mat <- if(p > 0) sample_matrix(n, p, "uniform") else NULL
  z_mat <- if(q > 0) sample_matrix(m, q, "uniform") else NULL

  # the following is to fix a signal mismatch. shared_beta requires centered Z and vice-versa.
  if (p>0 && q>0 && shared_gamma) x_mat <- sweep(x_mat, 2, colMeans(x_mat), "-")
  if (q>0 && q>0 && shared_beta) z_mat <- sweep(z_mat, 2, colMeans(z_mat), "-")

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
    gamma_vars <-  runif(q, 0.5, 1)^2
    gamma_mat <- if(shared_gamma) matrix(gamma_means, 1, q) else sample_matrix(n, q, "mvnorm", gamma_means, gamma_vars)
  }else
    gamma_mat <- NULL


  # Sparsity in coefficients matrices beta and Gamma --------------------
  if(p > 0)
    beta_mat <- zero_out(beta_mat, sparsity_beta, rowwise = sparse_by_variable)
  if(q > 0)
    gamma_mat <- t(zero_out(t(gamma_mat), sparsity_gamma, rowwise = sparse_by_variable))


  # Signal distribution between M, Xbeta, GammaZ and construct theta -----------
  parts <- list(M = m_mat)
  if(p > 0)
    parts[["beta"]] <-  x_mat %*% beta_mat
  if(q > 0)
    parts[["gamma"]] <-  gamma_mat %*% t(z_mat)

  if(p > 0 && shared_beta) parts[["beta"]] <- matrix(as.vector(parts[["beta"]]), n, m, byrow = FALSE)
 if(q > 0 && shared_gamma) parts[["gamma"]] <- matrix(as.vector(t(parts[["gamma"]])), n, m, byrow = TRUE)

  parts_norm <- vapply(parts, norm, 1, type="F")

  if(! is.null(signal_share)){

    signal_share <- signal_share[names(parts)]
    signal_share <- signal_share / sum(signal_share)
    scale_fac <- sqrt(signal_share * n * m) / parts_norm
  } else{
    comp <- c()
    if(r > 0)  comp <- c(comp, "M")
    if(p > 0) comp <- c(comp, "beta")
    if(q > 0) comp <- c(comp, "gamma")
    scale_fac <- rep(1, length(comp))
    names(scale_fac) <- comp
  }

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
      outlier_idx <- select_outlier_cells(mask, outlier_prop, outlier_structure)

      outlier_shift <- draw_outlier_shifts(length(outlier_idx), noise_sd, outlier_mag,  outlier_sign)

      y_mat[outlier_idx] <- y_mat[outlier_idx] + outlier_shift
      outlier_mask <- matrix(0L, n, m)
      outlier_mask[outlier_idx] <- 1L
    }else{
      outlier_mask <- NULL
      outlier_shift <- NULL
      outlier_idx <- NULL
    }
    # Done. Diagnostic information -----------------------------
    rel <- function(a, b) a / b

    theta_fro <- norm(theta, "F")
    diagnostics <- list(
      signal_share_target = signal_share,
      signal_share_realised = c(
        M = if(r > 0) rel(norm(parts$M, "F")^2, theta_fro^2) else NA_real_,
        beta = if(p > 0) rel(norm(parts$beta, "F")^2, theta_fro^2) else NA_real_,
        gamma = if(q > 0) rel(norm(parts$gamma, "F")^2, theta_fro^2) else NA_real_
      ),
      theta_rms = sqrt(sum(theta^2) / (n * m)),
      theta_sd = sd(theta),
      rank_M = if (r > 0) imr_compute_rank(m_mat) else 0L,
      rank_theta = imr_compute_rank(theta),
      singular_values_M = d_vec,

      orth_XM = if (r > 0 && p > 0) rel(norm(crossprod(x_mat, m_mat), "F"), theta_fro) else NA_real_,
      orth_MZ = if (r > 0 && q > 0) rel(norm(m_mat %*% z_mat, "F"), theta_fro) else NA_real_,

      noise_sd = noise_sd,

      missing_rate_realised = 1 - mean(mask),
      mean_prob_obs = mean(prob_obs),
      outlier_prop_realised = rel(length(outlier_idx), sum(mask))
    )

  # Output ------------------------------------------
    out <- list(
      theta = theta,
      M = m_mat,
      mask = mask,
      prob_obs = prob_obs,
      Y = y_mat,
      diagnostics = diagnostics,
      call = match.call()
    )
    if (p > 0) {
      out$X <- x_mat
      out$beta <- beta_mat
    }
    if (q > 0) {
      out$Z <- z_mat
      out$gamma <- gamma_mat
    }
    if (r > 0) {
      out$U <- u_mat
      out$V <- v_mat
      out$d <- d_vec
    }
    if (!is.null(outlier_mask)) {
      out$Y_clean <- y_clean
      out$outlier_mask <- outlier_mask
      out$outlier_idx <- outlier_idx
      out$outlier_shift <- outlier_shift
    }
    out
}



o = generate_simulated_data(
    # Dimensions
  n = 800,
  m = 900,
  r = 10,
  p = 3,
  q = 2,
  # settings for beta and gamma
  shared_beta = F,
  shared_gamma = F,
  sparsity_beta = 0,
  sparsity_gamma = 0.3,
  sparse_by_variable = TRUE,
  # settings for orthogonality between M, Xbeta, and Xgamma
  orthogonalise = TRUE,
  # signal proportion between the components
  signal_share = c(M = 4/3, beta = 1/3, gamma=1/3),
  # noise for Y
  snr = 1.5, # Reduced from 0.6 / 0.4
  # missing values in Y
  missing_mechanism  = "mar",
  missing_rate = .8,
  # settings for outliers
  outlier_prop = 0.7,
  outlier_mag = 10,
  outlier_structure = "rowwise", #c("cellwise", "rowwise", "colwise"),
  outlier_sign = "positive", #c("symmetric", "positive", "negative"),
  seed = NULL)
o$diagnostics


