generate_simulated_data <- function(
    n = 300,
    m = 400,
    r = 10,
    p = 6,
    q = 6,
    sparsity = 0.8,
    sparsity_beta = 1,
    sparsity_gamma = 1,
    snr = 1.5, # Reduced from 0.6 / 0.4
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  #  checks ---------------------------------------------------
  if (sparsity <= 0 || sparsity > 1) stop("`sparsity` must be in (0, 1].")
  if (sparsity_beta < 0 || sparsity_beta > 1) stop("`sparsity_beta` must be in [0, 1].")
  if (sparsity_gamma < 0 || sparsity_gamma > 1) stop("`sparsity_gamma` must be in [0, 1].")

  # Simulate covariates X and Z ---------------------------------------------
  x_mat <- if (p > 0) matrix(runif(n * p), nrow = n, ncol = p) else NULL
  z_mat <- if (q > 0) matrix(runif(m * q), nrow = m, ncol = q) else NULL

  # Simulate coefficient matrices  ------------------------
  beta_mat <- NULL
  if (p > 0) {
    beta_means <- runif(p, 0.1, 1) * sample(c(-1, 1), p, replace = TRUE)
    beta_vars <- runif(p, 0.5, 1)^2
    beta_mat <- t(MASS::mvrnorm(
      n = m,
      mu = beta_means,
      Sigma = diag(beta_vars, nrow = p)
    ))
  }

  gamma_mat <- NULL
  if (q > 0) {
    gamma_means <- runif(q, 0.1, 1) * sample(c(-1, 1), q, replace = TRUE)
    gamma_vars <- runif(q, 0.5, 1)^2
    gamma_mat <- MASS::mvrnorm(
      n = n,
      mu = gamma_means,
      Sigma = diag(gamma_vars, nrow = q)
    )
  }

  # Low-rank structure M  ------------------------
  u_mat <- matrix(runif(n * r, -1, 1), nrow = n, ncol = r)
  v_mat <- matrix(runif(m * r, -1, 1), nrow = m, ncol = r)
  m_mat <- u_mat %*% t(v_mat)

  # Missingness mask (1 = observed, 0 = missing) ----------------------------
  mask <- matrix(rbinom(n * m, size = 1, prob = 1 - sparsity), nrow = n, ncol = m)

  # Enforce sparsity in covariate coefficients --------------------
  if (sparsity_beta > 0 && p > 0) {
    to_zero_b <- sample(seq_len(length(beta_mat)), size = round(sparsity_beta * length(beta_mat)))
    beta_mat[to_zero_b] <- 0
  }

  if (sparsity_gamma > 0 && q > 0) {
    to_zero_g <- sample(seq_len(length(gamma_mat)), size = round(sparsity_gamma * length(gamma_mat)))
    gamma_mat[to_zero_g] <- 0
  }

  # Combine components and generate noise ---------------------------------------
  theta <- m_mat
  if (p > 0) theta <- theta + (x_mat %*% beta_mat)
  if (q > 0) theta <- theta + (gamma_mat %*% t(z_mat))

  noise_sd <- sqrt((sum((theta - mean(theta))^2) / (n * m - 1)) / (snr^2))
  e_mat <- matrix(rnorm(n * m, mean = 0, sd = noise_sd), nrow = n, ncol = m)

  y_mat <- (theta + e_mat) * mask
  if(100 < max(n/2, m/2)){
    rank_theta <- sum(irlba::irlba(theta, nv = 100)$d)
  }else
  rank_theta <- qr(theta)$rank

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
