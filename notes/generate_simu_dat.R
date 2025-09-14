generate_simulated_data <- function(
    n                    = 300,
    m                    = 400,
    r                    = 10,
    p                    = 6,
    q                    = 6,
    missp                = 0.8,
    collinear            = FALSE,
    #half_discrete        = FALSE,
    sparsity_beta        = 1,
    sparsity_gamma        = 1,
    prepare_for_fitting  = FALSE,
    cov_sparsity_mar      = TRUE,
    mv_coeffs              = TRUE,
    seed                 = NULL
) {
  require(utils)
  require(tidyverse)
  # 1) Reproducibility ---------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  Z <- gamma <- NULL
  # 2) Basic argument checks (non-intrusive guards) ----------------------------
  if (missp < 0 || missp > 1) {
    stop("`missp` must be in [0, 1].")
  }
  if (sparsity_beta < 0 || sparsity_beta > 1) {
    stop("`sparsity_beta` must be in [0, 1].")
  }
  if (sparsity_gamma < 0 || sparsity_gamma > 1) {
    stop("`sparsity_gamma` must be in [0, 1].")
  }
  if (collinear && p < 2) {
    stop("`collinear = TRUE` requires k >= 2.")
  }

  # 3) Simulate covariates X and Z ---------------------------------------------
  if(p > 0)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if(q >0)
    Z <- matrix(rnorm(m * q), nrow = m, ncol = q)

  # if (half_discrete) {
  #   n_disc    <- ceiling(p / 2)
  #   disc_cols <- sample(seq_len(p), n_disc)
  #   X[, disc_cols] <- matrix(
  #     rbinom(n * n_disc, size = 1, prob = 0.3),
  #     nrow = n, ncol = n_disc
  #   )
  # }

  # if (collinear) {
  #   # Make column 2 a noisy near-copy of column 1 to induce high correlation
  #   X[, 2] <- X[, 1] + rnorm(n, mean = 0, sd = 1.5)
  #   if(q >0)
  #     Z[, 2] <- Z[, 1] + rnorm(m, mean = 0, sd = 1.5)
  # }

  # 4) Simulate coefficient matrices beta and gamma ----------------------------
  if (mv_coeffs) {
    # simulate the means from unif[1,3] U unif[-1,-3]
    if(p>0) beta_means <- runif(p, 1, 3) * sample(c(-1, 1), p, replace = TRUE)
    if(q>0) gamma_means <- runif(q, 1, 3) * sample(c(-1,1), q, replace = TRUE)

    # simulate standard deviation from  unif[0.5,1]
    if(p>0) beta_vars  <- runif(p, 0.5, 1)^2
    if(q>0) gamma_vars <- runif(q, 0.5, 1)^2

    # simulate from mv normal distributions.
    if(p>0)
    beta <- t(MASS::mvrnorm(
      n     = m,
      mu    = beta_means,
      Sigma = diag(beta_vars, nrow = p)
    ))
    if(q>0)
      gamma <- MASS::mvrnorm(
        n     = n,
        mu    = gamma_means,
        Sigma = diag(gamma_vars, nrow = q)
      )

  } else {
    if(p>0)
      beta  <- matrix(runif(p * m, 1, 2), nrow = p, ncol = m)
    if(q>0)
      gamma <- matrix(runif(n * q, 1, 2), nrow = n, ncol = q)
  }
  #=====================================================================
  # 5) Low-rank structure M orthogonal to col(X) -------------------------------
  # Projector onto col(X): P_X = X (X'X)^{-1} X'
  # Then P_perp = I - P_X, and M = P_perp U V with U ∈ R^{n×r}, V ∈ R^{r×m}
  U     <- matrix(runif(n * r), nrow = n, ncol = r)
  V     <- matrix(runif(m * r), nrow = m, ncol = r)
  # we now make sure that the column space of X and U are othogonal
  # and the row space of Z and V are orthogonal.
  if(p>0){
    qrx <- qr(X)
    qrx.Q <- qr.Q(qrx)
    U <- (diag(1,n,n) - qrx.Q%*%t(qrx.Q)) %*% U
  }
  if(q>0){
    qrz <- qr(Z)
    qrz.Q <- qr.Q(qrz)
    V <- (diag(1,m,m) - qrz.Q%*%t(qrz.Q)) %*% V
  }
  M <- U %*% t(V)
  # P_X   <- X %*% solve(crossprod(X)) %*% t(X)
  # P_perp <- diag(n) - P_X
  # M     <- P_perp %*% U %*% V

  # 6) Missingness mask (1 = observed, 0 = missing) --------------------------
  mask <- matrix(
    rbinom(n * m, size = 1, prob = 1 - missp),
    nrow = n, ncol = m
  )

  # 7) Enforce proportion of informative covariates  ---------------------------
  if (sparsity_beta < 1 & p > 0) {
    total_beta <- length(beta)

    if (mar_sparse) {
      # Entrywise sparsity: randomly zero individual beta entries
      to_zero <- sample(seq_len(total_beta),
                        size = round(sparsity_beta * total_beta))
      beta[to_zero] <- 0
    } else {
      # Rowwise sparsity: keep only a proportion of informative covariate rows
      n_keep <- round( (1-sparsity_beta) * p)
      if (n_keep <= 0) {
        beta[] <- 0
      } else if (n_keep < p) {
        remove_idx <- sample(seq_len(p), p - n_keep)
        beta[remove_idx, ] <- 0
      }
    }
  }
  # fill it later >>
  if (sparsity_gamma < 1 & q>0) {
    total_gamma <- length(gamma)

    if (mar_sparse) {
      # Entrywise sparsity: randomly zero individual beta entries
      to_zero <- sample(seq_len(total_gamma),
                        size = round((1 - sparsity_gamma) * total_gamma))
      gamma[to_zero] <- 0
    } else {
      # Rowwise sparsity: keep only a proportion of informative covariate rows
      n_keep <- round( (1-sparsity_gamma) * q)
      if (n_keep <= 0) {
        gamma[] <- 0
      } else if (n_keep < q) {
        remove_idx <- sample(seq_len(q), q - n_keep)
        gamma[remove_idx, ] <- 0
      }
    }
  }

  # 8) Observed data with Gaussian noise ---------------------------------------
  THETA        <- M
  if(p>0)  THETA <- THETA + X %*% beta
  if(q>0) THETA <- THETA + gamma %*% t(Z)
  # generate noise
  noise_sd <- 1
  E        <- matrix(rnorm(n * m, mean = 0, sd = noise_sd), nrow = n, ncol = m)
  Y        <- (THETA + E) * mask
  rank_theta   <- qr(THETA)$rank

  # 9) Optional train/validation split & helpers -------------------------------
  fit_data <- NULL
  if (prepare_for_fitting) {
    # W_valid: 1 = train (kept), 0 = held-out validation entries
    W_valid <- MC_train_test_split(mask, testp = 0.2)
    Y_train <- Y * W_valid

    fit_data <- list(
      train   = as(Y_train, "Incomplete"),
      valid   = Y[W_valid == 0],
      Y_full  = as(Y, "Incomplete"),
      W_valid = W_valid
      # X       = list(Q = qr.Q(Xq), R = qr.R(Xq)),
      # Rbeta   = qr.R(Xq) %*% beta
    )
    if(p>0){
      qrx.R <- qr.R(qrx)
      fit_data$X <- list(Q=qrx.Q, R=qrx.R)
      fit_data$Rbeta <- qrx.R %*% beta
    }
    if(q>0)
    {
      qrz.R <- qr.R(qrz)
      fit_data$Z       = list(Q = qrz.Q, R = qrz.R)
      fit_data$gammaRt   =  gamma %*% t(qrz.R)
    }
  }

  # 10) Return -----------------------------------------------------------------
  out <- list(
    theta    = THETA,
    mask        = mask,
    # X        = X,
    # Z        = Z,
    # gamma    = gamma,
    Y        = Y,
    # beta     = beta,
    M        = M,
    rank     = rank_theta,
    fit_data = fit_data
  )
  if(p>0){
    out$X <- X
    out$beta <- beta
  }
  if(q>0){
    out$Z <- Z
    out$gamma <- gamma
  }
  return(out)
}
MC_train_test_split <-
  function(obs_mask,
           testp = 0.2,
           seed = NULL) {
    #' returns a new mask similar to mask with a new train-test sets.
    #' testp: proportion of test out of the nonzero cells in mask
    #' new_mask:
    #------------------------------------------------------------------
    #  1 -> train    (obs_mask=1) |
    #  0 -> test     (obs_mask=1) |
    #  1 -> missing  (obs_mask=0)
    #-----------------------------------------------------------------
    if(!is.null(seed)) set.seed(seed)
    n_rows <- dim(obs_mask)[1]
    n_cols <- dim(obs_mask)[2]
    # Create a data frame of all matrix indices
    # we only consider non-missing data (ie, with mask_ij=1)
    indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)[obs_mask == 1, ]
    # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
    indices <-  indices[sample(1:nrow(indices)),]
    row.names(indices) <- NULL

    test.idx = sample(1:nrow(indices),
                      size = nrow(indices) * testp,
                      replace = FALSE)
    test.indices = indices[test.idx, ]

    new_mask = matrix(1, nrow = n_rows, ncol = n_cols)
    new_mask[obs_mask == 0] <- 1
    new_mask[as.matrix(test.indices[, c("row", "col")])] <- 0

    return(new_mask)
  }

quick_camc_simu_res <- function(
    dat,
    fit,
    test_error  = error_metric$rmse
) {
  # Split train vs test
  M.estimates <- fit$u %*% (fit$d * t(fit$v))
  beta.estim <- tryCatch( solve(dat$fit_data$X$R) %*% fit$beta,
                          error = function(e) 0)
  gamma.estim <- tryCatch( fit$gamma %*% solve(t(dat$fit_data$Z$R)),
                           error = function(e) 0)



  estimates <- M.estimates
  estimates <- tryCatch(estimates + dat$fit_data$X$Q %*% fit$beta,
                        error = function(e) estimates)
  estimates <- tryCatch(estimates + fit$gamma %*% t(dat$fit_data$Z$Q),
                        error = function(e) estimates)
  estimates <- tryCatch(estimates + fit$phi.a %*% matrix(1,1,ncol(M.estimates)),
                        error = function(e) estimates)
  estimates <- tryCatch(estimates + matrix(1,nrow(M.estimates),1) %*% fit$phi.b,
                        error = function(e) estimates)

  mask <- as.matrix(dat$fit_data$Y_full == 0)
  estim.test  <- estimates[mask]
  estim.train <- estimates[mask == 0]
  obs.test    <- dat$theta[mask]
  obs.train   <- dat$theta[mask == 0]

  # Core metrics
  results <- list(
    error.test = test_error(estim.test, obs.test),
    corr.test  = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    error.M    = tryCatch(
      test_error(M.estimates, dat$M),
      error = function(e) NA
    ),
    error.beta = tryCatch(
      test_error(beta.estim, dat$beta),
      error = function(e) NA
    ),
    error.gamma = tryCatch(
      test_error(gamma.estim, dat$gamma),
      error = function(e) NA
    ),
    rank_M     = tryCatch(
      qr(M.estimates)$rank,
      error = function(e) NA
    ),
    rank_beta  = tryCatch(
      qr(beta.estim)$rank,
      error = function(e) NA
    ),
    rank_gamma  = tryCatch(
      qr(gamma.estim)$rank,
      error = function(e) NA
    ),
    sparse_in_sparse = tryCatch(
      sum(dat$beta == 0 & beta.estim == 0) /
        (sum(dat$beta == 0) + 1e-17),
      error = function(e) NA
    ),
    nonsparse_in_nonsparse = tryCatch(
      sum(dat$beta != 0 & beta.estim != 0) /
        (sum(dat$beta != 0) + 1e-17),
      error = function(e) NA
    ),
    sparse_all_beta = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    ),
    sparse_all_gamma = tryCatch(
      sum(gamma.estim == 0) / length(gamma.estim),
      error = function(e) NA
    )
  ) %>% as.data.frame() %>% round(4)
  results
}
