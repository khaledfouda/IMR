generate_simulated_data <- function(
    n                    = 300,
    m                    = 400,
    J                    = 10,
    p                    = 6,
    q                    = 6,
    missing_prob         = 0.8,
    collinear            = FALSE,
    half_discrete        = FALSE,
    prop_inform_x        = 1,
    prop_inform_z        = 1,
    prepare_for_fitting  = FALSE,
    mar_sparse           = FALSE,
    mv_coeffs              = TRUE,
    seed                 = NULL
) {
  require(utils)
  require(tidyverse)
  # 1) Reproducibility ---------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  Z <- gamma <- NULL
  # 2) Basic argument checks (non-intrusive guards) ----------------------------
  if (missing_prob < 0 || missing_prob > 1) {
    stop("`missing_prob` must be in [0, 1].")
  }
  if (prop_inform_x < 0 || prop_inform_x > 1) {
    stop("`informative_cov_prop` must be in [0, 1].")
  }
  if (collinear && p < 2) {
    stop("`collinear = TRUE` requires k >= 2.")
  }

  # 3) Simulate covariates X and Z ---------------------------------------------
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if(q >0)
    Z <- matrix(rnorm(m * q), nrow = m, ncol = q)

  if (half_discrete) {
    n_disc    <- ceiling(p / 2)
    disc_cols <- sample(seq_len(p), n_disc)
    X[, disc_cols] <- matrix(
      rbinom(n * n_disc, size = 1, prob = 0.3),
      nrow = n, ncol = n_disc
    )
  }

  if (collinear) {
    # Make column 2 a noisy near-copy of column 1 to induce high correlation
    X[, 2] <- X[, 1] + rnorm(n, mean = 0, sd = 1.5)
    if(q >0)
      Z[, 2] <- Z[, 1] + rnorm(m, mean = 0, sd = 1.5)
  }

  # 4) Simulate coefficient matrices beta and gamma ----------------------------
  if (mv_coeffs) {
    beta_means <- runif(p, 1, 3) * sample(c(-1, 1), p, replace = TRUE)
    gamma_means <- runif(q, 1, 3) * sample(c(-1,1), p, replace = TRUE)

    beta_vars  <- runif(p, 0.5, 1)^2
    gamma_vars <- runif(q, 0.5, 1)^2

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
    beta  <- matrix(runif(p * m, 1, 2), nrow = p, ncol = m)
    if(q>0)
      gamma <- matrix(runif(n * q, 1, 2), nrow = n, ncol = q)
  }

  # 5) Low-rank structure M orthogonal to col(X) -------------------------------
  # Projector onto col(X): P_X = X (X'X)^{-1} X'
  # Then P_perp = I - P_X, and M = P_perp U V with U ∈ R^{n×r}, V ∈ R^{r×m}
  U     <- matrix(runif(n * J), nrow = n, ncol = J)
  V     <- matrix(runif(J * m), nrow = J, ncol = m)
  P_X   <- X %*% solve(crossprod(X)) %*% t(X)
  P_perp <- diag(n) - P_X
  M     <- P_perp %*% U %*% V

  # 6) Missingness mask W (1 = observed, 0 = missing) --------------------------
  W <- matrix(
    rbinom(n * m, size = 1, prob = 1 - missing_prob),
    nrow = n, ncol = m
  )

  # 7) Enforce proportion of informative covariates  ---------------------------
  if (prop_inform_x < 1) {
    total_beta <- length(beta)

    if (mar_sparse) {
      # Entrywise sparsity: randomly zero individual beta entries
      to_zero <- sample(seq_len(total_beta),
                        size = round((1 - prop_inform_x) * total_beta))
      beta[to_zero] <- 0
    } else {
      # Rowwise sparsity: keep only a proportion of informative covariate rows
      n_keep <- round(prop_inform_x * p)
      if (n_keep <= 0) {
        beta[] <- 0
      } else if (n_keep < p) {
        remove_idx <- sample(seq_len(p), p - n_keep)
        beta[remove_idx, ] <- 0
      }
    }
  }
  # fill it later >>
  if (prop_inform_z < 100) {
    total_beta <- length(beta)

    if (mar_sparse) {
      # Entrywise sparsity: randomly zero individual beta entries
      to_zero <- sample(seq_len(total_beta),
                        size = round((1 - prop_inform_x) * total_beta))
      beta[to_zero] <- 0
    } else {
      # Rowwise sparsity: keep only a proportion of informative covariate rows
      n_keep <- round(prop_inform_x * p)
      if (n_keep <= 0) {
        beta[] <- 0
      } else if (n_keep < p) {
        remove_idx <- sample(seq_len(p), p - n_keep)
        beta[remove_idx, ] <- 0
      }
    }
  }

  # 8) Observed data with Gaussian noise ---------------------------------------
  THETA        <- X %*% beta  + M
  if(q >0)
    THETA <- THETA + gamma %*% t(Z)
  noise_sd <- 1
  E        <- matrix(rnorm(n * m, mean = 0, sd = noise_sd), nrow = n, ncol = m)
  Y        <- (THETA + E) * W
  rank_theta   <- qr(THETA)$rank

  # 9) Optional train/validation split & helpers -------------------------------
  fit_data <- NULL
  if (prepare_for_fitting) {
    # W_valid: 1 = train (kept), 0 = held-out validation entries
    W_valid <- MC_train_test_split(W, testp = 0.2)
    Y_train <- Y * W_valid

    Xq <- qr(X)
    if(q>0)
      Zq <- qr(Z)
    fit_data <- list(
      train   = as(Y_train, "Incomplete"),
      valid   = Y[W_valid == 0],
      Y_full  = as(Y, "Incomplete"),
      W_valid = W_valid,
      X       = list(Q = qr.Q(Xq), R = qr.R(Xq)),
      Rbeta   = qr.R(Xq) %*% beta
    )
    if(q>0)
      fit_data$Z       = list(Q = qr.Q(Zq), R = qr.R(Zq))
    if(q>0)
      fit_data$gammaRt   =  gamma %*% t(qr.R(Zq))

  }

  # 10) Return -----------------------------------------------------------------
  list(
    theta    = THETA,
    W        = W,
    X        = X,
    Z        = Z,
    gamma    = gamma,
    Y        = Y,
    beta     = beta,
    M        = M,
    rank     = rank_theta,
    fit_data = fit_data
  )
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
