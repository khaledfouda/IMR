#' The following function, $MCCI.fit()$, estimates $\hat\beta$ and $\hat B$ as above with fixed (given) hyperparameters.
#'  I will try to use the same notations as above to avoid confusion.


# EDIT: These functions return the 1 / (probability of inclusion) NOT missingness.
MCCI_weights <- list(
  binomial = function(X, W, ...) {
    # using logistic regression as indicated in (a1)
    n1 = dim(W)[1]
    n2 = dim(W)[2]
    theta_hat = matrix(NA, n1, n2)
    for (j in 1:n2) {
      model_data = data.frame(cbind(W[, j], X))
      model_fit = glm(X1 ~ ., family = binomial(), data = model_data)
      theta_hat[, j] = 1 / predict(model_fit, type = "response")
    }
    theta_hat[is.na(theta_hat) | is.infinite(theta_hat)] <- 0
    return(theta_hat)
  },
  column_avg = function(W, ...) {
    # A theta estimation function that selects the proportion of missing data within the same column
    # using formula (a2)
    n1 = dim(W)[1]
    n2 = dim(W)[2]
    theta_hat = matrix(NA, n1, n2)
    for (j in 1:n2) {
      theta_hat[, j] = n1 / sum(W[, j] == 1)
    }
    theta_hat[is.na(theta_hat) | is.infinite(theta_hat)] <- 0
    return(theta_hat)
  },
  uniform = function(W, ...) {
    # A theta estimation function that selects the proportion of missing data in the matrix
    # using formula (a3)
    n1 = dim(W)[1]
    n2 = dim(W)[2]
    theta_hat = matrix((n1 * n2) / sum(W == 1) , n1, n2)
    return(theta_hat)
  }

)

MCCI_fit <-
  function(Y,
           X,
           W,
           lambda_1,
           lambda_2,
           alpha,
           n1n2_optimized = TRUE,
           return_rank = TRUE,
           theta_estimator = MCCI_weights$binomial) {
    #
    #' ----------------------------------------------
    #' Input: Y: corrupted, partially observed A, (Y is assumed to be the product of Y*W)
    #'         Important: set missing values in Y to 0
    #'        W: Wij=1 if Aij is observed (ie, Yij !=0), and 0 otherwise
    #'         X: covariate matrix
    #'         lambda_1, lambda_2, alpha: hyperparameters
    #' ----------------------------------------------
    #' output: list of  A, Beta_hat, B_hat
    #' ----------------------------------------------
    stopifnot(is.matrix(Y))
    Y[is.na(Y)] <- 0
    n1 = dim(Y)[1]
    n2 = dim(Y)[2]
    m  = dim(X)[2]
    #yobs = W==1
    # The following two lines are as shown in (c) and (d)
    X.X = t(X) %*% X
    P_X = X %*% MASS::ginv(X.X) %*% t(X)
    P_bar_X = diag(1, n1, n1) - P_X

    if (n1n2_optimized == TRUE) {
      # we define the factor that will be used later:
      n1n2 = svd(X.X)$d[1]
    } else{
      n1n2 = n1 * n2 / 2
    }

    # The following part estimates theta (missingness probabilities)
    theta_hat = theta_estimator(W = W, X = X)
    # the following is the product of W * theta_hat * Y
    W_theta_Y = Y * theta_hat # * W

    # beta hat as (8)
    beta_hat = MASS::ginv(X.X + diag(n1n2 * lambda_1, m, m)) %*% t(X) %*% W_theta_Y
    # SVD decomposition to be used in (b)
    svdd = svd(P_bar_X %*% W_theta_Y)
    if (n1n2_optimized == TRUE) {
      # evaluation of  (b)
      n1n2 = svdd$d[1]
    } else{
      n1n2 = n1 * n2 / 2
    }
    T_c_D = svdd$u %*% (pmax(svdd$d - alpha * n1n2 * lambda_2, 0) * t(svdd$v))
    # B hat as in (11)
    B_hat = T_c_D / (1 + 2 * (1 - alpha) * n1n2 * lambda_2)
    # computing the rank of B [Copied from MCCI's code; Don't understand how it works.]
    # EQUIVALENT to  qr(B_hat)$rank + m   or   qr(A_hat)$rank
    # B is a low rank matrix
    #rank = sum(pmax(svdd$d - n1n2 * lambda_2 * alpha, 0) > 0) + m

    # Estimate the matrix as given in the model at the top
    A_hat = X %*% beta_hat + B_hat
    #A_hat[yobs] <- Y[yobs]
    rank = NULL
    if (return_rank)
      rank = qr(A_hat)$rank

    return(list(
      estimates = A_hat,
      M = B_hat,
      beta = beta_hat,
      rank = rank
    ))
  }

#' Hyperparameter optimization for $\lambda_1,\lambda_2,\alpha$ is done using k-fold (k=5 by default) and grid search.
#'  MCCI's paper optimizes for each parameter separately while fixing the other two.
#'   The explored/recommended the range of (0-2) for $\lambda_1$, (0.1,0.9) for $\lambda_2$, and (0.992,1) for $\alpha$.

prepare_fold_data <-
  function(Y_train,
           Y_valid,
           W_fold,
           W,
           X,
           n1n2_optimized,
           theta_estimator) {
    n1 = dim(Y_train)[1]
    n2 = dim(Y_train)[2]
    m  = dim(X)[2]

    # The following two lines are as shown in (c) and (d)
    X.X = t(X) %*% X
    P_X = X %*% MASS::ginv(X.X) %*% t(X)
    P_bar_X = diag(1, n1, n1) - P_X

    theta_hat = theta_estimator(W = W_fold, X = X)

    #---------
    # The following are partial parts of equations 8 and 11 that don't involve the hyperparameters.
    # this is useful to avoid unneccessary matrix multiplications.
    #----------
    if (n1n2_optimized == TRUE) {
      # this one is for equation 8, the product n1n2 is replace with the Eigen value
      n1n2Im = svd(X.X)$d[1]  * diag(1, m, m)
    } else{
      n1n2Im = n1 * n2  * diag(1, m, m) * 0.5
    }
    # the following is the product of W * theta_hat * Y
    W_theta_Y = Y_train * theta_hat
    X.W.theta.Y = t(X) %*% W_theta_Y
    svdd = svd(P_bar_X %*% W_theta_Y)
    if (n1n2_optimized == TRUE) {
      # this one is for equation 11, the product is also replace with the Eigen value of the SVD
      n1n2 = svdd$d[1]
    } else{
      n1n2 = n1 * n2 / 2
    }

    return(
      list(
        Y_valid = Y_valid,
        W_fold = W_fold,
        W = W,
        X = X,
        X.X = X.X,
        n1n2Im = n1n2Im,
        n1n2 = n1n2,
        X.W.theta.Y = X.W.theta.Y,
        svdd = svdd
      )
    )
  }

MCCI.fit_optimized <- function(data, lambda_1, lambda_2, alpha) {
  beta_hat = MASS::ginv(data$X.X + data$n1n2Im * lambda_1) %*% data$X.W.theta.Y
  T_c_D = data$svdd$u %*% (pmax(data$svdd$d - alpha * data$n1n2 * lambda_2, 0) * t(data$svdd$v))
  # B hat as in (11)
  B_hat = T_c_D / (1 + 2 * (1 - alpha) * data$n1n2 * lambda_2)
  # Estimate the matrix as given in the model at the top
  A_hat = data$X %*% beta_hat + B_hat

  return(A_hat[data$W_fold == 0 & data$W == 1])
}


MCCI.fit_optimized_part1 <- function(data, lambda_1) {
  # returns Xbeta only. Not used.
  beta_hat = MASS::ginv(data$X.X + data$n1n2Im * lambda_1) %*% data$X.W.theta.Y
  Xbeta = data$X %*% beta_hat
  return(Xbeta[data$W_fold == 0 & data$W == 1])
}

MCCI.fit_optimized_part2 <- function(data, lambda_2, alpha) {
  # returns Bhat only. Not used.
  T_c_D = data$svdd$u %*% (pmax(data$svdd$d - alpha * data$n1n2 * lambda_2, 0) * t(data$svdd$v))
  # B hat as in (11)
  B_hat = T_c_D / (1 + 2 * (1 - alpha) * data$n1n2 * lambda_2)
  return(B_hat[data$W_fold == 0 & data$W == 1])
}


MCCI.cv <-
  function(Y,
           X,
           W,
           n_folds = 5,
           lambda_1_grid = seq(0, 1, length = 20),
           lambda_2_grid = seq(0.9, 0.1, length = 20),
           alpha_grid = seq(0.992, 1, length = 20),
           seed = NULL,
           numCores = 1,
           n1n2_optimized = FALSE,
           test_error =IMR::error_metric$rmse,
           theta_estimator = MCCI_weights$binomial,
           sequential = FALSE) {
    start_time <- Sys.time()
    fit_counter <- 0
    #' -------------------------------------------------------------------
    #' Input :
    #' X :  Covariate matrix of size  n1 by m
    #' W : Binary matrix representing the mask. wij=1 if yij is observed. size similar to A
    #' The rest are cross validation parameters
    #' --------------------------------------------------------------------
    #' Output:
    #' list of best parameters and best score (minimum average MSE across folds)
    #' --------------------------------------------------------------------
    Y[is.na(Y)] <- 0 # insure that missing values are set 0
    if (!is.null(seed))
      set.seed(seed = seed)
    #indices = sample(cut(seq(1, nrow(A)), breaks=n_folds, labels=FALSE))
    best_score = Inf
    best_params = list(alpha = NA,
                       lambda_1 = NA,
                       lambda_2 = NA)

    folds <- MC_Kfold_split(nrow(Y), ncol(Y), n_folds, W, seed)

    fold_data = lapply(1:n_folds, function(i) {
      #train_indices = which(indices != i, arr.ind = TRUE)
      W_fold = folds[[i]] #W[train_indices,]
      #---------------------------------------------------------------
      # EDIT: I implemented this above in k_fold_cells() no longer needed
      #W_fold[W==0] = 1 # This to avoid having the missing data as test set.
      # Note that we don't have their original values so if they're passed to the validation step,
      # their original will be equal to 0. We hope we have enough W_fold = 0 while W = 1.
      #---------------------------------------------------
      Y_train = Y
      Y_train[W_fold == 0] <- 0
      Y_valid = Y[W_fold == 0 & W == 1]
      prepare_fold_data(Y_train,
                        Y_valid,
                        W_fold,
                        W,
                        X,
                        n1n2_optimized,
                        theta_estimator)
    })

    # ************************************************************
    if (numCores == 1 & sequential == FALSE) {
      results <-
        foreach(alpha = alpha_grid, .combine = rbind) %:%
        foreach(lambda_2 = lambda_2_grid, .combine = rbind) %do% {
          test_error =IMR::error_metric$rmse
          lambda_1 = 0
          score = 0
          for (i in 1:n_folds) {
            data = fold_data[[i]]
            A_hat_test = MCCI.fit_optimized(data, lambda_1, lambda_2, alpha)
            # Compute the test error using the provided formula
            score = score + test_error(A_hat_test, data$Y_valid)
            # counting the number of fits
            fit_counter = fit_counter + 1
          }
          score = score / n_folds
          c(alpha, lambda_2, score)
        }

      # Process results to find the best parameters
      min_score <- min(results[, 3])

      # Subset to only include results with the minimum score
      min_results <-
        results[results[, 3] == min_score, , drop = FALSE] # Keep it as a dataframe

      # Find the one with the highest lambda_2 in case of multiple results with the same score
      if (nrow(min_results) > 1) {
        best_result <- min_results[which.max(min_results[, 2]),]
      } else {
        best_result <-
          min_results  # If only one row, it's already the best result
      }
      best_params <-
        list(alpha = best_result[1],
             lambda_1 = 0,
             lambda_2 = best_result[2])
      best_score <- best_result[3]


    } else if (numCores == 1 & sequential) {
      # fixing optimal values of lambda 1 and alpha and optimizing for alpha separately
      lambda_1 = 0
      alpha = 1
      best_score = Inf
      for (lambda_2 in lambda_2_grid) {
        score = 0
        for (i in 1:n_folds) {
          data = fold_data[[i]]
          # compute the estimates with a modified fit function
          A_hat_test = MCCI.fit_optimized(data, lambda_1, lambda_2, alpha)
          # -- EDIT: Using MCCI's formula in page 205 to compute the test error
          score = score + test_error(A_hat_test, data$Y_valid)
          # counting the number of fits
          fit_counter = fit_counter + 1
        }
        score = score / n_folds

        if (score < best_score) {
          best_score = score
          best_params$lambda_2 = lambda_2
          #print(paste(score, "lambda_2", lambda_2))
        }
      }
      # fixing optimal values of lambda 2 and lambda 1 and optimizing for alpha separately
      lambda_2 = best_params$lambda_2
      lambda_1 = 0
      best_score = Inf
      for (alpha in alpha_grid) {
        score = 0
        for (i in 1:n_folds) {
          data = fold_data[[i]]
          # compute the estimates with a modified fit function
          A_hat_test = MCCI.fit_optimized(data, lambda_1, lambda_2, alpha)
          # -- EDIT: Using MCCI's formula in page 205 to compute the test error
          score = score + test_error(A_hat_test, data$Y_valid)
          # counting the number of fits
          fit_counter = fit_counter + 1
        }
        score = score / n_folds

        if (score < best_score) {
          best_score = score
          best_params$alpha = alpha
          #print(paste(score, "alpha", alpha))
        }
      }

    } else{
      # Run on multiple cores
      # prepare the cluster
      require(parallel)
      require(doParallel)
      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      # fixing lambda 1 at 0 and optimizing for lambda 2 and alpha using a grid
      # Export the MCCI.fit_optimized function and any other necessary objects to each worker
      clusterExport(cl, varlist = c("MCCI.fit_optimized"))
      results <-
        foreach(alpha = alpha_grid, .combine = rbind) %:%
        foreach(lambda_2 = lambda_2_grid, .combine = rbind) %dopar% {
          test_error =IMR::error_metric$rmse
          lambda_1 = 0
          score = 0
          for (i in 1:n_folds) {
            data = fold_data[[i]]
            A_hat_test = MCCI.fit_optimized(data, lambda_1, lambda_2, alpha)
            # scores[i] = mean((data$A.test - A_hat_test)^2)
            # -- EDIT: Using MCCI's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$Y_valid)
            # counting the number of fits
            fit_counter = fit_counter + 1
          }
          score = score / n_folds
          c(alpha, lambda_2, score)
        }
      # Process results to find the best parameters
      # Edited on Dec 1st to pick the minimum score with highest lambda_2 value.
      min_score <- min(results[, 3])
      # Subset to only include results with the minimum score
      min_results <-
        results[results[, 3] == min_score, , drop = FALSE] # drop to keep it as df
      # In case of multiple results with the same score, find the one with the highest lambda_2
      if (nrow(min_results) > 1) {
        best_result <- min_results[which.max(min_results[, 2]),]
      } else {
        best_result <-
          min_results  # If only one row, it's already the best result
      }
      #best_result <- results[which.min(results[, 3]), ] # old line
      # Extract the best parameters
      best_params <-
        list(alpha = best_result[1],
             lambda_1 = 0,
             lambda_2 = best_result[2])
      best_score <- best_result[3]
      # close the cluster
      stopCluster(cl)
    }
    #--------------------------------------------
    # fixing optimal values of lambda 2 and alpha and optimizing for lambda 1 separately
    lambda_2 = best_params$lambda_2
    alpha = best_params$alpha
    best_score = Inf
    for (lambda_1 in lambda_1_grid) {
      score = 0
      for (i in 1:n_folds) {
        data = fold_data[[i]]
        # compute the estimates with a modified fit function
        A_hat_test = MCCI.fit_optimized(data, lambda_1, lambda_2, alpha)
        # -- EDIT: Using MCCI's formula in page 205 to compute the test error
        score = score + test_error(A_hat_test, data$Y_valid)
        # counting the number of fits
        fit_counter = fit_counter + 1
      }
      score = score / n_folds

      if (score < best_score) {
        best_score = score
        best_params$lambda_1 = lambda_1
      }
    }
    #---------------------------------------------------
    best_fit <- MCCI_fit(Y = Y,
                        X = X,
                        W = W,
                        lambda_1 = best_params$lambda_1,
                        lambda_2 = best_params$lambda_2,
                        alpha = best_params$alpha,
                        n1n2_optimized = n1n2_optimized,
                        return_rank = TRUE,
                        theta_estimator = theta_estimator)
    # counting the number of fits
    fit_counter = fit_counter + 1
    #---------------------------------------------------
    obj = list(
      best_parameters = best_params,
      best_score = best_score,
      fit=best_fit,
      total_num_fits = fit_counter,
      time = round(
        as.numeric(
          difftime(
            Sys.time(),
            start_time,
            units = "secs")
        )))
    obj$time_per_fit = obj$time / obj$total_num_fits

    return(obj)

  }
MC_Kfold_split <-
  function(n_rows,
           n_cols,
           n_folds,
           obs_mask,
           seed = NULL) {

    if(! is.null(seed)) set.seed(seed)
    # Create a data frame of all matrix indices
    indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)
    # we only consider non-missing data (ie, with obs_mask=1)
    indices <- indices[obs_mask == 1, ]
    # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
    indices <- indices[sample(1:nrow(indices)),]

    # Assign each observed index to one of k groups,
    # ensuring that the number of validation cells in each row are equal (or close)
    indices <- indices %>%
      mutate(row = as.numeric(row),
             col = as.numeric(col)) %>%
      group_by(col) %>%
      # the following is to shuffle within each row
      do(sample_n(., size = nrow(.))) %>%
      mutate(fold = rep(1:n_folds, length.out = n())) %>%
      ungroup()

    # Assign each index to one of k groups
    # Create a list to hold each fold
    folds <- vector("list", n_folds)
    for (i in 1:n_folds) {
      # Create a mask for the test cells in this fold
      #  1 -> train    (obs_mask=1) |
      #  0 -> valid    (obs_mask=1) |
      #  1 -> missing  (obs_mask=0) (missing)
      valid_mask <- matrix(1, nrow = n_rows, ncol = n_cols)
      test_indices <- indices[indices$fold == i,]
      valid_mask[obs_mask == 0] <- 1
      valid_mask[as.matrix(test_indices[, c("row", "col")])] <- 0
      # Store the mask
      folds[[i]] <- valid_mask
    }
    return(folds)
  }
