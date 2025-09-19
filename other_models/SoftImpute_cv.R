simpute.cv <- function(y_train,
                       y_valid,
                       # mask_valid,
                       y_full,
                       n.lambda = 20,
                       lambda0_fun = softImpute::lambda0,
                       trace = FALSE,
                       print.best = TRUE,
                       tol = 5,
                       thresh = 1e-6,
                       rank.init = 10,
                       rank.limit = 50,
                       rank.step = 2,
                       maxit = 300,
                       test_error = IMR::error_metric$rmse,
                       seed = NULL) {
  # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
  start_time <- Sys.time()
  if(!is.null(seed))
    set.seed(seed)
  # valid_ind <- mask_valid == 0
  #y_full[y_full == 0] = NA
  #y_train[y_train == 0] = NA

  stopifnot(IMR::is.Incomplete(y_train))
  stopifnot(IMR::is.Incomplete(y_valid))
  stopifnot(IMR::is.Incomplete(y_full))

  lam0 <- lambda0_fun(y_full)

  y_train <- as.matrix(y_train)
  y_train[y_train == 0] = NA


  irow <- y_valid@i
  pcol <- y_valid@p

  lamseq <- seq(from = lam0,
                to = 0,
                length = n.lambda)


  rank.max <- rank.init
  warm <- NULL
  best_fit <-
    list(
      error = Inf,
      rank_M = NA,
      lambda = NA,
      rank.max = NA
    )
  counter <- 1
  loop_size <- 0

  for (i in seq(along = lamseq)) {
    loop_size <- loop_size + 1
    fiti <-
      softImpute::softImpute(
        as.matrix(y_train),
        type = "als",
        lambda = lamseq[i],
        rank.max = rank.max,
        warm.start = warm,
        thresh = thresh,
        maxit = maxit
      )

    # compute rank.max for next iteration
    rank <-
      sum(round(fiti$d, 4) > 0) # number of positive sing.values
    rank.max <- min(rank + rank.step, rank.limit)

    if(is.null(dim(fiti$u))) fiti$u <- as.matrix(fiti$u, ncol=1)
    vestim <- IMR:::partial_crossprod(fiti$u, fiti$d * t(fiti$v), irow, pcol)
    # soft_estim = fiti$u %*% (fiti$d * t(fiti$v))
    err = test_error(vestim, y_valid@x)
    # err = test_error(soft_estim[valid_ind], y_valid)
    #----------------------------
    warm <- fiti # warm start for next
    if (trace == TRUE)
      cat(
        sprintf(
          "%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
          i,
          lamseq[i],
          rank.max,
          rank,
          err
        )
      )
    #-------------------------
    # register best fir
    if (err <= best_fit$error) {
      best_fit$error = err
      best_fit$rank_M = rank
      best_fit$lambda = lamseq[i]
      best_fit$rank.max = rank.max
      counter = 1
    } else
      counter = counter + 1
    if (counter >= tol) {
      if (trace || print.best)
        cat(sprintf(
          "Performance didn't improve for the last %d iterations.",
          counter
        ))
      break
    }
  }
  #----------------------------------------
  if (print.best == TRUE)
    print(best_fit)
  #----------------------------------
  # one final fit on the whole data:
  y_full <- as.matrix(y_full)
  y_full[y_full == 0] = NA

  fiti <-
    softImpute::softImpute(
      y_full,
      type = "als",
      lambda = best_fit$lambda,
      rank.max = best_fit$rank.max,
      warm.start = warm,
      thresh = thresh,
      maxit = maxit
    )
  best_fit$fit <- fiti
  best_fit$total_num_fits = loop_size + 1
  # best_fit$estimates =  fiti$u %*% (fiti$d * t(fiti$v))
  # best_fit$rank_M    = qr(best_fit$estimates)$rank
  best_fit$rank = sum(fiti$d > 0)
  best_fit$time      = round(
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs")
    ))
  best_fit$time_per_fit = best_fit$time / best_fit$total_num_fits
  #---------------------------------------------------------------
  return(best_fit)
}
