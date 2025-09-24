library(MASS)
MA25.fit <- function(Y,
                     X,
                     max_iter = 30L,
                     tol=1e-6,
                     C_h=0.2,
                     delta_h=0.1,
                     r_bar = 8,
                     save_to_file=FALSE,
                     file_location = NULL){
  M <- MC_alt_LS_rank_solver(Y=Y, X=X, intercept_val=1,
                             max_iter = max_iter, .tol = tol, missing_model = 'logistic')
  rm(X, Y)
  M$.make_obs()
  # estimate missing model
  M$estimate_pi()
  #R <- summary(M$pi_est$fit)
  #Tab <- R$coefficients
  #rownames(Tab) <- c("(Intercept)", colnames(X))

  # Include only main effects for model fitting (best performing based on the article)
  M$set_fitting_cov(index = 1:(ncol(X)+1) )
  print(M$summary())
  # Rank estimation
  M$rank_estimation(C_h=C_h, delta_h=delta_h, penalty_names = c('h'), r_bar=r_bar)
  rhat <- M$rank_est$est['h'] # get the rank estimation result
  # Fitting with the estimated rank
  M$fitting(with_rs = rhat)
  if(save_to_file){
    if(is.null(file_location)){
      message("file location must be provided")
    }else{
      message(paste("Saving file at :", file_location))
      saveRDS(M, file_location)
    }
  }
  return(M)
}
# rank estimation for MC problem
#-------------------------------------------------------------
ER <- function(d = c(1,1), r.min = 1, r.max = length(d)-1){
  return(d[r.min:r.max]/d[r.min:r.max+1])
}

IC <- function(n = 2, m = 2, d = c(1,1), r.min = 1, r.max = length(d), pfun = NULL, ...){
  d.sq <- d^2
  if (is.null(pfun)){
    pfun <- function(n, m, alpha_n=1, delta_h = 0.2, C_h = 1, ...){
      u <- sqrt(alpha_n)*m*n/(m+n)
      pen <- c(outer(C_h, n^delta_h, '*'))/u
      return(pen)
    }
  }
  pen <- pfun(n,m, ...)
  ic <- log( sum(d.sq) - cumsum(d.sq[1:r.max])[r.min:r.max] ) - log(m*n) +  outer((r.min:r.max), pen, '*')
  rownames(ic)<- r.min:r.max
  return(ic)
}

eIC <- function(n = 2, m = 2, mse = 1, r = 1:length(mse), pfun = NULL, ...){
  if (is.null(pfun)){
    pfun <- function(n, m, alpha_n=1, obr=1, delta_h = 0.1, C_h = 0.9, ...){
      u <- alpha_n*m*n/(m+n)
      pen <- c(outer(C_h, n^delta_h, '*'))/sqrt(u)
      return(pen)
    }
  }
  pen <- pfun(n,m,...)

  eic <- log(mse) + outer(r,pen, '*')
  rownames(eic)<- r
  return(eic)
}
#-------------------------------------------------------------

# solving missing value problem with linear model
#-------------------------------------------------------------
lowrank.mc.iterative.svd <- function(Z, rank, miss = NULL, obs.index = NULL, pi.hat = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL,
                                     Z.ini = NULL, tol = 1e-6, max.iter = 50){

  # low rank approximation via iterative SVD
  find.missing <- function(){
    if (is.null(miss)) {
      reurn(Z.est[-obs.index$all])
    } else {
      return(Z.est[miss])
    }
  }

  n <- dim(Z)[1]
  m <- dim(Z)[2]

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat

  Z.est <- Z.ini
  err <- tol + 1
  iter <- 0
  ans <- list(Z = Z.est, L = NULL, R = NULL)
  while (err>tol & iter<max.iter){
    if(is.null(Z.est)) {
      Z.comp <- NULL
    } else {
      Z.comp <- find.missing()
    }
    ans <- svd.missing(Z = Z, rank = rank, miss = miss, obs.index = obs.index, miss.type = miss.type, pi.hat = pi.hat, Z.comp = Z.comp)
    Z.new <- ans$L %*% t(ans$R)
    err <- ifelse(is.null(Z.est), tol + 1, mean((Z.est - Z.new)^2))
    Z.est <- Z.new
    iter <- iter + 1
  }

  return(list(Z = Z.est, L = ans$L, R = ans$R, rank = rank,
              converge = list(iter = iter, reach.max = iter==max.iter, err = err)
  ))
}

lowrank.mc.iterative.reg <- function(Z, miss = NULL, obs.index = NULL, L.ini, R.ini = NULL,
                                     tol = 1e-6, max.iter = 50, epsl=1e-16){

  # low rank approximation via iterative regression

  L <- L.ini
  R <- R.ini
  Z.est <- if(is.null(R)) NULL  else L %*% t(R)
  iter <- 0
  err <- tol + 1
  while (iter<max.iter && err>tol){
    R <- reg.missing(Y = Z, X = L, miss = miss, obs.index_byj = obs.index$byj, epsl = epsl)
    L <- reg.missing(Y = t(Z), X = R, miss = if(is.null(miss)) NULL else t(miss), obs.index_byj = obs.index$byi, epsl = epsl)

    Z.new <- L %*% t(R)
    if (is.null(Z.est)){
      err <- tol + 1
    } else {
      err <- if (!is.null(miss)) mean((Z.est[!miss] - Z.new[!miss])^2) else mean((Z.est[obs.index$all] - Z.new[obs.index$all])^2)
    }
    Z.est <- Z.new
    iter <- iter + 1
  }
  return(list(Z = Z.est, L = L, R = R, rank = dim(L)[2],
              converge = list(iter = iter, reach.max = iter==max.iter, err = err)
  ))
}

pi.est <- function(miss = NULL, obs.index = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL, detail=F){

  if (!is.null(X)){
    vars <- apply(X, MARGIN = 2, var)
    X <- X[,vars!=0]
  }

  n <- if (!is.null(miss)) nrow(miss) else length(obs.index$byi)
  m <- if (!is.null(miss)) ncol(miss) else length(obs.index$byj)

  if (miss.type == 'MCR'){

    pi.hat <- if (!is.null(miss)) mean(!miss) else sum(lengths(obs.index$byi))/(n*m)
    return(list(pi.hat = pi.hat))

  } else if (miss.type == 'logistic'){
    if (is.null(X) || ncol(X)==0) stop('No X for logistic model')
    n.obs <- if (!is.null(miss)) rowSums(!miss) else lengths(obs.index$byi)
    n.obs <- cbind(n.obs,m-n.obs)

    fit <- glm(n.obs~1+X, binomial(link = "logit"))
    res <- list(pi.hat = fit$fitted.values, gam.hat = fit$coefficients)
    if (detail) res$fit <- fit
    return(res)
  } else {
    stop('Invalid missing_type')
  }
}

svd.missing <- function(Z, miss = NULL, obs.index = NULL, Z.comp = NULL, rank,
                        pi.hat = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL){

  n <- dim(Z)[1]

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat


  if (miss.type=='MCR' && pi.hat==1) {
    svd.result <- svd(Z, nu = rank, nv= rank)
    return(
      list( L = svd.result$u * sqrt(n),
            R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n),
            D = svd.result$d,
            rank = rank)
    )
  }

  # we pad the provided value, or we pad 0 and divide the matrix by pi-hat
  if (!is.null(Z.comp)){
    if (is.null(miss)) {
      Z[-obs.index$all] <- Z.comp
    } else { Z[miss] <-  Z.comp }
  } else {
    if (is.null(miss)){
      Z[-obs.index$all] <- 0
      Z <- Z/pi.hat
    } else {
      Z[miss] <- 0
      Z <- Z/pi.hat
    }
  }

  svd.result <- svd(Z, nu = rank, nv= rank)
  if (rank==0){
    return(
      list( D = svd.result$d,
            rank = rank)
    )
  }
  return(
    list( L = svd.result$u * sqrt(n),
          R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n),
          D = svd.result$d,
          rank = rank)
  )
}

reg.missing <- function(Y, X, miss = NULL, obs.index_byj = NULL, penalty = 0, epsl = 1e-16){

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  q <- if (!is.null(miss)) mean(!miss) else sum(lengths(obs.index_byj))/(n*m)

  if (q ==1) return( t(ginv(t(X) %*% X) %*% (t(X) %*% Y)) )

  p <- dim(X)[2]

  # solve each beta_j by using only those subject with non-missing response on j-th item

  B <- matrix(0, nr=m,nc=p)
  for (j in 1:m){
    if (is.null(miss)) {
      i<- obs.index_byj[[j]]
      k <- length(i)
    } else {
      i <- !miss[,j]
      k <- sum(i)
    }
    if (k==0) next
    Xr <- X[i, , drop = F]
    Yr <- Y[i,j, drop = F]
    H <- t(Xr) %*% Xr/k
    if (penalty>0) H <- H + diag(penalty, p)
    b <- as.vector( ginv(H) %*% (t(Xr) %*% Yr)/k )
    b[abs(b)<epsl] <- 0
    B[j,] <- b
  }

  return(B)
}

# The argument "miss" is the logic matrix that indicates the missing position
# obs.index is a list contains tow long lists obs.index$byi and obs.index$byj st obs.index$byi[[i]] or obs.index$byj[[j]] is a vector of observed index

# methods for MC problem
#-------------------------------------------------------------
MC.complete <- function(Y, X = NULL, rank = NULL, return.complete.hidden = T){
  ans <- list()
  n <- nrow(Y)
  if (!is.null(X)) {
    B <- ginv(t(X) %*% X, tol = 1e-6) %*% (t(X) %*% Y)
    Y <- Y - X %*% B
    ans$B <- t(B)
  }
  svd.result <- svd(Y, nu = rank, nv= rank)
  L = svd.result$u * sqrt(n)
  R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n)
  if (return.complete.hidden) ans$Z <- L %*% t(R)
  return( c(ans,
            list(L = L, R = R, D = svd.result$d, rank = rank, converge =list( iter = 0, reach.max = NA, err = NA) )
  )
  )
}



MC.covariate.iterative <- function(Y, X, miss = NULL, obs.index = NULL, rank = NULL,
                                   B.ini = NULL, L.ini = NULL, R.ini = NULL,
                                   max.iter = 50, tol = 1e-6, iter.method = c('reg', 'svd')[1],
                                   pi.hat = NULL, miss.type = c('MCR', 'logistic')[1]){

  # alternatively solving B and Z until converge
  # solving Z via regression or SVD

  # the given missing information could be a full matrix, 'miss', with (T,F) missing indicator
  # or a list of observed index, obs.index, including indices based on row ($byi), based on cloumn ($byj) and based on the matrix ($all)

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat

  if (miss.type=='MCR' && pi.hat==1) return(MC.complete(Y = Y, X = X, rank = rank))

  residul.matrix <- function(X=NULL,B=NULL,Z=NULL){
    W <- if (is.null(Z)) Y-X %*% t(B) else Y-Z
    if (!is.null(miss)){
      W[miss] <- 0
    } else {
      W[-obs.index$all] <- 0
    }
    if (inherits(Y, 'sparseMatrix')) W <- Matrix(W, sparse = TRUE)
    return(W)
  }

  if (is.null(B.ini) && (is.null(L.ini) || is.null(R.ini))) B.ini <- reg.missing(Y, X, miss,  obs.index_byj = obs.index$byj)

  ans <- list(L=L.ini, R=R.ini)
  B <- B.ini

  if (is.null(L.ini) || is.null(R.ini)){
    W <- residul.matrix(X = X , B = B.ini)
    if (iter.method=='reg' && (!is.null(L.ini) && !is.null(R.ini))){
      if (!is.null(L.ini)) ans$R <- reg.missing(W, L.ini, miss = miss, obs.index_byj = obs.index$byj)
      if (!is.null(R.ini)) ans$L <- reg.missing(t(W), R.ini, miss = if (is.null(miss)) NULL else t(miss), obs.index_byj = obs.index$byi)
    }
    if ( (is.null(L.ini) && is.null(R.ini)) || iter.method=='svd') {
      ans <- svd.missing(Z = W, miss = miss, obs.index = obs.index, rank = rank, pi.hat = pi.hat, miss.type = miss.type)
    }
  }

  if (is.null(B.ini)) B <- 0

  if (rank==0){
    ans$B <- B
    return(ans)
  }
  ans$Z <- ans$L %*% t(ans$R)
  Yh <- X%*%t(B) + ans$Z
  err.B <- NA
  err.Z <- NA
  err <- tol + 1
  iter <- 0
  while (err>tol && iter<max.iter) {
    W <- residul.matrix(Z = ans$Z)
    B <- reg.missing(Y = W, X = X, miss = miss, obs.index_byj = obs.index$byj)
    W <- residul.matrix(X = X, B = B)
    if (iter.method == 'reg') ans <- lowrank.mc.iterative.reg(Z = W, miss = miss, obs.index = obs.index, L.ini = ans$L, R.ini = ans$R, tol = tol, max.iter = 1)
    if (iter.method == 'svd') ans <- lowrank.mc.iterative.svd(Z = W, miss = miss, obs.index = obs.index, miss.type = miss.type, Z.ini = ans$Z, tol = tol, max.iter = 1, rank = rank)
    Yh.new <- X%*%t(B) + ans$Z
    err <- if(!is.null(miss)) mean((Yh[!miss] - Yh.new[!miss])^2) else mean((Yh[obs.index$all] - Yh.new[obs.index$all])^2)
    Yh <- Yh.new
    iter <- iter + 1
  }
  ans$B <- B
  ans$converge <- list( iter = iter, reach.max = (iter==max.iter), err = err)
  return(ans)
}

#-------------------------------------------------------------

# Hypothesis testing
#-------------------------------------------------------------
MC.bs.ht <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                     A.list, a.list, n.boots = 1000, alpha = 0.05, track = T){

  # Hypothesis testing on
  #
  #  H0: A*vec(t(B))=a
  #
  # w/ normal multiplier bootstrap method
  # Do the above test for given pairs A=A.list[i], a=a.list[i] i=1, 2, 3,...

  pivot <- function(A,a,B){ max(abs(A %*% matrix(t(B), ncol = 1) - a)) }

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  r <- dim(L.est)[2]
  p <- dim(B.est)[2]

  n.test <- length(A.list)

  test.names <- paste0('test_', 1:n.test)
  test.stat <- numeric(n.test)
  for (l in 1:n.test){
    test.stat[l] <- pivot(A.list[[l]], a.list[[l]], B.est)
  }
  bs.stat <- array(dim=c(n.test, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(n.test)
  names(pv) <- test.names
  names(test.stat) <- test.names

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))

  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)

  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est

  for (k in 1:n.boots){
    for (l in 1:n.test){
      tmp <- Sig.X.inv.Xt %*% (rnorm(n)*W)/n
      bs.stat[l,k] <- pivot(A.list[[l]], 0, t(tmp))
    }
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat > test.stat)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.list, a = a.list)
  )# 'FALSE' means reject null
}

MC.bs.ht.eachj <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                           A.list, a.list, n.boots = 1000, alpha = 0.05, track = T){

  # A: q-by-p matrix, a: q-by-m vector
  # Hypothesis testing on
  #
  #   H0_j: A*b_j=a
  #
  # for each b_j=t(B[j,]) w/ normal multiplier bootstrap method
  # Do the above test for given pairs A=A.list[i], a=a.list[i] i=1, 2, 3,...

  pivot <- function(A,a,B){ apply(abs(A %*% t(B) - c(a)), MARGIN = 2, max)  }

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  r <- dim(L.est)[2]
  p <- dim(B.est)[2]

  n.test <- length(A.list)
  test.names <- list(paste0('beta_',1:m), paste0('test_',1:n.test))
  test.stat <- array(dim=c(m, n.test), dimnames = test.names)
  for (l in 1:n.test){
    test.stat[,l] <- pivot(A.list[[l]], a.list[[l]], B.est)
  }
  pv <- array(dim=c(m, n.test), dimnames = test.names)
  bs.stat <- array(dim=c(m, n.test, n.boots), dimnames = c(test.names, list(1:n.boots)))

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))

  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)

  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est

  for (k in 1:n.boots){
    for (l in 1:n.test){
      tmp <- Sig.X.inv.Xt %*% (rnorm(n*m)*W)/n
      bs.stat[,l,k] <- pivot(A.list[[l]], 0, t(tmp))
    }
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- apply(bs.stat>c(test.stat), MARGIN = c(1,2), FUN = mean)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.list, a = a.list)
  )# 'FALSE' means reject null
}

MC.bs.ht.allj <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                          A.matrix, a.matrix = 0, n.boots = 1000, alpha = 0.05, track = T){

  # A.matrix: n.test x p
  # a.matrix: n.test x m
  #
  # H0: <v, b_j>=a_j for all j
  #
  # where <,> is the inner-product of vectors, v=A.matrix[k,], a_j=a[k,j] for k = 1, 2, ..., n.test.
  #
  # Best use for testing contrasts, combining groups.

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  r <- dim(L.est)[2]
  p <- dim(B.est)[2]

  n.test <- nrow(A.matrix)

  test.names <- paste0('test_', 1:n.test)
  test.stat <- apply( abs(A.matrix %*% t(B.est) - a.matrix), MARGIN = 1, max)
  bs.stat <- array(dim=c(n.test, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(n.test)
  names(pv) <- test.names
  names(test.stat) <- test.names

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  A.Sig.X.inv.Xt <- A.matrix %*% solve(t(pi.hat*X)%*%X/n, t(X))

  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)

  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est

  N <- prod(dim(A.Sig.X.inv.Xt))
  for (k in 1:n.boots){
    tmp <- abs((rnorm(N)*A.Sig.X.inv.Xt) %*% W/n)
    bs.stat[,k] <- apply(tmp, MARGIN = 1, max)
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat>test.stat)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.matrix, a = a.matrix)
  )# 'FALSE' means reject null
}

MC.bs.ht.simple <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                            n.boots = 1000, alpha = 0.05, track = T, selec=NULL){

  # For each p in the 'selec', test H0: B[, p]=0
  # Use MC.bs.ht() for more general tests

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  r <- dim(L.est)[2]
  p <- dim(B.est)[2]

  if(is.null(selec)) selec <- (1:p)
  selec <- sort(unique(selec[selec<=p]))
  p <- length(selec) # this is now the number of vectors of B in the tests


  test.names <- sprintf('beta_%d', selec)

  test.stat <- apply(abs(B.est[, selec, drop=F]), MARGIN = 2, max)
  bs.stat <- array(dim=c(p, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(p);
  names(pv) <- test.names
  names(test.stat) <- test.names

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))
  Sig.X.inv.Xt <- Sig.X.inv.Xt[selec,]

  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)

  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est

  for (k in 1:n.boots){
    tmp <- (Sig.X.inv.Xt*rnorm(p*n)) %*% W/n
    bs.stat[,k] <- apply(abs(tmp), MARGIN = 1, max)
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }

  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat>test.stat)

  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat)
  )# 'FALSE' means reject null
}

MC.bs.ht.all.zero <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                              n.boots = 1000, alpha = 0.05, track = T){

  # A simple testing on H0: B=0. This means all the covariates have no effect
  # Use MC.bs.ht() for more general tests

  n <- dim(Y)[1]
  m <- dim(Y)[2]

  r <- dim(L.est)[2]
  p <- dim(B.est)[2]

  test.stat <- max(abs(B.est))
  bs.stat <- numeric(n.boots)
  pv <- numeric(1)

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))

  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)

  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est

  pv <- 0
  for (k in 1:n.boots){
    tmp <- Sig.X.inv.Xt %*% (rnorm(n)*W/n)
    bs.stat[k] <- max(abs(tmp))
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- mean(bs.stat>test.stat)

  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat)
  )# 'FALSE' means reject null
}

#-------------------------------------------------------------

# Asymptotic covariance matrix
#-------------------------------------------------------------
MC.asymp.cov.cond.LFX <- function(X, Y, miss = NULL, obs.index = NULL, miss.type = c('MCR', 'logistic')[1],
                                  B.est, L.est, R.est, pi.hat=NULL, Zi = NULL, Zj = NULL, Bj = NULL, Bp = NULL){

  # input:
  # The data, the estimators and the desired index vectors: Zi, Zj, Bj, Bp. e.g. Zi=c(1,2,3), Zj=c(1,3,5), Bj=c(1:10), Bp=c(1:3)
  # Usually Bp is 1:p (or leave it NULL), so you'll get the full cov matrix of B_j

  # output:
  # the estimated variance of the estimator Z_ij for any i in Zi and j in Zj will be in out$Z.var
  # the estimated variance of the estimator EY_ij for any i in Zi and j in Zj will be in out$EY.var
  # the estimated covariance matrix of the estimator B_j[Bp] for any j in Bj will be in out$B.cov

  p <- dim(X)[2]
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  r <- dim(L.est)[2]

  # make sure each index is valid and unique
  Zi <- sort(unique(Zi[Zi<=n]))
  Zj <- sort(unique(Zj[Zj<=m]))
  Bj <- sort(unique(Bj[Bj<=m]))
  Bp <- sort(unique(Bp[Bp<=p]))

  # count how many we need, and get all if null
  zn <- ifelse(is.null(Zi), n, length(Zi))
  zm <- ifelse(is.null(Zj), m, length(Zj))
  bm <- ifelse(is.null(Bj), m, length(Bj))
  bp <- ifelse(is.null(Bp), p, length(Bp))

  if (bm*bp==0){
    B.cov <- NULL
    bm <- 0
  } else {
    B.cov <- array(0, dim = c(bp, bp, bm))
    p.names <- if (bp<p) Bp else 1:bp
    j.names <- if (bm<m) Bj else 1:bm
    dimnames(B.cov) <- list(p.names, p.names, j.names)
  }

  if (zn*zm==0){
    Z.var <- NULL
    EY.var <- NULL
    zm <- 0
  } else {
    Z.var <- matrix(0, nr=zn, nc=zm)
    i.names <- if (zn<n) Zi else 1:zn
    j.names <- if (zm<m) Zj else 1:zm
    dimnames(Z.var) <- list(i.names, j.names)
  }

  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat

  Z.est <- L.est%*%t(R.est)

  # error variance estimation
  if (!is.null(miss)){
    W <- Y - X %*% t(B.est)
    W[miss] <- Z.est[miss]
  } else {
    W <- Z.est
    U <- X %*% t(B.est)
    W[obs.index$all] <- Y[obs.index$all] - U[obs.index$all]
  }

  d <- svd(W, nu = 0, nv = 0)$d^2
  df <- m*n*mean(pi.hat) - m*p - (m+n-r)*r - 2*sum( d[(r+1):min(n,m)]/outer(-d[(r+1):min(n,m)], d[1:r], '+')  )
  sigE <- sum(d[(r+1):min(n,m)])/df
  rm(d, W, df)

  # constructing var and cov matrix
  sigX.inv <- ginv(t(X*pi.hat)%*%X/n)
  sigR.inv <- ginv(t(R.est) %*% R.est/m)
  sigL.inv <- ginv(t(L.est*pi.hat) %*% L.est/n)

  js <- if (is.null(Bj) || is.null(Zj)) 1:m else sort(unique(c(Bj, Zj)))

  zj <- 1
  bj <- 1
  for (j in js){
    tmp <-  sigX.inv %*% ( t(X) %*% (X * (pi.hat*Z.est[,j])^2 )) %*% sigX.inv/n #middle part of zeta_ij
    if (bj <= bm && (is.null(Bj) || Bj[bj]==j)) {
      B.cov[,,bj] <- if (is.null(Bp)) (sigE*sigX.inv + tmp)/n else (sigE*sigX.inv[Bp, Bp] + tmp[Bp, Bp])/n
      bj <- bj +1
    }
    if (zj <= zm && (is.null(Zj) || Zj[zj]==j)){
      Z.var[,zj] <- if (is.null(Zi)) rowSums( (X %*% tmp) * X )/n else rowSums( (X[Zi,,drop=F] %*% tmp) * X[Zi,,drop=F] )/n
      zj <- zj +1
    }
  }
  if (zn*zm==0)  return(list( B.cov = B.cov , Z.var = Z.var, EY.var = EY.var))
  if (!is.null(Zi)) {
    L.est <- L.est[Zi,,drop=F]
    X <- X[Zi,,drop=F]
    if (miss.type == 'MCR') {
      pi.hat <- rep(pi.hat, length(Zi))
    } else {
      pi.hat <- pi.hat[Zi]
    }
  }
  if (!is.null(Zj)) R.est <- R.est[Zj,,drop=F]
  tmp <- sigE*(outer( 1/pi.hat, rowSums((R.est %*% sigR.inv) * R.est )/m, '*') + rowSums((L.est %*% sigL.inv) * L.est )/n)
  Z.var <- Z.var + tmp
  EY.var <- tmp + sigE * rowSums( (X %*% sigX.inv) * X )/n
  dimnames(EY.var) <- dimnames(Z.var)
  return(list( B.cov = B.cov , Z.var = Z.var, EY.var = EY.var))
}

#-------------------------------------------------------------
#' @importClassesFrom Matrix Matrix
#' @importFrom methods setClassUnion setOldClass setRefClass
library(MASS)
library(Matrix)
library(methods)
library(magrittr)

# Update: 29 JAN 2024 ----- Niu, PoYao

setClassUnion("any_matrix", c("matrix", "Matrix"))

MC_model <- setRefClass('MC_model',
                        # ===== Fields =====
                        fields = list(
                          # main data -----
                          Y = 'any_matrix',
                          X = 'any_matrix',

                          # related to covariate X -----
                          .colmeans_X = 'vector',
                          .colsds_X = 'vector',
                          intercept = 'numeric',
                          p = 'integer',

                          # related to response Y -----
                          dim = 'vector',
                          .missing_type = 'character',
                          missing_model = 'character',
                          .obs = 'list',
                          pi_est = 'list',
                          obs_count = 'integer',

                          # fitting parameters -----
                          r = 'integer',
                          fitting_cov='integer',
                          pi_fitting_cov='integer',
                          max_iter = 'integer',
                          .tol = 'numeric',

                          # fitting results -----
                          ini_val = 'list',
                          fit = 'ANY',
                          rank_est = 'list',
                          .solved = 'ANY',

                          # user define -----
                          other_info = 'list'

                          # -----
                        ),

                        # ===== Methods =====
                        methods = list(
                          initialize = function(X=NULL, .missing_type = 'zero', missing_model = c('MCR', 'logistic')[1],
                                                intercept_in_X = F, intercept_val = ifelse(intercept_in_X>0,1,0),
                                                max_iter = 30L, .tol = 1e-6, ...){

                            .self$initFields(..., .missing_type = .missing_type, missing_model = missing_model, max_iter = max_iter, .tol = .tol)

                            .self$dim <- dim(.self$Y)

                            if(.self$.missing_type!='zero' && .self$.missing_type!='na') .self$.missing_type <- 'zero'

                            .self$set_X(new_X = X, intercept_in_X = intercept_in_X, intercept_val = intercept_val)
                          },

                          # setting, modifying -----
                          .make_obs = function(force = F){
                            if (length(.obs)!=0 && !force){
                              cat(sprintf('$.obs already exists ... (%s)\n', Sys.time()))
                            } else {
                              cat(sprintf('\nMaking obsevation index... (%s)\n', Sys.time()))
                              if (.missing_type == 'na') tmp <- !is.na(Y)
                              if (.missing_type == 'zero') tmp <- Y!=0

                              .obs$byi <<- tmp %>% apply(1, which)
                              .obs$byj <<- tmp %>% apply(2, which)
                              .obs$all <<- which(tmp)
                              cat(sprintf('$.obs done! ... (%s)\n', Sys.time()))
                            }


                            if (is.null(.obs$all)){
                              if (.missing_type == 'na') tmp <- !is.na(Y)
                              if (.missing_type == 'zero') tmp <- Y!=0
                              .obs$all <<- which(tmp)
                            }
                            obs_count <<- length(.obs$all) %>% as.integer()
                          },

                          estimate_pi = function(forced = F){
                            if (!is.null(pi_est$pi.hat) && !forced) {
                              cat(sprintf('$pi_est$pi.hat already exists ... (%s)\n', Sys.time()))
                              return(invisible(NULL))
                            }
                            .make_obs()
                            cat(sprintf('\nEstimating pi with %s model ... (%s)\n', missing_model, Sys.time()))
                            pi_est <<- pi.est(obs.index = .obs, miss.type = missing_model, X = X[,pi_fitting_cov,drop=F], detail = T)
                            cat(sprintf('$pi_est done! ... (%s)\n', Sys.time()))
                          },

                          set_X = function(new_X=NULL, intercept_in_X = F, intercept_val = ifelse(intercept_in_X>0, 1, 0)){
                            fitting_cov <<- integer(0)
                            pi_fitting_cov <<- integer(0)
                            if (is.null(new_X)) {
                              .self$X <- matrix(0, nr = .self$dim[1], nc = 0)
                            } else {
                              k <- as.integer(intercept_in_X)
                              if (k>0) .self$X <- new_X[,-k, drop=F] else .self$X <- new_X
                              if (is.null(colnames(.self$X))) colnames(.self$X) <- paste0('V',1:ncol(.self$X))
                            }
                            intercept <<- 0
                            .self$p <- dim(.self$X)[2]
                            if (intercept_val!= intercept) .self$set_intercept(value = intercept_val, clear_results = F)
                            .self$standardized_X <- F
                            .self$.colmeans_X <- numeric()
                            .self$.colsds_X <- numeric()
                            if (p>0){
                              W <- if(intercept==0) X else X[,-1, drop=F]
                              .colmeans_X <<- colMeans(W)
                              .colsds_X <<- apply(W, MARGIN = 2,FUN = sd)
                            }

                            set_fitting_cov()
                            set_pi_fitting_cov()
                          },

                          set_fitting_cov = function(index=NULL){
                            if (ncol(X)==0){
                              fitting_cov <<- integer(0)
                              p <<- 0L
                              clear_all_results(F)
                              return(invisible(NULL))
                            }
                            if (is.null(index)) index <- 1:ncol(X)

                            index <- sort(unique(as.integer(index[index<= ncol(X) & index>0])))
                            if (length(index)==length(fitting_cov) && all(index==fitting_cov)) {
                              cat('Same index. Nothing has changed.\n')
                              return(invisible(NULL))
                            }
                            fitting_cov <<- index
                            p <<- length(index)

                            clear_all_results(F)
                          },

                          set_pi_fitting_cov = function(index=NULL){
                            if (ncol(X)==0 || missing_model=='MCR'){
                              pi_fitting_cov <<- integer(0)
                              clear_pi_est()
                              return(invisible(NULL))
                            }
                            if (is.null(index)) index <- 1:ncol(X)

                            index <- sort(unique(as.integer(index[index<= ncol(X) & index>0])))
                            if (length(index)==length(pi_fitting_cov) && all(index==pi_fitting_cov)) {
                              cat('Same index. Nothing has changed.\n')
                              return(invisible(NULL))
                            }
                            pi_fitting_cov <<- index

                            clear_pi_est()
                          },

                          set_intercept  = function(value=1, clear_results = T){
                            if (intercept == value) {
                              cat('Same intercept value. Nothing has changed.\n')
                              return(invisible(NULL))
                            }
                            if (intercept != 0 && value!=0) {
                              X <<- cbind(value, X[,-1])
                              if (1 %in% fitting_cov){
                                if (length(fit)>0){
                                  for (k in 1:length(fit)){
                                    if (is.null(fit[[k]])) next
                                    fit[[k]]$B[,1] <<- fit[[k]]$B[,1]*intercept/value
                                  }
                                }
                                if (length(ini_val)>0){
                                  ini_val$B[,1] <<- ini_val$B[,1]*intercept/value
                                }
                              }
                            } else if (intercept==0) {
                              X <<- cbind(intercept=value, X)
                              fitting_cov <<- fitting_cov + 1L
                              pi_fitting_cov <<- pi_fitting_cov + 1L
                            }else if (value==0) {
                              X <<- X[,-1]
                              if (1 %in% fitting_cov){
                                clear_all_results(F)
                              }
                              fitting_cov <<- fitting_cov[fitting_cov>1] - 1L
                              pi_fitting_cov <<- pi_fitting_cov[pi_fitting_cov>1] - 1L
                            }
                            intercept <<- value
                          },

                          # organizing -----
                          clear_all_results = function(clear_pi_est=T){
                            clear_initial()
                            clear_rank_est()
                            clear_fit()
                            if (clear_pi_est) clear_pi_est()
                          },

                          clear_initial = function(){
                            .self$ini_val <- list()
                            cat('$ini_val: ')
                            print(.self$ini_val)
                          },

                          clear_fit = function(clear_index = NULL){
                            if (!is.null(clear_index)){
                              fit[clear_index] <<- NULL
                              keep_idx <- setdiff(1:length(r), clear_index)
                              .solved <<- .solved[keep_idx]
                              r <<- r[keep_idx]
                            } else {
                              fit <<- list()
                              .solved <<- logical()
                              r <<- integer()
                              cat('$fit: ')
                              print(.self$fit)
                            }
                          },

                          clear_rank_est = function(){
                            .self$rank_est <- list()
                            cat('$rank_est: ')
                            print(.self$rank_est)
                          },

                          clear_pi_est = function(){
                            .self$pi_est <- list()
                            cat('$pi_est: ')
                            print(.self$pi_est)
                          },

                          rearrange_fit = function(ix = NULL){
                            if (is.null(ix)) ix <- sort(r, index.return = T)$ix
                            r <<- r[ix]
                            fit <<- fit[ix]
                            .solved <<- .solved[ix]
                          }
                          # -----
                        )
)


# Basic class ----------
MC_alt_LS_solver <- setRefClass('MC_alt_LS_solver',
                                contains = 'MC_model',

                                # ===== Fields =====
                                fields = list(
                                  # related to covariate X -----
                                  standardized_X = 'logical',

                                  # related to response Y -----
                                  scaled_Y = 'logical',
                                  .multiplier_Y = 'numeric',
                                  .offset_Y = 'numeric'
                                  # -----
                                ),

                                # ===== Methods =====
                                methods = list(

                                  initialize = function(max_iter = 3L, .tol = 1e-6, ...){
                                    callSuper(max_iter = max_iter, .tol = .tol, ...)
                                    .self$scaled_Y <- F
                                    .self$.solved <- rep(F, length(.self$r))
                                  },

                                  # pre processing -----
                                  scale_Y = function(stdize = T, colwise = T, off = NULL, mult = NULL){
                                    if (scaled_Y) stop('Already scaled: .$scaled_Y = TRUE')

                                    cat(sprintf('\nScaling Y matrix... (%s)\n', Sys.time()))
                                    m <- .self$dim[2]
                                    if (!is.null(off)) .offset_Y <<- if (colwise) rep(off, len=m) else off[1]
                                    if (!is.null(mult)) .multiplier_Y <<- if (colwise) rep(mult, len = m) else mult[1]


                                    find_offset <- function(v){
                                      if (diff(range(v))==0) {
                                        return(0)
                                      } else {
                                        return( ifelse(stdize, mean(v), min(v)) )
                                      }
                                    }

                                    find_multiplier <- function(v){
                                      drv <- diff(range(v))
                                      if (drv==0) {
                                        return(v[1])
                                      } else {
                                        return( ifelse(stdize, sd(v), drv) )
                                      }
                                    }

                                    .make_obs()

                                    if (!colwise){
                                      if (is.null(off)) .offset_Y <<- find_offset(Y[.obs$all])
                                      if (is.null(mult)) .multiplier_Y <<- find_multiplier(Y[.obs$all])
                                      Y[.obs$all] <<- (Y[.obs$all]-.offset_Y)/.multiplier_Y
                                    } else {
                                      if (is.null(off)) .offset_Y <<- numeric(m)
                                      if (is.null(mult)) .multiplier_Y <<- numeric(m)
                                      for (j in 1:m){
                                        i <- .obs$byj[[j]]
                                        if (length(i)==0) next
                                        if (is.null(off)) .offset_Y[j] <<- find_offset(Y[i, j])
                                        if (is.null(mult)) .multiplier_Y[j] <<- find_multiplier(Y[i, j])
                                        Y[i, j] <<- (Y[i, j]-.offset_Y[j])/.multiplier_Y[j]
                                      }
                                    }
                                    clear_all_results()
                                    scaled_Y <<- T

                                  },

                                  standardize_X = function(na.rm = T){
                                    if (p==0) stop('No X matrix\n')
                                    if (standardized_X) stop('Already standardized: .$standardized_X = TRUE\n')
                                    cat(sprintf('\nStandardizing X matrix... (%s)\n', Sys.time()))
                                    W <- if(intercept==0) X else X[,-1]
                                    W <- t((t(W)-.colmeans_X)/.colsds_X)
                                    if (intercept != 0) {
                                      if (length(ini_val)>0){
                                        ini_val$B[,1] <<- ini_val$B %*% matrix(c(intercept, .colmeans_X), nc=1)
                                        ini_val$B[,-1] <<- ini_val$B[,-1] %*% diag(.colsds_X)
                                      }
                                      if (length(fit)>0){
                                        for (k in 1:length(fit)){
                                          fit[[k]]$B[,1] <<- fit[[k]]$B %*% matrix(c(intercept, .colmeans_X), nc=1)
                                          fit[[k]]$B[,-1] <<- fit[[k]]$B[,-1] %*% diag(.colsds_X)
                                        }
                                      }
                                      intercept <<- 1
                                      X <<- cbind(intercept,W)
                                    } else {
                                      X <<- W
                                      clear_all_results()
                                    }
                                    standardized_X <<- T
                                  },

                                  # fitting -----
                                  nllkh = function(L = NULL, R = NULL, B = NULL, Y_hat = NULL, with_r = NULL, limits = c(-Inf, Inf)){
                                    if (length(.obs)==0) .make_obs()

                                    if (is.null(Y_hat)) {
                                      if (is.null(L) || is.null(R) || (is.null(B) & p>0)){
                                        k <- which(r==with_r)
                                        L <- fit[[k]]$L
                                        R <- fit[[k]]$R
                                        if (p>0){
                                          B <- fit[[k]]$B
                                        }
                                      }
                                    }
                                    val <- 0
                                    for (j in 1:dim[2]){
                                      i <- .obs$byj[[j]]
                                      if (length(i)==0) next
                                      if (is.null(Y_hat)) {
                                        y_hat <- L[i,,drop=F] %*% t(R[j,,drop=F])
                                        if (p>0) y_hat <- y_hat + X[i,fitting_cov,drop=F] %*% t(B[j,,drop=F])
                                      } else {
                                        y_hat <- Y_hat[i,j]
                                      }
                                      y <- Y[i,j]

                                      off_set <- 0
                                      if (scaled_Y) {
                                        multiplier <- if (length(.multiplier_Y)==1 ) .multiplier_Y else .multiplier_Y[j]
                                        off_set <- if (length(.off_set_Y)==1 ) .offset_Y else .offset_Y[j]
                                        y_hat <- y_hat*multiplier
                                        y <- y*multiplier
                                      }

                                      if (!is.infinite(limits[1])) y_hat[(y_hat + off_set)<limits[1]] <- limits[1] - off_set
                                      if (!is.infinite(limits[2])) y_hat[(y_hat + off_set)>limits[2]] <- limits[2] - off_set

                                      val <- val + sum((y - y_hat)^2)
                                    }
                                    return(val)
                                  },

                                  .make_initial = function(with_r=NULL, tol=.tol){
                                    if (is.null(with_r)) with_r <- max(r)
                                    if (length(ini_val)!=0 && ini_val$rank >= with_r){
                                      cat(sprintf('$ini_val already exists with r = %d... (%s)\n', ini_val$rank, Sys.time()))
                                      return(invisible(NULL))
                                    }
                                    .make_obs()
                                    estimate_pi()
                                    cat(sprintf('\nInitialization with rank %d ... (%s)\n', with_r, Sys.time()))

                                    if (p>0) {
                                      ans <- MC.covariate.iterative(Y = Y, X = X[,fitting_cov,drop=F], obs.index = .obs, pi.hat = pi_est$pi.hat,
                                                                    max.iter = 0, rank = with_r, iter.method = 'reg', miss.type = missing_model, tol = tol)
                                    } else {
                                      ans <- svd.missing(Z = Y, obs.index = .obs, rank = with_r, miss.type = missing_model, pi.hat = pi_est$pi.hat)
                                    }
                                    # ans$D <- NULL
                                    ans$Z <- NULL
                                    ans$converge <- NULL
                                    ini_val <<- ans
                                    cat(sprintf('$ini_val done! ... (%s)\n', Sys.time()))
                                  },

                                  .solve = function(with_rs, max_iter = .self$max_iter, tol = .tol){

                                    for (with_r in with_rs){
                                      k <- which(r==with_r)
                                      if (.solved[k]) next
                                      cat(sprintf('r = %02d, max_iter = %02d ... (%s)\n', with_r, max_iter,  Sys.time()))

                                      if (p!=0){

                                        fit[[k]] <<- MC.covariate.iterative(Y=Y, X=X[,fitting_cov,drop=F], obs.index=.obs, rank=with_r, pi.hat = pi_est$pi.hat,
                                                                            B.ini=ini_val$B, L.ini=ini_val$L[,1:with_r, drop = F], R.ini=ini_val$R[,1:with_r, drop = F],
                                                                            max.iter = max_iter, iter.method = 'reg', miss.type = missing_model, tol = tol
                                        )
                                      } else {
                                        fit[[k]] <<- lowrank.mc.iterative.reg(Z=Y, obs.index = .obs,
                                                                              L.ini=ini_val$L[,1:with_r, drop = F], R.ini=ini_val$R[,1:with_r, drop = F],
                                                                              max.iter = max_iter, tol = tol)
                                      }

                                      fit[[k]]$Z <<- NULL
                                      fit[[k]]$rmse <<- sqrt(nllkh(with_r = with_r)/obs_count)
                                      .solved[k] <<- T
                                    }
                                  },

                                  fitting = function(with_rs = r, max_iter = .self$max_iter, tol = .tol, force = F){
                                    if (length(.solved)==0) {
                                      .solved <<- rep(F, length(r))
                                      fit <<- list()
                                    }
                                    with_rs %<>% unique %>% sort

                                    new_r <- setdiff(with_rs, r)
                                    old_r <- setdiff(with_rs, new_r)

                                    fit_r <- new_r
                                    if (length(old_r)>0){
                                      if (force){
                                        index <- outer(r, old_r, '==') %>% apply(MARGIN = 2, which)
                                        clear_fit(index)
                                        new_r <- sort(c(new_r, old_r))
                                        fit_r <- new_r
                                      } else {
                                        fit_r <- c(fit_r, setdiff(old_r, r[.solved]))
                                      }
                                    }

                                    r <<- c(r, new_r)
                                    .solved <<- c(.solved, rep(F, length(new_r)))
                                    if (length(fit_r)==0){
                                      cat('\nNo new fitting needed\n')
                                      return(invisible(NULL))
                                    }
                                    .make_initial(with_r = max(fit_r))
                                    cat('\nfitting with rank(s):', fit_r, '\n' )
                                    .solve(with_rs = fit_r, max_iter = max_iter, tol = tol)
                                  },

                                  # post analysis (hypothesis testing, asyptotic dist., prediction) -----
                                  ht_beta = function(A_list, a_list = rep(list(0), length(A_list)), n_boots = 500, alpha = 0.05, with_r, track = T){
                                    if (p==0) return(invisible(NULL))
                                    k <- which(r==with_r)
                                    H <- MC.bs.ht(X = X[,fitting_cov,drop=F], Y = as.matrix(Y), obs.index = .obs, miss.type = missing_model,
                                                  B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                  A.list = A_list, a.list = a_list, n.boots = n_boots, alpha = alpha, track = track)
                                    H$rank <- with_r
                                    return(H)
                                  },

                                  ht_beta_simple = function(ps = NULL, n_boots = 500, alpha = 0.05, with_r, track = T){
                                    if (p==0) return(invisible(NULL))

                                    k <- which(r==with_r)
                                    H <- MC.bs.ht.simple(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                                         B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                         alpha = alpha, n.boots = n_boots, track = track, selec = ps)
                                    if (is.null(ps)) ps <- 1:p
                                    H$describe <- sprintf('beta_j%d = 0, for all j', ps-(intercept!=0))
                                    H$rank <- with_r
                                    return(H)
                                  },

                                  ht_beta_allj = function(A.matrix = matrix(c(-1,1, rep(0, p - 2 + (intercept!=0))), nr=1), a.matrix = 0, n_boots = 500, alpha = 0.05, with_r, track = T){

                                    if (p==0) return(invisible(NULL))
                                    if (ncol(A.matrix)!=p) stop(sprintf('Incorrect column of matrix A , should be %d', p))

                                    pos <- A.matrix!=0
                                    describ <- character()
                                    for (l in 1:nrow(A.matrix)){
                                      u <- which(pos[l,])
                                      v <- A.matrix[l,u]
                                      if (length(v)==1) coeffs <- paste0(c('-', '')[(v>0) + 1], sprintf('%.2g', abs(v)))
                                      if (length(v)>1) coeffs <- paste0(c(c('-', '')[(v[1]>0) + 1], c(' - ', ' + ')[(v[-1]>0) + 1]), sprintf('%.2g', abs(v)))
                                      if(is.matrix(a.matrix)) {
                                        describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = a[%d,j]', l), ' for all j'))
                                      } else {
                                        if(length(a.matrix)==1) describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = %.4f', a.matrix), ' for all j'))
                                        if(length(a.matrix)==nrow(A.matrix)) describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = %.4f', a.matrix[l]), ' for all j'))
                                      }
                                    }

                                    k <- which(r==with_r)
                                    H <- MC.bs.ht.allj(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                                       B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                       A.matrix = A.matrix, a.matrix = a.matrix, alpha = alpha, n.boots = n_boots, track = track)
                                    H$describe <- describ
                                    H$rank <- with_r
                                    return(H)
                                  },

                                  ht_beta_eachj = function(A_list, a_list = rep(list(0), length(A_list)), n_boots = 500, alpha = 0.05, with_r, track = T){
                                    if (p==0) return(invisible(NULL))
                                    k <- which(r==with_r)
                                    H <- MC.bs.ht.eachj(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                                        B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                        A.list = A_list, a.list = a_list, n.boots = n_boots, alpha = alpha, track = track)
                                    H$describe <- 'A*beta_j = a'
                                    H$rank <- with_r
                                    return(H)
                                  },

                                  ht_beta_zero = function(n_boots = 500, alpha = 0.05, with_r, track = T){
                                    if (p==0) return(invisible(NULL))
                                    k <- which(r==with_r)
                                    H <- MC.bs.ht.all.zero(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                                           B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                           n.boots = n_boots, alpha = alpha, track = track)
                                    H$describe <- 'beta_jp = 0 for all j,p'
                                    H$rank <- with_r
                                    return(H)
                                  },

                                  asymp_var = function(Zi = NULL, Zj = NULL, Bj = NULL, Bp = NULL, with_r){
                                    k <- which(r==with_r)
                                    var <- MC.asymp.cov.cond.LFX(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                                                 B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                                                 Zi = Zi, Zj = Zj, Bj = Bj, Bp = Bp)
                                    return(var)
                                  },

                                  predict = function(i = 1:.self$dim[1], j = 1:.self$dim[2], with_r = NULL, cross=T, limits = c(-Inf, Inf)){
                                    if (is.null(with_r)){
                                      print('please provide the rank. (use \'$r\' to find available ranks)')
                                      return(invisible(NULL))
                                    } else {
                                      k <- which(r==with_r)
                                    }

                                    if (cross){
                                      i <- unique(i) %>% sort
                                      j <- unique(j) %>% sort
                                    }

                                    if (p>0) {
                                      Xi <- X[i,fitting_cov,drop=F]
                                      Bj <- fit[[k]]$B[j,, drop = F]
                                    }


                                    Li <- fit[[k]]$L[i,, drop = F]
                                    Rj <- fit[[k]]$R[j,, drop = F]

                                    if (cross){
                                      ans <- Rj %*% t(Li)
                                      if (p>0) ans <- ans + Bj %*% t(Xi)
                                      if (scaled_Y)  ans <- ans*.multiplier_Y[j] + .offset_Y[j]
                                      ans <- t(ans)
                                      rownames(ans) <- i
                                      colnames(ans) <- j
                                    } else {
                                      ans <- if (with_r==1) Li*Rj else rowSums(Li*Rj)
                                      if (p>0) ans <- if(p==1) ans + Xi*Bj else ans + rowSums(Xi*Bj)
                                      if (scaled_Y){
                                        ans <- ans*.multiplier_Y[j] + .offset_Y[j]
                                      }
                                      names(ans) <- paste0( i, ',' , j)
                                    }
                                    if (!is.infinite(limits[1])) ans[ans<limits[1]] <- limits[1]
                                    if (!is.infinite(limits[2])) ans[ans>limits[2]] <- limits[2]
                                    return( ans)
                                  },

                                  # utilities -----
                                  summary = function(detail = F){
                                    cat('\nY')
                                    cat(sprintf('\n  %-16s: %d, %d', 'dim', .self$dim[1], .self$dim[2]))
                                    cat(sprintf('\n  %-16s: %s', 'scaled', ifelse(scaled_Y, 'Yes', 'No')))

                                    cat('\n\nX')
                                    if (ncol(X)==0) {
                                      cat(sprintf('\n  %-16s\n', 'No X matrix'))
                                      return(invisible(NULL))
                                    }
                                    cat(sprintf('\n  %-16s: %d, %d', 'dim', nrow(X), ncol(X)))
                                    cat(sprintf('\n  %-16s: %s', 'standardized', ifelse(standardized_X, 'Yes', 'No')))
                                    cat(sprintf('\n  %-16s: %s', 'intercept', ifelse(intercept!=0, round(intercept, 2), 'None')))
                                    cat(sprintf('\n  %-16s:', 'satats'))
                                    A <- apply(X, MARGIN = 2, quantile)
                                    rownames(A) <- c('Min', '1st', 'Med', '3rd', 'Max')
                                    cat(sprintf('\n  %-8s  %s', '', paste0(sprintf('%10s', colnames(A)), collapse = ' ')))
                                    for (i in 1:5) cat(sprintf('\n  %-8s  %s', rownames(A)[i], paste0(sprintf('%10.3f', A[i,]), collapse = ' ')))
                                    cat(sprintf('\n  %-8s  %s', 'Mean', paste(c(sprintf('%10s','')[intercept>0],sprintf('%10.3f', .colmeans_X)), collapse = ' ')))
                                    cat(sprintf('\n  %-8s  %s', 'SD', paste(c(sprintf('%10s','')[intercept>0],sprintf('%10.3f', .colsds_X)), collapse = ' ')))

                                    #cat(sprintf('\n  %-16s:\n', 'peek')); print(X[1:min(5, .self$dim[1]),,drop=F])

                                    cat('\n\nModel')
                                    cat(sprintf('\n  %-16s: %s', 'fitting with', paste0(colnames(X)[fitting_cov], collapse = ', ')))
                                    cat(sprintf('\n  %-16s: %s', 'missing model', missing_model))
                                    if (missing_model=='logistic') cat(sprintf('\n  %-16s: %s', 'pi fitting with', paste0(colnames(X)[pi_fitting_cov], collapse = ', ')))
                                    cat('\n')
                                  }
                                  # -----
                                )
)

MC_alt_LS_rank_solver <- setRefClass('MC_alt_LS_rank_solver',
                                     contains = 'MC_alt_LS_solver',

                                     # ===== Methods =====
                                     methods = list(
                                       initialize = function(...){
                                         callSuper(...)
                                       },
                                       # rank estimation -----
                                       .get_r_bar = function(C = 1.5, ic_penalty_fun = NULL, tol=.tol, ...){
                                         estimate_pi()
                                         cat(sprintf('\ngetting r_bar for rank estimation ... (%s)\n', Sys.time()))
                                         if (p>0) {
                                           ans <- MC.covariate.iterative(Y = Y, X = X[,fitting_cov,drop=F], obs.index = .obs, pi.hat = pi_est$pi.hat,
                                                                         max.iter = 0, rank = 0, iter.method = 'reg', miss.type = missing_model, tol = tol)
                                         } else {
                                           ans <- svd.missing(Z = Y, obs.index = .obs, rank = 0, miss.type = missing_model, pi.hat = pi_est$pi.hat)
                                         }

                                         if (missing_model == 'MCR'){
                                           alp <- .self$pi_est$pi.hat
                                         }

                                         if (missing_model == 'logistic'){
                                           alp <- exp(.self$pi_est$gam.hat['(Intercept)'])
                                         }

                                         ic <- IC(n = .self$dim[1], m = .self$dim[2], d = ans$D, pfun = ic_penalty_fun, alpha_n=alp, ...)
                                         d_ic <- diff(ic)
                                         rank_est$r_bar <<- max(as.integer(C*(match(T, d_ic>0))), 5)
                                         cat(sprintf('$rank_est$r_bar = %d ... (%s)\n', rank_est$r_bar, Sys.time()))
                                       },

                                       .solve_rank = function(r_bar, max_iter = .self$max_iter, tol = .tol){

                                         if (p>0){
                                           W <- Y - X[,fitting_cov,drop=F] %*% t(ini_val$B)
                                           W[-.obs$all] <- 0
                                           W <- Matrix(W, sparse = TRUE)
                                         }
                                         for (with_r in 1:r_bar){
                                           if (!is.na(rank_est$mse[with_r])) next
                                           cat(sprintf('r = %02d, max_iter = %02d ... (%s)\n', with_r, max_iter,  Sys.time()))
                                           if (p>0){
                                             ans <- lowrank.mc.iterative.reg(Z=W, obs.index = .obs,
                                                                             L.ini=ini_val$L[,1:with_r, drop = F],
                                                                             max.iter = max_iter, tol = tol)
                                           } else {
                                             ans <- lowrank.mc.iterative.reg(Z=Y, obs.index = .obs,
                                                                             L.ini=ini_val$L[,1:with_r, drop = F],
                                                                             max.iter = max_iter, tol = tol)

                                           }
                                           rank_est$mse[with_r] <<- nllkh(L = ans$L, R = ans$R, B = ini_val$B)/obs_count
                                         }
                                       },

                                       rank_estimation = function(r_bar = NULL, C_rbar = 1.5, penalty_fun = NULL, penalty_names = NULL, ...){

                                         if (is.null(r_bar) || r_bar<1){
                                           if (is.null(rank_est$r_bar)) .get_r_bar(C = C_rbar)
                                         } else {
                                           rank_est$r_bar <<- as.integer(r_bar)
                                         }

                                         len_pen_old <- 0
                                         r_bar_old <- 0
                                         if (!is.null(rank_est$eIC)){
                                           len_pen_old <- ncol(rank_est$eIC)
                                           r_bar_old <- nrow(rank_est$eIC)
                                           rank_est$r_bar <<- max(rank_est$r_bar, r_bar_old)
                                           if (rank_est$r_bar>r_bar_old) {
                                             rank_est$eIC <<- rbind(rank_est$eIC, matrix(nr=rank_est$r_bar-r_bar_old, nc=len_pen_old))
                                           }
                                         }

                                         rank_est$mse <<- c(rank_est$mse, rep(NA, rank_est$r_bar - length(rank_est$mse)))

                                         len_pen <- length(eIC(n = 1, m = 1, mse = 1, r = 1, pfun = penalty_fun, ...))
                                         n_old <- 0
                                         for (name in penalty_names){
                                           if (name %in% colnames(rank_est$eIC)) n_old <- n_old+1
                                         }

                                         if (length(penalty_names)<len_pen){
                                           ext <- 1:(len_pen - length(penalty_names)) + length(penalty_names)- n_old + len_pen_old
                                           penalty_names <- c(penalty_names, paste0('pen[', ext, ']'))
                                         }
                                         penalty_names <- penalty_names[1:len_pen]

                                         .make_initial(with_r = rank_est$r_bar)
                                         cat(sprintf('\nRank estimation with max rank = %d ... (%s)\n', rank_est$r_bar, Sys.time()))
                                         .solve_rank(r_bar=rank_est$r_bar)

                                         if (missing_model == 'MCR'){
                                           alp <- .self$pi_est$pi.hat
                                         }

                                         if (missing_model == 'logistic'){
                                           alp <- exp(.self$pi_est$gam.hat['(Intercept)'])
                                         }

                                         cat(sprintf('\nCalculating eIC with penalty %s ... (%s)\n', paste(penalty_names, collapse = ', '), Sys.time()))
                                         eic <- eIC(n = .self$dim[1], m = .self$dim[2], mse = rank_est$mse, pfun = penalty_fun, alpha_n=alp, ...)
                                         colnames(eic) <- penalty_names

                                         for (name in penalty_names){
                                           if (name %in% colnames(rank_est$eIC)){
                                             rank_est$eIC[,name] <<- eic[,name, drop=F]
                                           } else {
                                             rank_est$eIC <<- cbind(rank_est$eIC,eic[,name, drop=F])
                                           }
                                         }
                                         row.names(rank_est$eIC) <<- paste0('r=', 1:rank_est$r_bar)
                                         rank_est$est <<- apply(rank_est$eIC, MARGIN = 2, which.min)
                                         cat(sprintf('$rank_est with max rank = %d done! ... (%s)\n', rank_est$r_bar, Sys.time()))
                                       }
                                       #-----
                                     )
)
