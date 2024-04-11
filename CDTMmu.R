opt.mu <- function(gamma_prior, lambda, covar=NULL,
                   maxits=1000, time_slices, mod) {
  if(mod == "DTM"){
    mu <- matrix(colMeans(lambda), ncol = 1)
    return(list(mu = mu, gamma = NULL))
  }
  if(min(time_slices) == 0) time_slices <- time_slices+1
  # Variational Linear Regression with a Gamma hyperprior
  ### the part we use in cdtm now
  Ts <- length(unique(time_slices))
  gamma <- rep(list(matrix(0, ncol = ncol(lambda), nrow = dim(covar)[2])), Ts)
  gamma_tmean <- matrix(0, ncol = ncol(lambda), nrow = dim(covar)[2])
  gamma.se <- rep(list(matrix(0, ncol = ncol(lambda), nrow = dim(covar)[2])), Ts)
  #gamma <- vector(mode="list",length=ncol(lambda))
  #gamma.se <- vector(mode="list",length=ncol(lambda))
  Xcorr <- crossprod(covar)
  sigk <- rep(1, K-1)

  if(is.null(gamma_prior)){
    gamma_prior <- rep(0, K-1)
  }
  for (i in 1:ncol(lambda)) {
    for (tt in 1:Ts) {
      indtt <- which(time_slices==tt)
      tmp <- vb.variational.reg(Y=lambda[indtt,i], X=covar[indtt,], Xcorr=Xcorr, maxits=maxits,
                                prior = gamma_prior[i]) 
      gamma[[tt]][,i] <- tmp$w
      gamma.se[[tt]][,i] <- tmp$se
    }
  }
  gamma_tmp <- matrix(0, ncol = ncol(lambda), nrow = dim(covar)[2])
  for (tt in 1:Ts) {
    gamma_tmp <- gamma_tmp+gamma[[tt]]/Ts
  }
  mu <- matrix(0, ncol = nrow(lambda), nrow = K-1)
  for (i in 1:nrow(lambda)) {
    mu[,i] <- t(matrix(covar[i,], nrow = 1) %*% gamma[[time_slices[i]]])
  }
  #if its not a regular matrix,coerce it as it won't be sparse.
  if(!is.matrix(mu)) {
    mu <- as.matrix(mu)
  }
  return(list(mu=mu, gamma=gamma, gamma.se=gamma.se, gamma_prior=gamma_tmp))
}

# Variational Linear Regression with a Half-Cauchy hyperprior 
# (Implementation based off the various LMM examples from Matt Wand)
# This code is intended to be passed a Matrix object
vb.variational.reg <- function(Y,X, b0=1, d0=1, Xcorr=NULL, maxits=1000, prior) {
  if(is.null(Xcorr)) Xcorr <- crossprod(X)
  XYcorr <- crossprod(X,Y)
  
  an <- (1 + nrow(X))/2
  D <- ncol(X)
  N <- nrow(X)
  w <- rep(0, ncol(X))
  error.prec <- 1 #expectation of the error precision
  converge <- 1000
  cn <- ncol(X) # - 1 for the intercept and +1 in the update cancel
  dn <- 1
  Ea <- cn/dn #expectation of the precision on the weights
  ba <- 1
  
  ct <- 1
  while(converge>.0001) {
    w.old <- w
    
    #add the coefficient prior.  Form depends on whether X is a Matrix object or a regular matrix.
    if(is.matrix(X)) {
      ppmat <- diag(x=c(0, rep(as.numeric(Ea), (D-1))),nrow=D) 
    } else {
      ppmat <- Diagonal(n=D, x=c(0, rep(as.numeric(Ea), (D-1))))
    }
    invV <- error.prec*Xcorr + ppmat
    #if its a plain matrix its faster to use the cholesky, otherwise just use solve
    if(is.matrix(invV)) {
      V <- try(chol2inv(chol(invV)), silent = TRUE)
      if(inherits(V, "try-error")){
        V <- solve(invV)
      }
    } else {
      #Matrix package makes this faster even when its non-sparse
      V <- solve(invV)     
    }
    w <- error.prec*V%*%(XYcorr+prior/error.prec)
    
    # parameters of noise model (an remains constant)
    sse <- sum((X %*% w - Y)^ 2)
    bn <- .5*(sse + sum(diag(Xcorr%*%V))) + ba
    error.prec <- an/bn
    ba <- 1/(error.prec + b0)
    
    #subtract off the intercept while working out the hyperparameters
    # for the coefficients
    w0 <- w[1]
    w <- w[-1]
    da <- 2/(Ea + d0)
    dn <- 2*da + (crossprod(w) + sum(diag(V)[-1]))
    Ea <- cn / dn
    #now combine the intercept back in 
    w <- c(w0,w)
    ct <- ct + 1
    if(ct > maxits) {
      stop("Prevalence regression failing to converge within iteration limit.  May want to try gamma.prior='L1'. You can change max iterations using control.  See stm documentation")
    }
    converge <- sum(abs(w-w.old))
  }
  sig.hat <- sse/N
  se.beta <- diag(sig.hat*solve(Xcorr))
  return(list(w=w,se=se.beta))
}


