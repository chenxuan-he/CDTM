estep <- function(documents, update.mu, mod,
                  beta, lambda.old, mu, sigma, time_slices, step) {
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  Ts <- length(unique(time_slices))
  if(!update.mu) mu.i <- as.numeric(mu)
  
  if(mod=="DTM") sigma.ss <- diag(0, nrow=(K-1))
  else sigma.ss <- rep(list(diag(0, nrow=(K-1))), Ts)
  obs <- vector(mode="list", length=Ts)
  for(i in 1:Ts) {
    obs[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  if(mod=="CDTM"){
    sigmaentropy <- vector("list", length = Ts)
    siginv <- vector("list", length = Ts)
    for (tt in 1:Ts) {
      sigobj <- try(chol.default(sigma[[tt]]), silent=TRUE)
      if(inherits(sigobj,"try-error")) {
        sigmaentropy[[tt]] <- (.5*determinant(sigma[[tt]], logarithm=TRUE)$modulus[1])
        siginv[[tt]] <- solve(sigma[[t]])
      } else {
        sigmaentropy[[tt]] <- sum(log(diag(sigobj)))
        siginv[[tt]] <- chol2inv(sigobj)
      }
    }
  }else{
    sigobj <- try(chol.default(sigma), silent=TRUE)
    if(inherits(sigobj,"try-error")) {
      sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
      siginv <- solve(sigma)
    } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
    }
  }

  for(i in 1:N) {
    # update components
    if(min(time_slices) == 0) {time_slice <- time_slices[i]+1
    }else time_slice <- time_slices[i]
    doc <- documents[[i]]
    words <- doc[1,]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[time_slice]][,words,drop=FALSE]
    
    # core part of updating
    if (mod=="DTM") {
      doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                       doc=doc, sigmaentropy=sigmaentropy)
      # update sufficient statistics 
      sigma.ss <- sigma.ss + doc.results$eta$nu
      obs[[time_slice]][,words] <- doc.results$phis + obs[[time_slice]][,words]
      bound[i] <- doc.results$bound
      lambda[[i]] <- c(doc.results$eta$lambda)
      
    }else{
      doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv[[time_slice]], beta=beta.i, 
                                       doc=doc, sigmaentropy=sigmaentropy[[time_slice]])
      # update sufficient statistics 
      sigma.ss[[time_slice]] <- sigma.ss[[time_slice]] + doc.results$eta$nu
      obs[[time_slice]][,words] <- doc.results$phis + obs[[time_slice]][,words]
      bound[i] <- doc.results$bound
      lambda[[i]] <- c(doc.results$eta$lambda)
    }
  }
  
  lambda <- do.call(rbind, lambda)
  # NEW change: not update beta here
  return(list(sigma=sigma.ss, obs=obs, bound=bound, lambda=lambda))
  #else return(list(sigma=sigma.ss, bound=bound, lambda=lambda))
}

logisticnormalcpp <- function(eta, mu, siginv, beta, doc, sigmaentropy, 
                              method="BFGS", control=list(maxit=500),
                              hpbcpp=TRUE) {
  doc.ct <- doc[2,]
  Ndoc <- sum(doc.ct)
  #even at K=100, BFGS is faster than L-BFGS
  optim.out <- optim(par=eta, fn=lhoodcpp, gr=gradcpp,
                     method=method, control=control,
                     doc_ct=doc.ct, mu=mu,
                     siginv=siginv, beta=beta)
  
  if(!hpbcpp) return(list(eta=list(lambda=optim.out$par)))
  
  #Solve for Hessian/Phi/Bound returning the result
  hpbcpp(optim.out$par, doc_ct=doc.ct, mu=mu,
         siginv=siginv, beta=beta,
         sigmaentropy=sigmaentropy)
}

source("RcppExports.R")
sourceCpp("src/CDTMCfuns.cpp")
sourceCpp("src/RcppExports.cpp")
