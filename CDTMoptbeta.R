source("CDTMoptbetaFuns.R")
#Optimizing beta
opt.beta <- function(beta.ss, time_slices, settings) {
  K <- settings$dim$K
  Ts <- settings$dim$Ts
  V <- settings$dim$V
  verbose <- settings$verbose
  
  # contain it as the sufficient statistics
  ## and normalize one as obs
  for (i in 1:Ts) {
    # to normalize
    obs.global_byT[[i]] <<- beta.ss[[i]]/rowSums(beta.ss[[i]])
  }
  
  # re-arrange obs
  for (k in 1:K) {
    tmp <- NULL
    for (t in 1:Ts) {
      tmp <- rbind(tmp, obs.global_byT[[t]][k,])
    }
    obs.global_byK[[k]] <<- tmp
  }
  for (k in 1:K) {
    # fit the sequential model topic-wise
    beta.k <- NULL
    # get the beta.k(matrix,Ts\times V)
    for (i in 1:Ts) {
      beta.k <- rbind(beta.k, beta.ss[[i]][k,])
    }
    # before we update beta, we need to update zeta and mean
    for (v in 1:V) {
      update_post_mean(v, k, beta.k, settings)
    }
    
    update_zeta(k, settings)
    ## core function to update
    # attention: the dimension of the obs/beta are not the same
    if(verbose) print(paste0("M-step...update obs...topic_",as.character(k)))
    update_beta_k(k, beta.k,  time_slices, settings)
  }
  # need to re-arrange obs to beta.ss and return
  for (t in 1:Ts) {
    tmp <- NULL
    for (k in 1:K) {
      tmp <- rbind(tmp, obs.global_byK[[k]][t,])
    }
    obs.global_byT[[t]] <<- tmp
  }
  
  beta.ss <- obs.global_byT
  logbeta <- NULL
  beta <- NULL
  beta_byT <- NULL
  for (k in 1:K) {
    logbeta[[k]] <- matrix(rep(0, Ts*V), nrow=V)
    beta[[k]] <- matrix(rep(0, Ts*V), nrow=V)
    for (v in 1:V) {
      for (t in 1:Ts) {
        logbeta[[k]][v,t] <- mean.global$mean[[k]][v,t+1]-log(zeta.global[[k]][t])
      }
    }
    for (t in 1:Ts) {
      beta[[k]][,t] <- exp(logbeta[[k]][,t])/sum(exp(logbeta[[k]][,t]))
    }
  }
  for (t in 1:Ts) {
    tmp <- NULL
    for (k in 1:K) {
      tmp <- cbind(tmp, beta[[k]][,t])
    }
    beta_byT[[t]] <- t(tmp)
  }
  return(list(beta_byK=beta, beta=beta_byT, logbeta=logbeta))
}

update_beta_k <- function(k, beta.k,  time_slices, settings) {
  verbose <- settings$verbose
  V <- settings$dim$V
  K <- settings$dim$K
  N <- settings$dim$N
  Ts <- settings$dim$Ts
  bound <- 0
  old_bound <- 0
  sslm_fit_threshold <- 1e-6
  sslm_max_iter <- 2
  converged <- sslm_fit_threshold + 1
  iter_ = 0
  while (converged > sslm_fit_threshold & iter_ < sslm_max_iter) {
    iter_ <- iter_+1
    old_bound <- bound
    # the last place to update
    if(verbose) print(paste0("M-step...update obs...topic_",as.character(k),
                             ", times_", as.character(iter_)))
    update_obs(k, beta.k,  settings)
    
    bound <- compute_bound(k, beta.k, time_slices, settings)
    converged = abs((bound - old_bound) / old_bound)
  }
}

