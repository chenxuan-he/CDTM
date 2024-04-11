update_obs <- function(k, beta.k, settings) {
  OBS_NORM_CUTOFF <- 2
  STEP_SIZE <- 0.001
  TOL <- 1e-3
  V <- settings$dim$V
  Ts <- settings$dim$Ts
  runs <- 0
  verbose <- settings$verbose
  
  beta.new <- matrix(rep(0, Ts*V), nrow = Ts)
  totals <- apply(beta.k, 1, sum)
  # a matrix to get the derivative of the mean/beta
  mean_deriv_mtx <- matrix(rep(0, Ts*(Ts+1)), ncol = Ts+1)
  norm_cutoff_obs <- NULL
  
  totals <- apply(beta.k, 1, sum)
  tmpObs <- matrix(rep(0,Ts*V),nrow = Ts, ncol = V)
  for (v in 1:V) {
    w_counts <- beta.k[,v]
    counts_norm <- sqrt(sum(w_counts^2))
    if (counts_norm < OBS_NORM_CUTOFF & !is.null(norm_cutoff_obs)) {
      obs <- obs.global_byK[[k]][,v]
      norm_cutoff_obs <- obs
    }else {
      if (counts_norm<OBS_NORM_CUTOFF) {
        w_counts <- rep(0, length(w_counts))
      }
      for (t in 1:Ts) {
        mean_deriv_mtx[t,] <- compute_mean_deriv(k,v,t,mean_deriv_mtx[t,], settings)
      }
      deriv <- rep(0, Ts)
      tryCatch({
        tmp <- Rcgmin(fn = f_obs, gr=df_obs, par = obs.global_byK[[k]][,v], 
                      #upper = rep(1,Ts), lower = rep(0,Ts),
                      w_counts=w_counts, totals=totals, settings=settings,
                      mean_deriv_mtx=mean_deriv_mtx, v=v, deriv=deriv, beta.k=beta.k, k=k)
        tmpObs[,v] <- tmp$par
      }, error = function(e){
        print("An error occured in function Rcgmin")
      })
    }
  }
  # # softmax operation to obs.global_byK
  # for (tt in 1:Ts) {
  #   for (vv in 1:V) {
  #     obs.global_byK[[k]][tt,vv] <<- exp(tmpObs[tt,vv])/sum(exp(tmpObs[tt,]))
  #   }
  # }
}

f_obs <- function(x, w_counts, totals, settings, mean_deriv_mtx, v, deriv, beta.k, k){
  obs.global_byK[[k]][,v] <<- x
  obs_variance <- settings$beta$obs_variance
  chain_variance <- settings$beta$chain_variance
  zeta <- zeta.global[[k]]
  init_mult <- 1000
  Ts <- settings$dim$Ts
  val <- 0
  term1 <- 0
  term2 <- 0
  zeta <- zeta.global[[k]]
  update_post_mean(v, k, beta.k, settings)
  mean <- mean.global$mean[[k]][v,]
  variance <- var.global$var[[k]][v,]
  for (t in 2:(Ts+1)) {
    mean_t <- mean[t]
    mean_t_prev <- mean[t-1]
    
    val <- mean_t-mean_t_prev
    term1 = term1+val^2
    term2 = term2+w_counts[t-1]*mean_t-totals[t-1]*exp(mean_t+variance[t]/2)/zeta[t-1]
  }
  if(chain_variance > 0){
    term1 <- -(term1/(2*chain_variance))
    term1 <- term1-mean[1]*mean[1]/(2*init_mult*chain_variance)
  }
  else temr1 <- 0
  final <- -(term1+term2)
  return(final)
}

df_obs <- function(x, w_counts, totals, settings, mean_deriv_mtx, v, deriv, beta.k, k) {
  obs.global_byK[[k]][,v] <<- x
  update_post_mean(v, k, beta.k, settings)
  deriv <- compute_obs_deriv(v, k, w_counts, totals, mean_deriv_mtx, deriv, settings)
  return(-deriv)
}

compute_obs_deriv <- function(v, k, w_counts, totals, mean_deriv_mtx, deriv, settings){
  init_mult <- 1000
  Ts <- settings$dim$Ts
  mean <- mean.global$mean[[k]][v,]
  variance <- var.global$var[[k]][v,]
  chain_variance <- settings$beta$chain_variance
  # temporary value of zeta
  temp_vect <- rep(0, Ts)
  for (t in 1:Ts) {
    temp_vect[t] <- exp(mean[t+1]+variance[t+1]/2)
  }
  for (t in 1:Ts) {
    mean_deriv <- mean_deriv_mtx[t,]
    term1 <- 0
    term2 <- 0
    term3 <- 0
    term4 <- 0
    for (u in 2:(Ts+1)) {
      mean_u <- mean[u]
      mean_u_prev <- mean[u-1]
      dmean_u <- mean_deriv[u]
      dmean_u_prev <- mean_deriv[u-1]
      
      term1 <- term1+(mean_u-mean_u_prev)*(dmean_u-dmean_u_prev)
      term2 <- term2+(w_counts[u-1]-(totals[u-1]*temp_vect[u-1]/zeta.global[[k]][u-1]))*dmean_u
    }
    if(chain_variance) {
      term1 <- -(term1/chain_variance)
      term1 <- term1-(mean[1]*mean_deriv[1])/(init_mult*chain_variance)
    }
    deriv[t] <- term1+term2
  }
  return(deriv)
}

compute_mean_deriv <- function(k, v,time,deriv,settings) {
  Ts <- settings$dim$Ts
  obs_variance <- settings$beta$obs_variance
  chain_variance <- settings$beta$chain_variance
  fwd_var_ <- var.global$fwd_var[[k]][v,]
  deriv[1] <- 0
  for (t in 2:(Ts+1)) {
    if (obs_variance>0) {
      w <- obs_variance/(fwd_var_[t-1]+chain_variance+obs_variance)
    }
    else w <- 0
    val <- w*deriv[t-1]
    if (time == t-1) {
      val <- val+(1-w)
    }
    deriv[t] <- val
  }
  for (t in Ts:1) {
    if (chain_variance == 0) {
      w <- 0
    }
    else {
      w <- chain_variance/(fwd_var_[t]+chain_variance)
    }
    deriv[t] <- w*deriv[t]+(1-w)*deriv[t+1]
  }
  return(deriv)
}

compute_bound <- function(k, beta.k, time_slices, settings) {
  V <- settings$dim$V
  Ts <- settings$dim$Ts
  term_1 = 0
  term_2 = 0
  term_3 = 0
  
  val = 0
  ent = 0
  totals <- summary(factor(time_slices))
  
  chain_variance <- settings$beta$chain_variance
  # compute mean, fwd_mean; they'are not updating here
  #update_post_mean(k, beta.k, settings)
  #update_zeta(k, settings)
  
  var.k <- var.global$var[[k]]
  mean.k <- mean.global$mean[[k]]
  val <- 0
  for (v in 1:V) {
    val <- val+var.k[v,1]-var.k[v,Ts+1]
  }
  val <- val/2*chain_variance
  for (t in 2:(Ts+1)) {
    term_1 <- 0
    term_2 <- 0
    ent <- 0
    for (v in 1:V) {
      m <- mean.k[v,t]
      prev_m <- mean.k[v,t-1]
      
      v_ <- var.k[v,t]
      term_1 <- term_1 + (m-prev_m)^2/(2*chain_variance)-v/chain_variance-log(chain_variance)
      term_2 <- term_2 + beta.k[t-1,v]*m
      ent = ent + log(v_)/2
    }
    term_3 = -totals[t-1] * log(zeta.global[[k]][t-1])
    val = val + term_2 + term_3 + ent - term_1
  }
  return(val)
}

update_zeta <- function(k, settings) {
  Ts <- settings$dim$Ts
  V <- settings$dim$V
  K <- settings$dim$K
  zeta.k <- zeta.global[[k]]
  mean.k <- mean.global$mean[[k]]
  var.k <- var.global$var[[k]]
  for (t in 1:Ts) {
    zeta.k[t] <- sum(exp(mean.k[,t+1]+var.k[,t+1]/2))
  }
  zeta.global[[k]] <<- zeta.k
}

update_post_mean <- function(v, k, beta.k, settings) {
  mean.k <- mean.global$mean[[k]]
  fwd_mean.k <- mean.global$fwd_mean[[k]]
  var.k <- var.global$var[[k]]
  fwd_var.k <- var.global$fwd_var[[k]]
  
  Ts <- settings$dim$Ts
  V <- settings$dim$V
  obs_variance <- settings$beta$obs_variance
  chain_variance <- settings$beta$chain_variance
  #for (v in 1:V) {
  obs <- obs.global_byK[[k]][,v]
  fwd_mean.k[v,1] <- 0
  # forward
  for (t in 2:(Ts+1)) {
    c <- obs_variance/(fwd_var.k[v,t-1]+chain_variance+obs_variance)
    fwd_mean.k[v,t] <- c*fwd_mean.k[v, t-1]+(1-c)*obs[t-1]
  }
  # backward
  mean.k[v,Ts+1] <- fwd_mean.k[v, Ts+1]
  for (t in Ts:1) {
    c <- chain_variance/(fwd_var.k[v,t]+chain_variance)
    mean.k[v,t] <- c*fwd_mean.k[v,t]+(1-c)*mean.k[v, t+1]
  }
  #}
  
  mean.global$mean[[k]] <<- mean.k
  mean.global$fwd_mean[[k]] <<- fwd_mean.k
}



