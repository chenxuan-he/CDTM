update_var <- function(settings) {
  # get the number of words and update it one by one
  fwd_var <- var.global$fwd_var
  var <- var.global$var
  V <- settings$dim$V
  K <- settings$dim$K
  Ts <- settings$dim$Ts
  chain_variance <- settings$beta$chain_variance
  obs_variance <- settings$beta$obs_variance
  INIT_VARIANCE_CONST <- 1000
  for (k in 1:K) {
    for (v in 1:V) {
      fwd_var[[k]][v,1] <- chain_variance * INIT_VARIANCE_CONST
      for (t in 2:(Ts+1)) {
        c <- obs_variance/(fwd_var[[k]][v,t-1]+chain_variance+obs_variance)
        fwd_var[[k]][v,t] <- c*(fwd_var[[k]][v,t-1]+chain_variance)
      }
      # backward pass
      var[[k]][v,Ts+1] = fwd_var[[k]][v,Ts+1]
      for (t in Ts:1) {
        if (fwd_var[[k]][v,t]>0) {
          c <- (fwd_var[[k]][v,t] / (fwd_var[[k]][v,t] + chain_variance))^2
        }
        else c <- 0
        var[[k]][v,t] <- (c * (var[[k]][v,t+1] - chain_variance)) + ((1 - c) * fwd_var[[k]][v,t])
      }
    }
  }
  return(list(var=var, fwd_var=fwd_var))
}