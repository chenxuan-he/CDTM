#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu_, lambda_, mu_, sigprior, time_slices, mod) {  
  if(mod == "DTM"){
    nu <- nu_
    mu <- mu_
    lambda <- lambda_
    if(ncol(mu)==1) {
      covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
    } else {
      covariance <- crossprod(lambda-t(mu)) 
    }
    sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
    sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    return(sigma)
  }
  sigma_ <- nu_
  Ts <- length(unique(time_slices))
  if(min(time_slices)==0) time_slices=time_slices+1
  for (tt in 1:Ts) {
    nu <- nu_[[tt]]
    lambda <- lambda_[which(time_slices==tt),]
    mu <- mu_[,which(time_slices==tt)]
    #find the covariance
    if(ncol(mu)==1) {
      covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
    } else {
      covariance <- crossprod(lambda-t(mu)) 
    }
    sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
    sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    sigma_[[tt]] <- sigma
  }
  return(sigma_)
}


