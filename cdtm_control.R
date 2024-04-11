cdtm_control <- function(documents, vocab, time_slices, settings, model=NULL) {
  globaltime <- proc.time()
  verbose <- settings$verbose
  # Step 1: Initialize Parameters
  mu <- list(mu=model$mu)
  sigma <- model$sigma
  beta <- list(beta=model$beta)
  lambda <- model$lambda
  convergence <- NULL
  
  #Pull out some book keeping elements
  ntokens <- sum(settings$dim$wcounts$x)
  stopits <- FALSE
  suffstats <- vector(mode="list", length=1)
  
  step <- 0
  
  # Step 2: Run EM
  while(!stopits) {
    t1 <- proc.time()
    # run the model
    if(verbose) print("Execute E-step...")
    suffstats <- estep(documents=documents, mod=settings$gamma$mod,
                       update.mu=(!is.null(mu$gamma)), beta=beta$beta,
                       lambda.old=lambda, mu=mu$mu, sigma=sigma, time_slices, step=step)
    t1 <- proc.time()
    sigma.ss <- suffstats$sigma
    lambda <- suffstats$lambda
    obs.ss <- suffstats$obs
    bound.ss <- suffstats$bound
    
    # do the m-step: optimize mu/sigma/beta
    if(verbose) print("M-step...update mu&sigma")
    mu <- opt.mu(gamma_prior=mu$gamma_prior,
                 lambda=lambda,
                 covar=settings$covariates$X,
                 maxits=settings$gamma$maxits,
                 time_slices=time_slices, mod=settings$gamma$mod)
    sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, time_slices=time_slices,
                       mu=mu$mu, sigprior=settings$sigma$prior, mod=settings$gamma$mod)
    # before updating beta, update the variance/fwd_variance first
    var.global <<- update_var(settings)
    # after we got the variance, we can iterate to update beta until it converged
    if(verbose) print("M-step...update beta")
    beta <- opt.beta(obs.ss, time_slices, settings)
    #Convergence
    convergence <- convergence.check(bound.ss, convergence, settings)
    if(verbose) print(paste("Convergence bound is: ",convergence$bound))
    stopits <- convergence$stopits
    step <- step + 1
  }
  # Step 3: Construct Output
  time <- (proc.time() - globaltime)[3]
  lambda <- cbind(lambda,0)
  model <- list(mu=mu, sigma=sigma, beta=list(obs=obs.ss, beta_byT=beta$beta, beta_byK=beta$beta_byK, logbeta=beta$logbeta), 
                settings=settings, 
                vocab=vocab, convergence=convergence,
                theta=exp(lambda - log(rowSums(exp(lambda)))),
                eta=lambda[,-ncol(lambda), drop=FALSE], time=time)
  
  class(model) <- "CDTM"
  return(model)
}