#' @param documents The documents passed to the initialization
#' @param settings The initialized settings in settings_init.R
cdtm_init <- function(documents, settings) {
  K <- settings$dim$K
  V <- settings$dim$V
  N <- settings$dim$N
  Ts <- settings$dim$Ts
  alpha <- settings$init$alpha 
  eta <- settings$init$eta 
  verbose <- settings$verbose
  # Prepare the Gram matrix
  docs <- doc.to.ijv(documents)
  mat <- Matrix::sparseMatrix(docs$i,docs$j, x=docs$v)
  rm(docs)
  wprob <- Matrix::colSums(mat)
  wprob <- wprob/sum(wprob)
  keep <- NULL
  # initialize all the beta
  Q <- gram(mat)
  Qsums <- rowSums(Q)
  if(any(Qsums==0)) {
    # if there are zeroes, we want to remove them for just the anchor word procedure.
    temp.remove <- which(Qsums==0)
    if(is.null(keep)) {
      keep <- which(Qsums!=0)
    } else {
      keep <- keep[which(Qsums!=0)]
    }
    Q <- Q[-temp.remove,-temp.remove]
    Qsums <- Qsums[-temp.remove]
    wprob <- wprob[-temp.remove]
  }
  Q <- Q/Qsums
  anchor <- fastAnchor(Q, K=K, verbose=verbose)
  beta <- recoverL2(Q, anchor, wprob, verbose=verbose, recoverEG=settings$init$recoverEG)$A
  
  # generate other parameters
  mu <- matrix(0, nrow=(K-1),ncol=1)
  if(settings$gamma$mod=="CDTM"){  sigma <- rep(list(diag(20, nrow=(K-1))), Ts)
  }else sigma <- diag(20, nrow=(K-1))
  lambda <- matrix(0, nrow=N, ncol=(K-1))
  # repeat beta for time_slices
  beta <- rep(list(beta), Ts)
  ## generate the variance/fwd_variance
  var <- rep(list(matrix(rep(0, (Ts+1)*V), ncol = Ts+1)), K)
  fwd_var <- rep(list(matrix(rep(0, (Ts+1)*V), ncol = Ts+1)), K)
  # generate the mean/fwd_mean
  mean <- rep(list(matrix(rep(0, (Ts+1)*V), ncol = Ts+1)), K)
  fwd_mean <- rep(list(matrix(rep(0, (Ts+1)*V), ncol = Ts+1)), K)
  ## generate zeta (variational parameter)
  zeta <- rep(list(matrix(rep(0, Ts))), K)
  
  model <- list(mu=mu, sigma=sigma, beta=beta, lambda=lambda, 
                var=list(var=var, fwd_var=fwd_var), 
                mean=list(mean=mean, fwd_mean=fwd_mean), 
                zeta=zeta)
  
  return(model)
}

