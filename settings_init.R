#' @param out pre-processed documents
#' @param K number of topics
#' @param max.em.its maximum EM iteration 
#' @param obs_variance varsigma^2
#' @param chain_variance sigma^2
#' @param verbose verbose or not
#' @param seed seed
#' @param mod "CDTM" or "DTM"
settings_init <- function(out, K, prevalence=NULL, max.em.its=25, obs_variance=0.5, chain_variance=0.05, 
                          verbose=0, seed=1, mod="CDTM"){
  documents <- out$documents
  vocab <- out$vocab

  # dynamic model initialization
  time_slices <- out$meta$Time_slices-min(out$meta$Time_slices)
  num_time_slices <- length(unique(time_slices))
  data <- out$meta

  # Convergence tolerance. EM stops when relative change is below this level
  emtol <- 1e-5
  sigma.prior <- 0
  
  # turn to stm format as suggested in package stm
  args <- asSTMCorpus(documents, vocab, data)
  documents <- args$documents
  vocab <- args$vocab
  data <- args$data
  
  # length of Documents
  N <- length(documents)
  # to make it a vector: organize by document/word
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
  wcounts <- list(Group.1=sort(unique(wcountvec)))
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)
  # number of words
  V <- length(wcounts$Group.1)
  Ts <- length(unique(out$meta$Time_slices))
  if (mod=="CDTM" & !is.null(prevalence)){
    # deal with the prevalence covariate matrix
    termobj <- terms(prevalence, data = data)
    xmat <- sparse.model.matrix(termobj,data=data)
    xmat <- as.matrix(xmat)
    # to initialize all the settings
    covariates <- list(X=xmat, formula=prevalence)
  }else{
    covariates <- NULL
  }
  settings <- list(dim=list(K=K, V=V, N=N, Ts=Ts, wcounts=wcounts),
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol, 
                                    allow.neg.change=TRUE),
                   covariates=covariates,
                   gamma=list(maxits = 1000, mod=mod),
                   sigma=list(prior=sigma.prior),
                   init=list(alpha=(50/K), eta=.01, recoverEG=TRUE), 
                   beta=list(time_slices=time_slices, num_time_slices=num_time_slices,
                             obs_variance=obs_variance, chain_variance=chain_variance),
                   seed=seed,
                   verbose=verbose)
  return(settings)
}



