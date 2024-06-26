# Some functions borrowed from stm package.
## Roberts, Margaret E., Brandon M. Stewart, and Dustin Tingley. 
## "Stm: An r Package for Structural Topic Models."
## Journal of Statistical Software 91 (October 31, 2019): 1–40. 
## https://doi.org/10.18637/jss.v091.i02.


#' STM Corpus Coercion
#'
#' Convert a set of document term counts and associated metadata to
#' the form required for processing by the \code{\link{stm}} function.
#'
#' @param documents A documents-by-term matrix of counts, or a set of
#' counts in the format returned by \code{\link{prepDocuments}}. Supported
#' matrix formats include \pkg{quanteda} \link[quanteda]{dfm}
#' and \pkg{Matrix} sparse matrix objects in \code{"dgCMatrix"} or
#' \code{"dgTMatrix"} format.
#'
#' @param vocab Character vector specifying the words in the corpus in the
#' order of the vocab indices in documents. Each term in the vocabulary index
#' must appear at least once in the documents.  See \code{\link{prepDocuments}}
#' for dropping unused items in the vocabulary.  If \code{documents} is a
#' sparse matrix or \pkg{quanteda} \link[quanteda]{dfm} object, then \code{vocab} should not
#'  (and must not) be supplied.  It is contained already inside the column
#'  names of the matrix.
#'
#' @param data An optional data frame containing the prevalence and/or content
#' covariates.  If unspecified the variables are taken from the active
#' environment.
#'
#' @param \dots Additional arguments passed to or from other methods.
#'
#' @return A list with components \code{"documents"}, \code{"vocab"}, and
#' \code{"data"} in the form needed for further processing by the \code{stm}
#' function.
#'
#' @seealso \code{\link{prepDocuments}}, \code{\link{stm}}
#'
#' @examples
#' \donttest{
#' library(quanteda)
#' gadarian_corpus <- corpus(gadarian, text_field = "open.ended.response")
#' gadarian_dfm <- dfm(gadarian_corpus, 
#'                      remove = stopwords("english"),
#'                      stem = TRUE)
#' asSTMCorpus(gadarian_dfm)
#' }
asSTMCorpus <- function(documents, vocab, data = NULL, ...) {
  UseMethod("asSTMCorpus")
}

#' @method asSTMCorpus list
#' @keywords internal
asSTMCorpus.list <- function(documents, vocab=NULL, data = NULL, ...) {
  list(documents = documents, vocab = vocab, data = data)
}

#' @method asSTMCorpus dfm
#' @keywords internal
asSTMCorpus.dfm <- function(documents, vocab, data = NULL, ...) {
  if (!missing(vocab)) {
    # in case K was not specified by name, and it was confused with the
    # vocab argument (missing for dfm inputs)
    if (is.numeric(vocab) & length(vocab)==1) {
      stop("incorrect argument type for vocab, did you mean to specify K = ", vocab, "?")
    } else {
      stop("if documents is a dfm, do not specify vocab separately")
    }
  }
  
  # convert the dfm input as the first argument into the structure of the
  # older function where this is split into a list
  dfm_stm <- quanteda::convert(documents, to = "stm", docvars = data)
  if(is.null(data)) data <- dfm_stm[["meta"]]
  
  list(documents = dfm_stm[["documents"]], vocab = dfm_stm[["vocab"]],
       data = data)
}

#' @method asSTMCorpus dgCMatrix
#' @keywords internal
asSTMCorpus.dgCMatrix <- function(documents, vocab, data = NULL, ...) {
  if (!missing(vocab)) {
    # in case K was not specified by name, and it was confused with the
    # vocab argument (missing for dfm inputs)
    if (is.numeric(vocab) & length(vocab)==1) {
      stop("incorrect argument type for vocab, did you mean to specify K = ", vocab, "?")
    } else {
      stop("if documents is a matrix, do not specify vocab separately")
    }
  }
  
  # drop unused terms
  tot <- Matrix::colSums(documents)
  unseen <- which(tot == 0)
  if (length(unseen) > 0) {
    warning(sprintf("dropping %d unseen terms from the vocabulary", length(unseen)))
    documents <- documents[, -unseen, drop = FALSE]
  }
  
  # convert to docs-by-terms, one column per doc
  x <- as(Matrix::t(documents), "dgCMatrix")
  
  # get the vocab
  vocab <- rownames(x)
  
  # extract the term count information
  terms <- x@i + 1L
  counts <- as.integer(x@x)
  off <- x@p[-length(x@p)] + 1L
  len <- x@p[-1L] - x@p[-length(x@p)]
  
  # fill in the documents
  documents <- vector("list", ncol(x))
  names(documents) <- colnames(x)
  for (i in seq_along(documents)) {
    ix <- seq.int(off[[i]], length.out = len[[i]])
    documents[[i]] <- matrix(c(terms[ix], counts[ix]), nrow = 2L,
                             byrow = TRUE)
  }
  
  list(documents = documents, vocab = vocab, data = data)
}

#' @method asSTMCorpus dgTMatrix
#' @keywords internal
asSTMCorpus.dgTMatrix <- function(documents, vocab, data = NULL, ...) {
  documents <- methods::as(documents,"dgCMatrix")
  asSTMCorpus.dgCMatrix(documents, vocab, data, ...)
}

#' Make a B-spline Basis Function
#' 
#' This is a simple wrapper around the \code{\link[splines]{bs}} function in
#' the splines package.  It will default to a spline with 10 degrees of
#' freedom.
#' 
#' This is a simple wrapper written as users may find it easier to simply type
#' \code{s} rather than selecting parameters for a spline. We also include
#' \code{predict} and \code{makepredictcall} generic functions for the class
#' so it will work in settings where \code{\link{predict}} is called.
#' 
#' @param x The predictor value.
#' @param df Degrees of freedom.  Defaults to the minimum of 10 or one minus
#' the number of unique values in x.
#' @param \dots Arguments passed to the \code{\link[splines]{bs}} function.
#' @return A predictor matrix of the basis functions.
#' @seealso \code{\link[splines]{bs}} \code{\link[splines]{ns}}
#' @export
s <- function(x, df, ...) {
  if(inherits(x,"Date")) {
    warning("A Date object coerced to numeric. 
            Converting variable in advance will stop this warning in the future.
            Postprocessing tools may not work with dates.")
    x <- as.numeric(x)
  }
  
  nval <- length(unique(x))
  if(missing(df)) {
    df <- min(10, (nval-1))
  }
  obj <- splines::bs(x, df,...)
  attr(obj, "class") <- c("s", attr(obj, "class")) #we need this to ensure that our predict generics trigger
  return(obj)
}

#' @keywords internal
predict.s <- function (object, newx, ...) 
{
  if (missing(newx)) 
    return(object)
  a <- c(list(x = newx), attributes(object)[c("degree", "knots", 
                                              "Boundary.knots", "intercept")])
  do.call("splines::bs", a)
}

#' @keywords internal
makepredictcall.s <- function (var, call) 
{
  #if (as.character(call)[1L] != "bs") 
  #  return(call)
  at <- attributes(var)[c("degree", "knots", "Boundary.knots", 
                          "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

read.slam <- function(corpus) {
  #convert a simple triplet matrix to list format.
  if(!inherits(corpus, "simple_triplet_matrix")) stop("corpus is not a simple triplet matrix")
  if (inherits(corpus,"TermDocumentMatrix")) {
    non_empty_docs <- which(slam::col_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[,non_empty_docs]$j, corpus[,non_empty_docs]$i, corpus[,non_empty_docs]$v) 
    names(documents) <- corpus[,non_empty_docs]$dimnames$Docs
  } else {
    non_empty_docs <- which(slam::row_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[non_empty_docs,]$i, corpus[non_empty_docs,]$j, corpus[non_empty_docs,]$v) 
    names(documents) <- corpus[non_empty_docs,]$dimnames$Docs
  }
  vocab <- corpus$dimnames$Terms
  return(list(documents=documents,vocab=vocab))
}


#Our Format to Triplet Format
doc.to.ijv <- function(documents, fixzeroindex=TRUE) {
  #Turns our format into triplet format (can be zero indexed)
  indices <- unlist(lapply(documents, '[',1,)) #grab the first row
  if((0 %in% indices) & fixzeroindex) indices <- indices + 1 #if zero-indexed, fix it.
  counts <- lapply(documents, '[',2,)  #grab the second row but preserve the list structure for a moment
  VsubD <- unlist(lapply(counts,length)) #grab the number of unique words per document
  rowsums <- unlist(lapply(counts,sum)) #grab the number of tokens per documents
  docids <- rep(1:length(documents), times=VsubD) #add row numbers
  counts <- unlist(counts) #unlist the count structure
  #now we return using the convention for the simple triplet matrix,
  #plus the row sums which we use in DMR.
  return(list(i=as.integer(docids), j=as.integer(indices), v=as.integer(counts), rowsums=as.integer(rowsums)))
}

#Triplet Format to our Document Format
ijv.to.doc <- function(i,j,v) {
  index <- split(j,i)
  index <- lapply(index,as.integer)
  count <- split(v,i)
  count <- lapply(count,as.integer)
  mapply(rbind,index,count, SIMPLIFY=FALSE)
}

recoverL2 <- function(Qbar, anchor, wprob, verbose=TRUE, recoverEG=TRUE, ...) {
  #NB: I've edited the script to remove some of the calculations by commenting them
  #out.  This allows us to store only one copy of Q which is more memory efficient.
  #documentation for other pieces is below.
  
  #Qbar <- Q/rowSums(Q)
  X <- Qbar[anchor,]
  XtX <- tcrossprod(X)
  
  #In a minute we will do quadratic programming
  #these jointly define the conditions.  First column
  #is a sum to 1 constraint.  Remainder are each parameter
  #greater than 0.
  Amat <- cbind(1,diag(1,nrow=nrow(X)))
  bvec <- c(1,rep(0,nrow(X)))
  
  #Word by Word Solve For the Convex Combination
  condprob <- vector(mode="list", length=nrow(Qbar))
  for(i in 1:nrow(Qbar)) {
    if(i %in% anchor) { 
      #if its an anchor we create a dummy entry that is 1 by definition
      vec <- rep(0, nrow(XtX))
      vec[match(i,anchor)] <- 1
      condprob[[i]] <- vec
    } else {
      y <- Qbar[i,]
      
      if(recoverEG) {
        solution <- expgrad(X,y,XtX, ...)$par
      } else {
        #meq=1 means the sum is treated as an exact equality constraint
        #and the remainder are >=
        solution <- quadprog::solve.QP(Dmat=XtX, dvec=X%*%y, 
                                       Amat=Amat, bvec=bvec, meq=1)$solution  
      }
      
      if(any(solution <= 0)) {
        #we can get exact 0's or even slightly negative numbers from quadprog
        #replace with machine double epsilon
        solution[solution<=0] <- .Machine$double.eps
      } 
      condprob[[i]] <- solution
    }
  }
  #Recover Beta (A in this notation)
  #  Now we have p(z|w) but we want the inverse
  weights <- do.call(rbind, condprob)
  A <- weights*wprob
  A <- t(A)/colSums(A)
  
  #Recover The Topic-Topic Covariance Matrix
  #Adag <- mpinv(A)  
  #R <- t(Adag)%*%Q%*%Adag
  return(list(A=A))
  #return(list(A=A, R=R, condprob=condprob))
}

expgrad <- function(X, y, XtX=NULL, alpha=NULL, tol=1e-7, max.its=500) {
  if(is.null(alpha)) alpha <- 1/nrow(X) 
  alpha <- matrix(alpha, nrow=1, ncol=nrow(X))
  if(is.null(XtX)) XtX <- tcrossprod(X)
  
  ytX <- y%*%t(X)
  converged <- FALSE
  eta <- 50
  sse.old <- Inf
  its <- 1
  while(!converged) {
    #find the gradient (y'X - alpha'X'X)
    grad <- (ytX - alpha%*%XtX) #101-105
    sse <- sum(grad^2) #106, sumSquaredError
    grad <- 2*eta*grad
    maxderiv <- max(grad)    
    
    #update parameter 
    alpha <- alpha*exp(grad-maxderiv)
    #project parameter back to space
    alpha <- alpha/sum(alpha)
    
    converged <- abs(sqrt(sse.old)-sqrt(sse)) < tol
    if(its==max.its) break
    sse.old <- sse
    its <- its + 1
  } 
  entropy <- -1*sum(alpha*log(alpha))
  return(list(par=as.numeric(alpha), its=its, converged=converged,
              entropy=entropy, log.sse=log(sse)))
}


####
# Non-Exported Utility Functions
####

##
# Random Utilities

#function for collapsing a character vector to a comma separated list
#note that we eliminate things with zero length so that we can return
#fewer than n words and still have the lists look nice
commas <- function(text){  
  paste(text[nchar(text)>0], collapse=", ")
}

#from the R documentation for is.integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
} 

posint <- function(x) {
  all(is.wholenumber(x)) & all(x>0)
}

nonnegint <- function(x) {
  all(is.wholenumber(x)) & all(x>=0)
}


#' Draw from a Multivariate Normal
#' 
#' A basic function for doing multivariate normal simulations
#' via the cholesky decomposition of the covariance matrix. Function
#' is based on one by Peter Hoff.
#' 
#' This is a pretty standard multivariate normal generator. It could
#' almost certainly be faster if we ported it over to \pkg{RcppArmadillo}
#' but it isn't used a ton at the moment.
#' 
#' @param n number of draws
#' @param mu the K-dimensional mean
#' @param Sigma the K by K dimensional positive definite covariance matrix
#' @param chol.Sigma the cholesky decomposition of the Sigma matrix.
#' 
#' @keywords internal
rmvnorm<-function(n,mu,Sigma,chol.Sigma=chol(Sigma)) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol.Sigma) +c(mu))
}

#' Calculate FREX (FRequency and EXclusivity) Words
#' 
#' A primarily internal function for calculating FREX words.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' FREX attempts to find words which are both frequent in and exclusive to a topic of interest.
#' Balancing these two traits is important as frequent words are often by themselves simply functional
#' words necessary to discuss any topic.  While completely exclusive words can be so rare as to not
#' be informative. This accords with a long-running trend in natural language processing which is best exemplified
#' by the Term frequency-Inverse document frequency metric.  
#' 
#' Our notion of FREX comes from a paper by Bischof and Airoldi (2012) which proposed a Hierarchical
#' Poisson Deconvolution model.  It relies on a known hierarchical structure in the documents and requires
#' a rather complicated estimation scheme.  We wanted a metric that would capture their core insight but
#' still be fast to compute.
#' 
#' Bischof and Airoldi consider as a summary for a word's contribution to a topic the harmonic mean of the
#' word's rank in terms of exclusivity and frequency.  The harmonic mean is attractive here because it 
#' does not allow a high rank along one of the dimensions to compensate for the lower rank in another. Thus
#' words with a high score must be high along both dimensions.
#' 
#' The formula is ' 
#'\deqn{FREX = \left(\frac{w}{F} + \frac{1-w}{E}\right)^{-1}}{FREX = ((w/F) + ((1-w)/E))^-1} 
#' where F is the frequency score given by the empirical CDF of the word in it's topic distribution.  Exclusivity
#' is calculated by column-normalizing the beta matrix (thus representing the conditional probability of seeing
#' the topic given the word).  Then the empirical CDF of the word is computed within the topic.  Thus words with
#' high values are those where most of the mass for that word is assigned to the given topic.
#' 
#' For rare words exclusivity will always be very high because there simply aren't many instances of the word.
#' If \code{wordcounts} are passed, the function will calculate a regularized form of this distribution using a
#' James-Stein type estimator described in \code{\link{js.estimate}}.
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @param w a value between 0 and 1 indicating the proportion of the weight assigned to frequency 
#' @param wordcounts a vector of word counts.  If provided, a James-Stein type shrinkage estimator is 
#' applied to stabilize the exclusivity probabilities. This helps with the concern that the rarest words
#' will always be completely exclusive.
#' @references 
#' Bischof and Airoldi (2012) "Summarizing topical content with word frequency and exclusivity"
#' In Proceedings of the International Conference on Machine Learning.
#' @seealso \code{\link{labelTopics}} \code{\link{js.estimate}}
#' @export
#' @keywords internal
calcfrex <- function(logbeta, w=.5, wordcounts=NULL) {
  excl <- t(t(logbeta) - col.lse(logbeta))
  if(!is.null(wordcounts)) {
    #if word counts provided calculate the shrinkage estimator
    excl <- safelog(sapply(1:ncol(excl), function(x) js.estimate(exp(excl[,x]), wordcounts[x])))
  } 
  freqscore <- apply(logbeta,1,data.table::frank)/ncol(logbeta)
  exclscore <- apply(excl,1,data.table::frank)/ncol(logbeta)
  frex <- 1/(w/freqscore + (1-w)/exclscore)
  apply(frex,2,order,decreasing=TRUE)
}

#' A James-Stein Estimator Shrinking to a Uniform Distribution
#' 
#' A primarily internal function used in \code{\link{calcfrex}}.
#' 
#' This calculates a James-Stein type shrinkage estimator for a discrete probability
#' distribution regularizing towards a uniform distribution. The amount of shrinkage
#' is a function of the variance of MLE and the L2 norm distance from the uniform.
#' 
#' This function is based off the ideas in Hausser and Strimmer (2009)
#' 
#' @param prob the MLE estimate of the discrete probability distribution
#' @param ct the count of words observed to estimate that distribution
#' 
#' @references 
#' Hausser, Jean, and Korbinian Strimmer. "Entropy inference and the James-Stein estimator, 
#' with application to nonlinear gene association networks." Journal of Machine Learning Research 
#' 10.Jul (2009): 1469-1484.
#' @export
#' @keywords internal
js.estimate <- function(prob, ct) {
  if(ct<=1) {
    #basically if we only observe a count of 1
    #the variance goes to infinity and we get the uniform distribution.
    return(rep(1/length(prob), length(prob)))
  }
  # MLE of prob estimate
  mlvar <- prob*(1-prob)/(ct-1)
  unif <- rep(1/length(prob), length(prob)) 
  
  # Deviation from uniform
  deviation <- sum((prob-unif)^2)
  
  #take care of special case,if no difference it doesn't matter
  if(deviation==0) return(prob)
  
  lambda <- sum(mlvar)/deviation
  #if despite  our best efforts we ended up with an NaN number-just return the uniform distribution.
  if(is.nan(lambda)) return(unif)
  
  #truncate
  if(lambda>1) lambda <- 1
  if(lambda<0) lambda <- 0
  
  #Construct shrinkage estimator as convex combination of the two
  lambda*unif + (1 - lambda)*prob
}

#' Calculate Lift Words
#' 
#' A primarily internal function for calculating words according to the lift metric.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' Lift is the calculated by dividing the topic-word distribution by the empirical
#' word count probability distribution.  In other words the Lift for word v in topic
#' k can be calculated as:
#' 
#' \deqn{Lift = \beta_{k,v}/(w_v/\sum_v w_v)}{Lift = \beta/wbar} 
#' 
#' We include this after seeing it used effectively in Matt Taddy's work including his
#' excellent \pkg{maptpx} package. Definitions are given in Taddy(2012).
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @param wordcounts a V length vector indicating the number of times each word appears in the corpus. 
#' @references 
#' Taddy, Matthew. 2012. "On Estimation and Selection for Topic Models." AISTATS JMLR W&CP 22
#' 
#' @seealso \code{\link{labelTopics}}
#' @export
#' @keywords internal
calclift <- function(logbeta, wordcounts) {
  emp.prob <- log(wordcounts) - log(sum(wordcounts))
  lift <- logbeta - rep(emp.prob, each=nrow(logbeta)) 
  apply(lift, 1, order, decreasing=TRUE)
}

#' Calculate Score Words
#' 
#' A primarily internal function for calculating words according to the score metric.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' Score is a metric which we include because it is used effectively in the 
#' \pkg{lda} package by Jonathan Chang. It is calculated as:
#' \deqn{\beta_{v, k} (\log \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{v,k'})}
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @references 
#' Jonathan Chang (2015). lda: Collapsed Gibbs Sampling Methods for Topic Models. R package version 1.4.2.
#' https://CRAN.R-project.org/package=lda 
#' @seealso \code{\link{labelTopics}} 
#' @export
#' @keywords internal
calcscore <- function(logbeta) { 
  ldascore <- exp(logbeta)*(logbeta - rep(colMeans(logbeta), each=nrow(logbeta)))
  apply(ldascore, 1, order, decreasing=TRUE)
} 
##
#Document convertors
##

#Our Format to Triplet Format
doc.to.ijv <- function(documents, fixzeroindex=TRUE) {
  #Turns our format into triplet format (can be zero indexed)
  indices <- unlist(lapply(documents, '[',1,)) #grab the first row
  if((0 %in% indices) & fixzeroindex) indices <- indices + 1 #if zero-indexed, fix it.
  counts <- lapply(documents, '[',2,)  #grab the second row but preserve the list structure for a moment
  VsubD <- unlist(lapply(counts,length)) #grab the number of unique words per document
  rowsums <- unlist(lapply(counts,sum)) #grab the number of tokens per documents
  docids <- rep(1:length(documents), times=VsubD) #add row numbers
  counts <- unlist(counts) #unlist the count structure
  #now we return using the convention for the simple triplet matrix,
  #plus the row sums which we use in DMR.
  return(list(i=as.integer(docids), j=as.integer(indices), v=as.integer(counts), rowsums=as.integer(rowsums)))
}

#Triplet Format to our Document Format
ijv.to.doc <- function(i,j,v) {
  index <- split(j,i)
  index <- lapply(index,as.integer)
  count <- split(v,i)
  count <- lapply(count,as.integer)
  mapply(rbind,index,count, SIMPLIFY=FALSE)
}

##
#A series of fast softmax functions mostly wrappers around matrixStats package functions
##
logsoftmax <- function(x) {
  x - lse(x)
}

lse <- function(x) {
  matrixStats::logSumExp(x)
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}
col.lse <- function(mat) {
  matrixStats::colLogSumExps(mat)
}

softmax <- function(x) {
  exp(x - lse(x))
}

safelog <- function(x, min=-1000) {
  out <- log(x)
  out[which(out< min)] <- min
  out
}


# Note: I started documenting this but am not exporting because I would need to appropriately
# generalize.  It is currently only set up to do the mgaussian one I think- but I haven't looked
# at it in a while.

#' Unpack a \pkg{glmnet} object
#' 
#' A function to quickly unpack a \pkg{glmnet} model object and calculate an
#' optimal model from the regularization path.
#' 
#' This is a small utility we wrote to deal with the slow methods dispatch for S4
#' classes.  The more straightforward option is the \code{coef()} method for \pkg{glmnet}
#' objects but when trying to make thousands of calls a second, that can be very slow
#' 
#' @param mod the glmnet model
#' @param ic.k the information criterion value.  AIC is \code{ic.k=2} and BIC would be \code{ic.k=log n}
#' 
#' @return 
#' A list
#' \item{coef}{a matrix of coefficients}
#' \item{intercept}{the intercepts}
#' @keywords internal
unpack.glmnet <- function(mod, ic.k) {
  dev <- (1-mod$dev.ratio)*mod$nulldev
  df <- colSums(mod$dfmat)
  ic <- dev + ic.k*df
  lambda <- which.min(ic)
  
  #methods dispatch here is crazy, so define out own function
  subM <- function(x, p) {
    ind <- (x@p[p]+1):x@p[p+1]
    rn <- x@i[ind]+1
    y <- x@x[ind]
    out <- rep(0, length=nrow(x))
    out[rn] <- y
    out
  }
  coef <- lapply(mod$beta, subM, lambda) #grab non-zero coefs
  coef <- do.call(cbind,coef)
  intercept <- mod$a0[,lambda]
  return(list(coef=coef, intercept=intercept))
}

gram <- function(mat) {
  nd <- Matrix::rowSums(mat)
  mat <- mat[nd>=2,] #its undefined if we don't have docs of length 2
  nd <- nd[nd>=2]
  divisor <- nd*(nd-1)
  Q <- Matrix::crossprod(mat/sqrt(divisor)) - Matrix::Diagonal(x=Matrix::colSums(mat/divisor))
  return(as.matrix(Q))
}


fastAnchor <- function(Qbar, K, verbose=TRUE) {
  basis <- c()
  rowSquaredSums <- rowSums(Qbar^2) #StabilizedGS
  
  for(i in 1:K) {
    basis[i] <- which.max(rowSquaredSums) #83-94
    
    maxval <- rowSquaredSums[basis[i]]
    normalizer <- 1/sqrt(maxval) #100
    
    #101-103
    Qbar[basis[i],] <- Qbar[basis[i],]*normalizer 
    
    #For each row
    innerproducts <- Qbar%*%Qbar[basis[i],] #109-113
    
    #Each row gets multiplied out through the Basis
    project <- as.numeric(innerproducts)%o%Qbar[basis[i],] #118
    
    #Now we want to subtract off the projection but
    #first we should zero out the components for the basis
    #vectors which we weren't intended to calculate
    project[basis,] <- 0 #106 (if statement excluding basis vectors)
    Qbar <- Qbar - project #119
    rowSquaredSums <- rowSums(Qbar^2)
    rowSquaredSums[basis] <- 0 #here we cancel out the components we weren't calculating.
  }
  return(basis)
}

