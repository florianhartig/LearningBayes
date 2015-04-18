##' Latin Hypercube sampling
##' @param pars list. length ndim - one element per parameter
##' @param nseq scalar = number of sequences
##' @return x matrix nseq x ndim. parameter-wise range [xmin, xmax]
##
## depends on rand_utils.R
LHSInit <- function(pars, nseq){
  ## TODO: R approach?
  ## lapply(pars, function(r)
  ##        sample(seq(min(r), max(r), length = nseq))
  ##        )

  ## dimensions
  ## xmin,xmax. vectors of length ndim
  ## ran. matrix nseq x ndim. range [0,1]
  ## idx. vector length nseq range [1,nseq]
  ## P. vector length nseq. range [0,1]
  
  xmin <- sapply(pars, function(x) min(x[[1]]))
  xmax <- sapply(pars, function(x) max(x[[1]]))

  ## Number of parameters
  ndim <- length(pars)
  ## Initialize array ran with random numbers
  ran <- rand(nseq,ndim)
  
  ## Initialize array s with zeros
  s <- matrix(0,nseq,ndim)
  
  ## Fill s by iterating through parameters
  for (j in 1:ndim){
    ## Random permutation
    idx <- sample(nseq)
    P <- (idx-ran[,j])/nseq
    s[,j] <- xmin[j]+P*(xmax[j]-xmin[j])
  } ##for pars

  return(s)
} ## LHSInit
