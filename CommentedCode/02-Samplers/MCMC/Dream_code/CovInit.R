##' Alternative sampling method
##'  used in examples, may not be useful in practice
##'  @param pars list of parameter vectors
##'  @param nseq scalar
##'  @param muX vector of length ndim
##'  @param qcov
##'  @param bound.handling. character one of: none,reflect,bound,fold,rand
##'  @return x
##'    matrix dim nseq x ndim. parameter wise range [xmin,xmax] assured by handleBounds
##
## depends on:
##   handleBounds.R
##   rand_utils.R
##
##' When used as input to dream, muX, qcov and bound.handling must be passed as extra parameterrs
CovInit <- function(pars,nseq,muX,qcov,bound.handling)
{

  ##[x] = repmat(Extra.muX,MCMCPar.seq,1) + randn(MCMCPar.seq,MCMCPar.n) * chol(Extra.qcov);

  ## Components verified to match matlab
  ## print(t(matrix(rep(muX,nseq),length(muX))))
  ## print(chol(qcov))
  
  x <- t(matrix(rep(muX,nseq),length(muX)))+randn(nseq,length(pars)) %*% chol(qcov)
  
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  x <- handleBounds(x,lower,upper,bound.handling)
  return(x)
}
