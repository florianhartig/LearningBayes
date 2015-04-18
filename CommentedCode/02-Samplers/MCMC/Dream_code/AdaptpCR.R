##' Updates the probabilities of the various crossover values
##'
##' @param CR matrix nseq x steps
##' @param delta.tot vector of length nCR
##' @param lCR.old vector of length nCR
##' @param control. list needs elements: nCR,nseq
##'
##' @return ... list with elements
##'   pCR vector of length nCR
##'   lCR vector of length nCR
AdaptpCR <- function(CR,delta.tot,lCR.old,control){

  if(!any(delta.tot>0)) stop("AdaptpCR: no changes in X, would cause NaN in pCR")
  
  ## dimensions:
  ##  CR vector of length nseq*steps
  ##  zz iter. 1:nCR
  ##  cr.count. scalar. range [0,nseq*steps]
  
  ## Make CR to be a single vector
  CR <- c(CR)
  
  ## Determine lCR
  lCR <- rep(NA,control$nCR)
  for (zz in 1:control$nCR){
    
    ## Determine how many times a particular CR value is used
    ## TODO: shouldn't this be which(CR==CR[zz]])?
    cr.count <- length(which(CR==zz/control$nCR))
    
    ## This is used to weight delta.tot
    lCR[zz] <- lCR.old[zz]+cr.count
  }                                     #for CRs
  
  ## Adapt pCR using information from averaged normalized jumping distance
  pCR <- control$nseq * (delta.tot / lCR) / sum(delta.tot)
    
  ## Normalize pCR
  pCR <- pCR/sum(pCR)

  return(list(pCR=pCR,lCR=lCR))
} ##AdaptpCR
