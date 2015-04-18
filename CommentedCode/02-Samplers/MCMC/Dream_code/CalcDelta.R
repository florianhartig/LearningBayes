##' Calculate total normalized Euclidean distance for each crossover value

##' @param nCR. scalar
##' @param delta.tot vector of length nCR
##' @param delta.normX vector of length nseq
##' @param CR vector of length nseq

##' @return delta.tot vector of length nCR
CalcDelta <- function(nCR,delta.tot,delta.normX,CR){

  ##stopifnot(sum(delta.tot)>0 || sum(delta.normX)>0)
  
  ## Dimensions:
  ##  zz. iter 1:nCR
  ##  idx. vector. length [0,nseq]. range [1,nseq]
  
  ## Derive sum_p2 for each different CR value
  for (zz in 1:nCR){
    ## Find which chains are updated with zz/MCMCPar.nCR
    ## TODO: possible that floating point error prevents exact comparison?
    idx <- which(CR==zz/nCR)
    
    ## Add the normalized squared distance tot the current delta_tot;
    delta.tot[zz] <- delta.tot[zz]+sum(delta.normX[idx])
    
  } ## for CRs

  ##stopifnot(!any(is.na(delta.tot)) && sum(delta.tot)>0)
  return(delta.tot)
} ##CalcDelta
