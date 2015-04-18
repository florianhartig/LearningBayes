##' Generates CR values based on current probabilities
##' @param pCR vector of length nCR, summing to 1 (?)
##' @param control must have pars nseq,steps,nCR
##' @return ...
##'   CR matrix nseq x steps. range [1/nCR,1]
GenCR <- function(pCR,control){

  ##dimensions:
  ## L vector of length nCR
  ## L2 vector of length nCR+1
  ## r vector of length nseq*steps of range [1,nseq*steps]
  ## idx vector of variable length. range of r
  ## CR vector of length nseq*steps. range [1/nCR,1]
  
  ## How many candidate points for each crossover value?
  ## TODO: verify result matches matlab
  ## MATLAB: [L] = multrnd(MCMCPar.seq * MCMCPar.steps,pCR);
  L <- as.numeric(rmultinom(1,size=control$nseq*control$steps,prob=pCR))
  L2 <- c(0,cumsum(L))
  
  ## Then select which candidate points are selected with what CR
  r <- sample(control$nseq*control$steps)

  CR <- rep(NA,control$nseq*control$steps)
  
  ## Then generate CR values for each chain
  for (zz in 1:control$nCR){
    ## Define start and end
    i.start <- L2[zz]+1
    i.end <- L2[zz+1]
    
    ## Take the appropriate elements of r
    idx <- r[i.start:i.end]
    
    ## Assign these indices control$CR(zz)
    CR[idx] <- zz/control$nCR
    
  } ## for nCR
    
  ## Now reshape CR
  CR <- matrix(CR,control$nseq,control$steps)

  return(CR)
} ## GenCR
