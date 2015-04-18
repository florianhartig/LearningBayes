##' Determine which sequences to evolve with what DE strategy
##' @param control list with elements DEpairs, nseq
##' @return DEversion 
##'   vector of length control$nseq. range [1,DEpairs]
##
##  ? maximises cumulative sum of 1/DEpairs < random number
##  TODO: alternative implementation?
DEStrategy<-function(control){

  ## dimensions:
  ##  p.pair vector length DEpairs+1. range [0,1]
  ##  Z vector length nseq. range [0,1]
  ##  qq. iter 1:nseq
  
  ## Determine probability of selecting a given number of pairs
  p.pair <- (1/control$DEpairs)*rep(1,control$DEpairs)
  p.pair <- c(0,cumsum(p.pair))
  
  ## Generate a random number between 0 and 1
  Z <- runif(control$nseq)
  
  ## Select number of pairs
  DEversion<-rep(NA,control$nseq)
  for (qq in 1:control$nseq){
    DEversion[qq]<-tail(which(p.pair<Z[qq]),1)
  }
  return(DEversion)
}#DEStrategy
