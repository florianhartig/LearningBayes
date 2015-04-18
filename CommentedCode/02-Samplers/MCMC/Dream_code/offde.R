##' Generate offspring using METROPOLIS HASTINGS monte-carlo markov chain
##'
##' @param x.old matrix nseq x ndim. parameter-wise range [xmin,xmax]
##' @param control needs nseq,ndim,eps
##' @param CR vector. length nseq. range [1/nCR,1]
##' @param lower see handleBounds
##' @param upper see handleBounds
##' @param Table.JumpRate. defined in dream. matrix ndim x DEpairs. range (0,~1.683]
##' @return list with elments
##   x.new, CR. same dim and range as input
##
## depends:
##  handleBounds.R
offde<-function(x.old,control,CR,
                lower,upper,
		Table.JumpRate #TODO: may be better in another accessible scope rather than as a parameter
                ){

  nseq <- control$nseq
  ndim <- control$ndim

  ## dimensions
  ##  eps. matrix nseq x ndim
  ##  DEversion. vector of length nseq. range [0,DEpairs]
  ##  tt. matrix nseq-1 x nseq. range [1,nseq-1]
  ##  D. matrix nseq x ndim. range [0,1]
  ##  noise.x. matrix nseq x ndim
  ##  delta.x. matrix nseq x ndim
  ##  qq. iter. range [1,nseq]
  ##  ii. vector. length nseq-1. range [1,nseq]
  ##  rr. vector. length [0,2*DEpairs]. range [1,nseq]
  ##  i. vector. length [1,ndim]. range [1,ndim]
  ##  JumpRate. scalar. range (0,~1.683]
  ##  delta.  vector. length ndim
  ##  R. matrix ndim x ndim
 
  #Generate ergodicity term
  eps <- 1e-6 * randn(nseq,ndim)

  #Not a delayed rejection step -> generate proposal with DE
  ## Determine which sequences to evolve with what DE strategy
  DEversion <- DEStrategy(control)
  
  ## Generate series of permutations of chains
  ## Get indices of column-wise sort
  # MATLAB: tt<-sort(rand(nseq-1,nseq))
  tt <- apply(rand(nseq-1,nseq),2,order)
  
  ## Generate uniform random numbers for each chain to determine which dimension to update
  D <-rand(nseq,ndim)

  ## Ergodicity for each individual chain
  noise.x<-control$eps * (2*rand(nseq,ndim)-1)
    
  ## Initialize the delta update
  delta.x<-matrix(0,nseq,ndim)
  
  ## Each chain evolves using information from other chains to create offspring
  ## TODO: can this loop be parallelised?
  ##  is an expensive part of dream. 36.6% of time
  ##  but depends on state of other loops
  for (qq in 1:nseq){
    
    ## Define ii and remove current member as an option
    ii <- (1:nseq)[-qq]
    
    ## randomly select two members of ii
    rr <- ii[tt[1:(2*DEversion[qq]),qq]]
    
    ## --- WHICH DIMENSIONS TO UPDATE? DO SOMETHING WITH CROSSOVER ----
    i <- which(D[qq,]>(1-CR[qq]))

    ## Update at least one dimension
    if (length(i)==0) i <- sample(ndim,1)
    
    ## ----------------------------------------------------------------
    ## Determine the number of dimensions that are going to be updated
    NrDim <- length(i)
    
    ## Determine the associated JumpRate and compute the jump
    if (runif(1)<4/5){
      ## Lookup Table
      JumpRate <- Table.JumpRate[NrDim,DEversion[qq]]
      
      ## Produce the difference of the pairs used for population evolution
      delta <- colSums(x.old[rr[1:DEversion[qq]],,drop=FALSE]-x.old[rr[(DEversion[qq]+1):(2*DEversion[qq])],,drop=FALSE])
      
      ## Then fill update the dimension
      delta.x[qq,i] <- (1+noise.x[qq,i])*JumpRate*delta[i]

    } else {
      ## Set the JumpRate to 1 and overwrite CR and DEversion
      JumpRate <- 1
      CR[qq] <- -1
      
      ## Compute delta from one pair
      delta <- x.old[rr[1],]-x.old[rr[2],]
                       
      ## Now jumprate to facilitate jumping from one mode to the other in all dimensions
      delta.x[qq,] <- JumpRate*delta
    }##runif
    
    ## Avoid that jump = 0 and xnew is similar to xold
    if (sum(delta.x[qq,]^2)==0){
      ## Compute the Cholesky Decomposition of x.old
      R <- (2.38/sqrt(nseq))*chol(cov(x.old)+1e-5*diag(ndim))
        
      ## Generate jump using multinormal distribution
      delta.x[qq,] <- c( rnorm(ndim) %*% R )
    }

  }#for qq
  
  ## TODO?:
  ## If delayed rejection step --> generate proposal with DR
  ## Loop over all chains -- all dimensions are updated
  ## Generate a new proposal distance using standard procedure
  
  ## Update x.old with delta_x and eps;
  x.new <- x.old+delta.x+eps
  x.new <- handleBounds(x.new,lower,upper,control$boundHandling)

  stopifnot(!any(is.na(x.new)))
  stopifnot(!any(is.na(CR)))
  return(list(x.new=x.new,CR=CR))
}#function offde
