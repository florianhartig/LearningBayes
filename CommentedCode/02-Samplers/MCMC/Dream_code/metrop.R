##' Metropolis rule for acceptance or rejection
##'
##' @param x matrix nseq x ndim
##' @param p.x vector of length nseq
##' @param logp.x vector of length nseq
##' @param x.old matrix nseq x ndim
##' @param p.old vector of length nseq
##' @param logp.old vector of length nseq
##' @param func.type Type of function output. One of:
##'   posterior.density, logposterior.density,
##'   calc.loglik. requires optional parameter measurement with elements data & sigma
##'   calc.rmse, calc.weighted.rmse.  requires measurement$data
##' @param control list with elements:
##' @param measurement. needs N and sigma
##'   gamma
##'   metrop.opt range [1,5] 
##' @return ... list with elements
##'   newgen matrix nseq x ndim+2 (same as X in dream.R)
##'   alpha scalar probability of acceptance. range [0,1]
##'   accept vector indicating whether each sequences was accepted. length nseq
metrop<-function(x,p.x,logp.x,
                 x.old,p.old,logp.old,
                 func.type,control,
                 measurement=NULL
                 ){

  stopifnot(!is.null(measurement) || func.type%in% c("posterior.density","logposterior.density"))
  stopifnot(!any(is.na(p.x)))
  
  ## dimensions:
  ##  nr.chains scalar. should be = control$nseq
  ##  Z vector of length nseq range [0,1]
  ##  idx vector length [0,nseq] range [1,nseq]
  
  ## Calculate the number of Chains
  ## TODO: redundant because always equal to control$nseq?
  nr.chains <- nrow(x)
  
  ## First set newgen to the old positions in X
  newgen <- cbind(x.old,p.old,logp.old)
  
  ## And initialize accept with false
  accept <- rep(FALSE,nr.chains)
  
  switch(func.type,
         posterior.density = {
           alpha <- pmin(p.x/p.old,1)
         },
         calc.loglik = { ## Lnp probability evaluation
           alpha <- pmin(exp(p.x-p.old),1)
         },
         calc.rmse = { ## SSE probability evaluation
           alpha <- pmin((p.x/p.old)^(-measurement$N*(1+control$gamma)/2),1)
         },
         logposterior.density = { ## Lnp probability evaluation
           alpha <- pmin(exp(p.x-p.old),1)
         },
         calc.weighted.rmse = { ## Similar to 3 but now weighted with Measurement.Sigma
           ## signs are different because we write -SSR
           alpha <- pmin(exp(-0.5*(-p.x + p.old)/measurement$sigma^2),1);
         },
         stop("Unrecognised value of func.type")
         )
  
  ## Generate random numbers
  Z <- runif(nr.chains)
  ## Find which alpha's are greater than Z
  idx <- which(alpha>Z)

  ##stopifnot(length(idx)>0) ##Unlikely, but possible
  
  ## And update these chains
  newgen[idx,] <- c(x[idx,],p.x[idx],logp.x[idx])
         
  ## And indicate that these chains have been accepted
  accept[idx] <- TRUE

  return(list(newgen=newgen,alpha=alpha,accept=accept))
} ##metrop
