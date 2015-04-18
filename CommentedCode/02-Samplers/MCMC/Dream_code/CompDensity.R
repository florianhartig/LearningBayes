##' Computes the density of each pars value
##'
##' @param pars matrix nseq x ndim
##' @param control. list containing gamma,Wb,Cb,nseq
##' @param FUN - the model to run
##'   R function with first argument a vector of length ndim.
##'   returns a scalar or vector corresponding to one of the options below.
##' @param func.type Type of function output. One of:
##'   posterior.density, logposterior.density,
##'   calc.loglik. requires optional parameter measurement with elements data & sigma
##'   calc.rmse, calc.weighted.rmse.  requires measurement$data
##' @param measurement list containing TODO: not sure
##'   data: vector of observations corresponding to model output
##'   sigma: scalar
##' @param ... additional arguments to FUN
##' @return list with elements
##'   p vector of length nseq
##'   logp vector of length nseq
##
## TODO: p may be erroneously equal to logp?
## TODO: more appropriate naming of options?
## TODO: allow shortenings of option?
CompDensity <- function(pars,control,FUN,func.type,
                        measurement=NULL,FUN.pars=NULL){

  ## Should be guaranteed by dream
  ## stopifnot(!is.null(measurement) || func.type%in% c("posterior.density","logposterior.density"))

  stopifnot(!any(is.na(pars)))

  if (control$parallel=="snow.chains"){
      ## Custom code for John Joseph 6/8/2011
      ## FUN should be logp=f(cluster instance identifier, list of pars)
      ## func.type,measurement,FUN.pars should all be NULL
      logp <- as.numeric(clusterApply(cl=cl,x=1:nrow(pars),fun=FUN,pars=pars))
      if (any(is.na(logp))) {
          stop("likelihood function produced invalid probabilities (NA/NaN)")
      }
      ##stopifnot(!any(is.na(logp))) ##Not used anyway
      return(list(p=logp,logp=logp))
  }

  ## dimensions:
  ##  i. iter 1:nseq
  ##  modpred. scalar or vector commensurate to measurement$data
  ##  err. vector of same length as modpred
  ##  SSR scalar
  ## temp. list of length nseq with elements of length 2: p and logp

  do.calc <- function (pp,control,MFUN,func.type,measurement,FUN.pars){
    ## Call model to generate simulated data
    FUN.pars[[names(formals(FUN))[1]]] <- pp
    modpred <- do.call(MFUN,FUN.pars)
    ##stopifnot(inherits(modpred,"numeric"))

    switch(func.type,
           ## Model directly computes posterior density
           posterior.density={
             p <- modpred
             if (any(modpred<0)) stop("Posterior density returned by FUN should be positive. Otherwise use logposterior.density?")
             logp <- log(modpred)
           },
           ## Model computes output simulation
           calc.loglik={
             err <- as.numeric(measurement$data-modpred)

             ## Derive the log likelihood
             logp <- measurement$N*log(control$Wb/measurement$sigma)-
               control$Cb*(sum((abs(err/measurement$sigma))^(2/(1+control$gamma))))
             ## And retain in memory
             p <- logp
           },
           ## Model computes output simulation
           ## TODO: may need as.numeric
           calc.rmse={
             err <- as.numeric(measurement$data-modpred)
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p <- -SSR
             logp <- -0.5*SSR

           },
           ## Model directly computes log posterior density
           logposterior.density={
             p <- modpred
             logp <- modpred
             stopifnot(all(logp<=0))
           },
           ## Similar as 3, but now weights with the Measurement Sigma
           ## TODO: identical to rmse because difference is in metrop
           calc.weighted.rmse={
             ## Define the error
             err <- as.numeric(measurement$data-modpred)
             ## Derive the sum of squared error
             SSR <- sum(abs(err)^(2/(1+control$gamma)))
             ## And retain in memory
             p <- -SSR
             logp <- -0.5*SSR
           }) ##switch
    c(p,logp)
  }

  pars <- lapply(1:nrow(pars),function(x) pars[x,])
  switch(control$parallel,
         "multicore"={
           temp <- mclapply(pars,do.calc,control=control,MFUN=FUN,func.type=func.type,
                            measurement=measurement,FUN.pars=FUN.pars,
                            mc.preschedule=FALSE)
         },
         "foreach"={
           ## foreach(pp=iter(pars,by="row")) %dopar%
           temp <- foreach(pp=pars) %dopar% {do.calc(pp,control=control,MFUN=FUN,func.type=func.type,
                             measurement=measurement,FUN.pars=FUN.pars)}
         },
         "snow"={
           temp <- clusterApplyLB(cl=cl,x=pars,fun=do.calc,
                                  control=control,MFUN=FUN,func.type=func.type,
                                  measurement=measurement,FUN.pars=FUN.pars)
         },
         temp <- lapply(pars,FUN=do.calc,control=control,MFUN=FUN,func.type=func.type,
                        measurement=measurement,FUN.pars=FUN.pars)
         )##switch parallel

  p <- sapply(temp,function(x) x[1])
  logp <- sapply(temp,function(x) x[2])

  if (!is.numeric(p)) {
    print(p)
    stop(sprintf("Expected class numeric, got class %s. Error with multicore? Set control$parallel='none' to not use parallelisation",class(p)))
  }
  if (any(is.na(p))) {
      stop("likelihood function produced invalid probabilities (NA/NaN)")
  }
  ##stopifnot(!any(is.na(logp))) ##Not used anyway
  return(list(p=p,logp=logp))
} ##CompDensity
