## Copyright (c) 2008, Los Alamos National Security, LLC
## All rights reserved.
##
## Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.
## Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is
## operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.
## Government has rights to use, reproduce, and distribute this software.
##
## NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR
## IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to
## produce derivative works, such modified software should be clearly marked, so as not to
## confuse it with the version available from LANL.
##
## Additionally, redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## * Redistributions of source code must retain the above copyright notice, this list of
##   conditions and the following disclaimer.
## * Redistributions in binary form must reproduce the above copyright notice, this list of
##   conditions and the following disclaimer in the documentation and/or other materials
##   provided with the distribution.
## * Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL
##   the U.S. Government, nor the names of its contributors may be used to endorse or promote
##   products derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
## OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
## ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
## HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
## EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## MATLAB code written by Jasper A. Vrugt, Center for NonLinear Studies (CNLS)
##
## Written by Jasper A. Vrugt: vrugt@lanl.gov


################################################################################################
## This R code has been converted and modified from the original MATLAB code by               ##
## Joseph Guillaume and Felix Andrews, 2010.                                                  ##
################################################################################################


dreamDefaults <- function()
    list(
         nCR = 3,                ## Crossover values used to generate proposals (geometric series)
         gamma = 0,              ## Kurtosis parameter Bayesian Inference Scheme
         steps = 10,             ## Number of steps in sem
         eps = 5e-2,             ## Random error for ergodicity
         outlierTest = 'IQR_test', ## What kind of test to detect outlier chains?
         pCR.Update = TRUE,      ## Adaptive tuning of crossover values
         boundHandling = 'reflect', ## Boundary handling: "reflect", "bound", "fold", "none"
         burnin.length=0.1, ## Proportion of iterations considered to be burnin. 0 to turn off.
### Termination criteria. TODO: are the 2nd two valid, given that ndraw is used in adaptive pcr
         ndraw = 1e5,            ## maximum number of function evaluations
         maxtime = Inf,           ## maximum duration of optimization in seconds
         Rthres=1.01,            ## R value at which to stop. Vrugt suggests 1.2
### Thinning
         thin.t=NA,            ## parameter for reduced sample collection
### Efficiency improvements
         REPORT = 1000,            ## approximate number of function evaluations between reports. >0. 0=none
         parallel = "none",   ##packages to use for parallel in order of preference: multicore,snow,foreach
### Parameters with auto-set values
         ndim=NA,			 ## number of parameters (automatically set from length of pars)
         DEpairs = NA,          ## Number of DEpairs. defaults to max val floor((nseq-1)/2)
         nseq = NA,              ## Number of Markov Chains / sequences (defaults to N)
         Cb=NA,Wb=NA
         ## Currently unused parameters
##         trace = 0,              ## level of user feedback
         )

library(coda)



##' @param FUN model function with first argument a vector of parameter values of length ndim
##' @param func.type type of value FUN returns.
##'  one of: posterior.density, logposterior.density,calc.loglik,calc.rmse,calc.weighted.rmse
##' @param pars a list of variable ranges
##' @param INIT f(pars,nseq,...) returns nseq x ndim matrix of initial parameter values
##' @param control
##' @param measurement list. must be included unless func.type=posterior.density or logposterior.density is selected
##'  for calc.rmse: must have element data

##' @return ...
##'   TODO
##'   X converged nseq points in parameter space. matrix nseq x ndim
##'   Sequences mcmc.list object. nseq mcmc elements of ndim variables
##'   Reduced.Seq mcmc.list object. nseq mcmc elements of ndim variables
##'   AR acceptance rate for each draw. matrix max.counter x 2
##'   outlier vector of variable length
##'   R.stat Gelman.Diag statistic for each variable at each step. matrix max.counter/steps x 1+ndim
##'   CR. Probability of crossover. matrix max.counter/steps x 1+length(pCR)
##'

##' Terminates either when control$ndraw or control$maxtime is reached
##'
##' MATLAB function:
##' function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option)

dream <- function(FUN, func.type,pars,
                  FUN.pars=list(),
                  INIT = LHSInit,
                  INIT.pars=list(),
                  control = list(),
                  measurement=NULL
                  )
{


  ## dimensions
  ##  hist.logp matrix. ndraw/nseq x nseq. length nearly ndraw.
  ##    TODO: removed counter.fun.evals for simplicity. should have been kept?
  ##  CR nseq x steps
  ##  pCR length nCR or scalar
  ##  lCR length nCR or scalar
  ##  Table.JumpRate ndim x DEpairs. range (0,~1.683]
  ##  delta.tot vector of length nCR

  ## Sequences. array max.counter*1.125 x ndim+2 x nseq
  ## Reduced.Seq array max.counter*1.125 x ndim+2 x nseq

  ##  counter.outloop [2,ndraw/nseq]. count number of outside loops
  ## counter.fun.evals [nseq,ndraw(+steps*nseq)],
  ## counter [2,ndraw/nseq] . number of generations - iterations of inner loop

############################
  ## Process parameters

  ## Check validity of parameters
  if (is.character(FUN))  FUN <- get(FUN, mode = "function")
  stopifnot(is.function(FUN))
  stopifnot(is.list(pars))
  stopifnot(length(pars) > 0)
  pars <- lapply(pars, function(x) if (is.list(x)) x else list(x))
  if (is.null(names(pars))){
    pad.length <- nchar(as.character(length(pars)))
    names(pars) <- sprintf(paste("p%0",pad.length,"d",sep=""),1:length(pars))
  }
  stopifnot(is.list(control))
  stopifnot(func.type %in% c("calc.rmse","calc.loglik","calc.weighted.rmse","posterior.density","logposterior.density"))
  stopifnot(!is.null(measurement) || func.type %in% c("posterior.density","logposterior.density"))
  stopifnot(!func.type %in% c("calc.rmse","calc.loglik","calc.weighted.rmse") || "data" %in% names(measurement))

  ## Check INIT and FUN have required extra parameters in INIT.pars & FUN.pars
  req.args.init <- names(formals(INIT))
  if(!all(req.args.init %in% c("pars","nseq",names(INIT.pars))))
    stop(paste(c("INIT Missing extra arguments:",
                 req.args.init[!req.args.init %in% c("pars","nseq",names(INIT.pars))]),
               sep=" ",collapse=" "))

  req.args.FUN <- names(formals(FUN))
  ## if (length(req.args.FUN)<length(FUN.pars)+1) stop("Some FUN.pars are not required by FUN")
  ## if (length(req.args.FUN)>1){
  ##   req.args.FUN <- req.args.FUN[2:length(req.args.FUN)] ##optional pars only
  ##   if(!all(req.args.FUN %in% names(FUN.pars))) stop(paste(c("FUN Missing extra arguments:",
  ##                                                            req.args.FUN[!req.args.FUN %in% c(names(FUN.pars))]),
  ##                                                          collapse=" "))
  ## }

  ## Update default settings with supplied settings

  control <- modifyList(dreamDefaults(), control)
  isValid <- names(control) %in% names(dreamDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))

  if (!is.null(measurement)){
    if (! "sigma" %in% names(measurement)) measurement$sigma <- sd(measurement$data)
    if (! "N" %in% names(measurement)) measurement$N <- length(measurement$data)
  }

  ## Set automatically determined values
  control$ndim<-length(pars)
  if (is.na(control$nseq)) control$nseq <- control$ndim
  if (is.na(control$DEpairs)) control$DEpairs <- floor((control$nseq-1)/2)

  ## Correct to match nseq
  control$REPORT <- (control$REPORT%/%control$nseq) * control$nseq

  if (control$burnin.length<1) control$burnin.length <- control$burnin.length*control$ndraw
  if (identical(tolower(control$outlierTest),'none')) control$burnin.length <- 0

  ##Choice of parallel backend
  if (!control$parallel %in% c("none","snow.chains")){
    parallel <- "none"
    for (p in control$parallel) {
      if (require(p,character.only=TRUE)) {
        parallel <- p
        break
      }
    }
    if (parallel=="none") warning(sprintf("Requested parallel backends not available (%s)",paste(control$parallel,collapse=",")))
    control$parallel <- parallel
  }

  ## Check validity of settings
  if (control$DEpairs==0) stop("control$DEpairs set to 0. Increase nseq?")
  stopifnot(control$DEpairs<=(control$nseq-1)/2) ## Requirement of offde
  stopifnot(control$boundHandling %in% c("reflect", "bound", "fold", "none"))
  if (control$boundHandling == 'none') warning("No bound handling in use, parameters may cause errors elsewhere")
  stopifnot(control$REPORT>=0)
  if (control$parallel=="snow.chains"){
      library(snow)
      if (func.type != "logposterior.density") stop("control$parallel=snow.chains only supports func.type=logposterior.density")
      if (!is.null(control$measurement)) stop("control$measurement is ignored for control$parallel=snow.chains")
      if (!identical(FUN.pars,list())) stop("FUN.pars is ignored for control$parallel=snow.chains")

  }

############################
  ## Initialize variables

  NDIM <- control$ndim
  NCR <- control$nCR
  NSEQ <- control$nseq

  ## Counters
  counter.fun.evals <- NSEQ
  counter <- 2
  counter.outloop <- 2
  counter.thin <- 1

  ## Max number of times through loops
  max.counter <- ceiling((control$ndraw+control$steps*NSEQ)/NSEQ)+1
  max.counter.outloop <- ceiling((control$ndraw+control$steps*NSEQ)/NSEQ/control$steps)+1

  if (!is.na(control$thin.t) && floor(max.counter/control$thin.t)<2) stop(sprintf("Thin parameter should be much smaller than number of iterations (currently %d,%d)",control$thin.t,max.counter))

  ## Calculate the parameters in the exponential power density function of Box and Tiao (1973)
  if (!func.type %in% c("posterior.density","logposterior.density")){
    cbwb <- CalcCbWb(control$gamma)
    if (is.na(control$Cb))  control$Cb <- cbwb$Cb
    if (is.na(control$Wb))  control$Wb <- cbwb$Wb
  }
  
  ## Generate the Table with JumpRates (dependent on number of dimensions and number of pairs)
  Table.JumpRate<-matrix(NA,NDIM,control$DEpairs)
  for (zz in 1:control$DEpairs) Table.JumpRate[,zz] <- 2.38/sqrt(2*zz*1:NDIM)

  ## Initialize the array that contains the history of the log_density of each chain
  hist.logp <- matrix(NA_real_,max.counter,NSEQ)
  real.hist.logp <- matrix(NA_real_,max.counter,NSEQ)

  if (control$pCR.Update){
    ## Calculate multinomial probabilities of each of the nCR CR values
    pCR <- rep(1/NCR,NCR)

    ## Calculate the actual CR values based on p
    CR <- GenCR(pCR,control)
    lCR <- rep(0,NCR)
  } else {
    pCR <- 1/NCR
    ## Define
    CR <- matrix(pCR,NSEQ,control$steps)
    lCR <- NSEQ*control$steps
  } ##pCR.Update

  end.burnin <- control$burnin.length



############################
  ## Initialise output object

  obj <- list()
  class(obj) <- c("dream", class(obj))
  obj$call <- match.call()
  obj$control <- control

  obj$in.burnin <- TRUE

  obj$EXITFLAG <- NA
  obj$EXITMSG <- NULL

  ## counter.fun.evals + AR at each step
  obj$AR<-matrix(NA,max.counter,2)
  obj$AR[1,2]<-NSEQ-1 ##Number if only one rejected
  colnames(obj$AR) <- c("fun.evals", "AR")

  ##counter.fun.evals + R statistic for each variable at each step
  ## TODO: now using counter.report
  obj$R.stat<-matrix(NA,max.counter.outloop,NDIM+1)
  ##  n<10 matlab: -2 * ones(1,MCMCPar.n);
  obj$R.stat[1,] <- c(counter.fun.evals,rep(-2,NDIM))
  colnames(obj$R.stat) <- c("fun.evals", names(pars))

  ##counter.fun.evals + pCR for each CR
  obj$CR <- matrix(NA,max.counter.outloop,length(pCR)+1)
  colnames(obj$CR) <- c("fun.evals",paste("CR",1:length(pCR),sep=""))

  obj$outlier<-NULL

  Sequences <- array(NA_real_, c(max.counter,NDIM+2,NSEQ))
  colnames(Sequences) <- c(names(pars), "p", "logp")
  ## Sequences[1,] <- sapply(pars, mean) ## TODO: include?

  ## Check whether will save a reduced sample
  if (!is.na(control$thin.t)){
    counter.redseq <- 0
    Reduced.Seq <- array(NA_real_,c(ceiling(max.counter/control$thin.t),NDIM+2,NSEQ))
  } else Reduced.Seq <- NULL

############################

  ## Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
  control$steps<-control$steps-1

  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0
  counter.report <- 1

################################

  ## Step 1: Sample s points in the parameter space

  x <- do.call(INIT,modifyList(INIT.pars,list(pars=pars,nseq=NSEQ)))

  ## Test that FUN returns numeric
if (control$parallel!="snow.chains"){ ## TODO: something equivalent with snow.chains?
  test.pars <- FUN.pars
  test.pars[[names(formals(FUN))[1]]] <- x[1,]
  modpred <- do.call(FUN,test.pars)
  if (!is.numeric(modpred))
      stop(sprintf("Result of FUN should be of class numeric, not %s", class(modpred)))
}
  ## make each element of pars a list and extract lower / upper
  lower <- sapply(pars, function(x) min(x[[1]]))
  upper <- sapply(pars, function(x) max(x[[1]]))

  ##Step 2: Calculate posterior density associated with each value in x
  tmp <- CompDensity(pars = x, control = control, FUN = FUN, func.type = func.type,
                     measurement = measurement, FUN.pars = FUN.pars)

  ##Save the initial population, density and log density in one list X
  X <- cbind(x = x, p = tmp$p, logp = tmp$logp)
  colnames(X) <- c(names(pars), "p", "logp")

  ##Initialise the sequences
  for (qq in 1:NSEQ){
    Sequences[1,,qq] <- X[qq,]
  }

  ##Save pCR in memory and initialize delta.tot
  obj$CR[1,] <- c(counter.fun.evals,pCR)
  delta.tot <- rep(0,NCR)

  ##Save history log density of individual chains
  hist.logp[1,] <- X[,"logp"]
  real.hist.logp[1,] <- X[,"logp"]

################################
  ##Start iteration
  ## max times (ndraw+steps*NSEQ)/(NSEQ*steps)
  while (counter.fun.evals < control$ndraw) {

    ## max times ceiling(ndraw+steps*NSEQ)/NSEQ
    for (gen.number in 1:control$steps) {

      ## Initialize DR properties
      counter.thin <- counter.thin + 1

      ## Define the current locations and associated posterior densities
      x.old <- X[,1:NDIM,drop=FALSE]
      p.old <- X[,NDIM+1]
      logp.old <- X[,NDIM+2]

      ## Now generate candidate in each sequence using current point and members of X
      ## Table.JumpRate appears to match matlab version
      tmp <- offde(x.old, control = control,
                   CR=CR[,gen.number],
                   lower = lower, upper = upper,
                   Table.JumpRate=Table.JumpRate)
      x.new <- tmp$x.new
      stopifnot(!identical(x.new,x.old))
      CR[,gen.number] <- tmp$CR

      ## Now compute the likelihood of the new points
      tmp <- CompDensity(pars = x.new, control = control, FUN = FUN, func.type = func.type,
                         measurement = measurement, FUN.pars = FUN.pars)
      p.new <- tmp$p
      logp.new <- tmp$logp

      ## Now apply the acceptance/rejectance rule for the chain itself
      tmp <- metrop(x.new,p.new,logp.new,
                    x.old,p.old,logp.old,
                    func.type,control,
                    measurement
                    )
      X <- tmp$newgen ##Update X using current members of sequences. N.B. used to be separate step after Delayed Rejection
      alpha12 <- tmp$alpha
      accept <- tmp$accept
      ## stopifnot(any(accept)) #Unlikely, but possible

      ## NOTE: original MATLAB code had option for DR Delayed Rejection here)
      ## accept2,ItExtra not required

      ## Update location in sequence and update the locations of the Sequences with the current locations
      Sequences[counter,,] <- t(X)

      ## Check reduced sample collection
      if (!is.na(control$thin.t) && counter.thin == control$thin.t){
        ## Update iloc_2 and counter.thin
        counter.redseq <- counter.redseq+1
        counter.thin <- 0
        ## Reduced sample collection
        Reduced.Seq[counter.redseq,,] <- t(X)
      }

      if (control$pCR.Update) {
        ## Calculate the standard deviation of each dimension of X
        ## TODO: matlab syntax is unclear - seems to be columnwise
        ## element-wise: sd(c(X[,1:NDIM]))
        r <- apply(X[,1:NDIM,drop=FALSE],2,sd)
        ## Compute the Euclidean distance between new X and old X
        delta.normX <- rowSums(((x.old-X[,1:NDIM,drop=FALSE])/r)^2)
        ## Use this information to update sum_p2 to update N_CR
        delta.tot <- CalcDelta(NCR,delta.tot,delta.normX,CR[,gen.number])
      }

      ## Update hist.logp
      hist.logp[counter,] <- X[,NDIM+2]
      real.hist.logp[counter,] <- X[,NDIM+2]

      ## Save Acceptance Rate
      obj$AR[counter,] <- c(counter.fun.evals,100 * sum(accept) / NSEQ)

      ## CompDensity executes function NSEQ times per loop
      counter.fun.evals <- counter.fun.evals + NSEQ
      counter <- counter + 1
    } ##for gen.number steps

    ## ---------------------------------------------------------------------

    ## Store Important Diagnostic information -- Probability of individual crossover values
    obj$CR[counter.outloop, ] <- c(counter.fun.evals,pCR)

    ## Do this to get rounded iteration numbers
    if (counter.outloop == 2) control$steps <- control$steps + 1

    if (control$burnin.length!=0) outliers <- RemOutlierChains(X,hist.logp[1:(counter-1),],control)

    if (counter.fun.evals <= end.burnin) {
      ## Check whether to update individual pCR values
      if (control$pCR.Update) {
        ## Update pCR values
        tmp <- AdaptpCR(CR, delta.tot, lCR, control)
        pCR <- tmp$pCR
        lCR <- tmp$lCR
      }

      ## Change any outlier chains to current best value of X
      ## TODO: matlab code didn't match paper. Outlier removal should be within burnin period

      ## Loop over each outlier chain (if length>0)
      for (out.id in outliers$chain.id){
        ## Draw random other chain -- cannot be the same as current chain
        r.idx <- which.max(outliers$mean.hist.logp)
        ## Added -- update hist_logp -- chain will not be considered as an outlier chain then
        hist.logp[1:(counter-1),out.id] <- hist.logp[1:(counter-1),r.idx]
        real.hist.logp[(counter-1),out.id] <- real.hist.logp[(counter-1),r.idx]
        ## Jump outlier chain to r_idx -- Sequences, X
        Sequences[(counter-1),1:(NDIM+2),out.id] <- X[r.idx,]
        X[out.id,1:(NDIM+2)] <- X[r.idx,]
        ## Add to outlier tracker
        obj$outlier <- rbind(obj$outlier,c(counter.fun.evals,out.id))
      } ##for remove outliers

    } else if (control$burnin.length!=0){
      obj$in.burnin <- FALSE
      if (length(outliers)>0){
        warning("Outliers detected outside burn-in period, returning to burn in")
        obj$in.burnin <- TRUE
        end.burnin <- counter.fun.evals+control$burnin.length
      }
    }   ##in burn in.


    if (control$pCR.Update) {
      ## Generate CR values based on current pCR values
      CR <- GenCR(pCR, control)
    } else {
      CR <- matrix(pCR, nrow = NSEQ, ncol = control$steps)
    }

    ## Calculate Gelman and Rubin convergence diagnostic
    ## Compute the R-statistic using 50% burn-in from Sequences
    ## TODO: alternatively, convert matlab implementation
    if (control$REPORT>0 && counter.fun.evals %% control$REPORT==0) {

      counter.report <- counter.report+1

      try({
        obj$R.stat[counter.report,] <-
          c(counter.fun.evals,
            gelman.diag(
                        as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[1:(counter-1),1:NDIM,i]))),
                        autoburnin=TRUE)$psrf[,1])
        if (counter.report == 2){
          message("R.stats:")
          message(format(colnames(obj$R.stat), width = 10,justify="right"))
        }
        message(format(obj$R.stat[counter.report,], width = 10, digits = 4))
      })

      if (all(!is.na(obj$R.stat[counter.report,])) &&
          all(obj$R.stat[counter.report,-1]<control$Rthres)) {
        obj$EXITMSG <- 'Convergence criteria reached'
        break
        ## obj$EXITFLAG <- 3
      }

    } ##counter.report

    ## Update the counter.outloop
    counter.outloop = counter.outloop + 1

    ## break if maximum time exceeded
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      obj$EXITMSG <- 'Exceeded maximum time.'
      obj$EXITFLAG <- 2
      break
    }

  } ##while

  toc <- as.numeric(Sys.time()) - tic

  if (counter.fun.evals>= control$ndraw){
    obj$EXITMSG <- "Maximum function evaluations reached"
   ## obj$EXITFLAG <- 4
  }



  ## Trim outputs to collected data - remove extra rows
  ## Convert sequences to mcmc objects
  ## TODO: remove all prior refs to Reduced.Seq - now not needed, given window.mcmc is used
  Sequences <- Sequences[1:(counter-1),,]
  obj$Sequences <- as.mcmc.list(lapply(1:NSEQ,function(i) as.mcmc(Sequences[,1:NDIM,i])))
  if(!is.na(control$thin.t)){
    obj$Sequences <- window(obj$Sequences,thin=control$thin.t)
  }

  ## TODO: make these 'ts' objects and sync with Reduced.Seq by thinning
  ##  Would it be better to keep all data, and sync when needed by matching start, end,thin?
  ##   See coef.dream for eg.
  obj$X <- X
  obj$R.stat <- obj$R.stat[1:counter.report,,drop=FALSE]
  obj$hist.logp <- real.hist.logp[1:(counter-1),,drop=FALSE]
  obj$AR <- obj$AR[1:(counter-1),,drop=FALSE]
  obj$CR <- obj$CR[1:(counter.outloop-1),,drop=FALSE]

  ## store number of iterations
  obj$iterations <- counter.outloop
  ## store number of function evaluations
  obj$fun.evals <- counter.fun.evals
  ## store the amount of time taken
  obj$time <- toc

  obj
} ##dream


