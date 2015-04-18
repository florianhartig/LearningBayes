##' Continue an existing dream MCMC set of chains
##' @param object. dream object
##' @param nsim. approximate number of function evaluations. default 1000
##' @param seed passed to set.seed before continuing
##' @return a dream object with approximately the requested number of function evaluations
## TODO. extra parameters to set in control?
## TODO: does not seem to yield stationary distribution?
## update method?
simulate.dream <- function(object,nsim=1000,seed=NULL,...){
  
  ## Generate more results from converged chains
  object$control$REPORT <- 0
  object$control$Rthres <- 0
  object$control$ndraw <- nsim
  object$control$burnin.length <- 0
  if (is.na(object$control$thin.t)) object$control$thin.t <- 10
  object$call$control <- object$control
  object$call$INIT <- function(pars,nseq) object$X[,1:object$control$ndim]

  message(sprintf("Will require %d function evaluations",nsim))

  if(!is.null(seed)) set.seed(seed)
  
  ee <- eval(object$call)
  return(ee)
  
}##simulate.dream
