## control requires gamma, (N)
calc.rmse <- function(predicted,observed,control=list(gamma=0)){
  if (! "N" %in% names(control)) control$N <- length(observed)
  req.control <- c("N","gamma")
  if (!all(req.control %in% names(control)))
    stop(sprintf("Missing arguments to control: %s",
                 paste(req.control[!req.control %in% names(control)],collapse=", ")
                 ))

  err <- as.numeric(observed-predicted)

  SSR <- sum(abs(err)^(2/(1+control$gamma)))
  (-control$N*(1+control$gamma)/2)*0.5*SSR
}## calc.rmse

calc.weighted.rmse <- function(predicted,observed,control=list(gamma=0)){
  ##if (! "sigma" %in% names(control)) control$sigma <- sd(observed)
  req.control <- c("gamma","sigma")
  if (!all(req.control %in% names(control)))
    stop(sprintf("Missing arguments to control: %s",
                 paste(req.control[!req.control %in% names(control)],collapse=", ")
                 ))

  err <- as.numeric(observed-predicted)
  SSR <- sum(abs(err)^(2/(1+control$gamma)))
  -0.5*SSR/control$sigma^2
}

## control: Wb,Cb,gamma,measurement.sigma
calc.loglik <- function(predicted,observed,control=list(gamma=0)){

  req.control <- c("gamma","sigma")
  if (!all(req.control %in% names(control)))
    stop(sprintf("Missing arguments to control: %s",
                 paste(req.control[!req.control %in% names(control)],collapse=", ")
                 ))

  if (!all(c("Cb","Wb") %in% names(control))) control <- modifyList(control,CalcCbWb(control$gamma))
  if (! "N" %in% names(control)) control$N <- length(observed)

  err <- as.numeric(observed-predicted)

  ## Derive the log likelihood
  with(control,N*log(Wb/sigma)-Cb*(sum((abs(err/sigma))^(2/(1+gamma)))))
}## calc.loglik

## Design decisions:
##  lik.fun is log.likelihood=f(predicted,observed,control=default.list)
##   must have default, even if it is an empty list
dreamCalibrate <- function(FUN,
                           pars,
                           obs,
                           lik.fun=calc.rmse,
                           lik.control=NULL,
                           FUN.pars=list(),
                           ... ##Extra arguments to dream
                           ){
  if (is.null(lik.control)) lik.control <- eval(formals(lik.fun)$control)
  ## TODO: remove ... from here and move do.call processing from dream to here.
  wrap.lik.fun <- function(pars,...) lik.fun(FUN(pars,...),obs,lik.control)

  ## TODO: need effective and efficient way of making lik.fun etc available inside function
  ## wrap.lik.fun <- eval(parse(text=paste("function(id,pars,...){",
  ##                             "LOClik.fun <- ",paste(deparse(lik.fun),collapse="\n"),
  ##                             "LOCFUN <-",paste(deparse(FUN),collapse="\n"),
  ##                             "LOCobs <-",paste(deparse(obs),collapse="\n"),
  ##                             "LOClik.control <-",paste(deparse(lik.control),collapse="\n"),
  ##                             "LOClik.fun(LOCFUN(id,pars,...),LOCobs,LOClik.control)}",
  ##                             collapse="\n",sep="\n"
  ##                             )))
  dd <- dream(FUN=wrap.lik.fun,
              pars=pars,
              func.type="logposterior.density",
              FUN.pars=FUN.pars,
              ... ##INIT,control
              )

  dd$call <- match.call()
  dd$FUN <- FUN
  dd$FUN.pars <- FUN.pars
  dd$lik.fun <- function(p) do.call(wrap.lik.fun,modifyList(FUN.pars,list(pars=p)))
  class(dd) <- c("dream_model",class(dd))
  dd
} ##dreamCalibrate
