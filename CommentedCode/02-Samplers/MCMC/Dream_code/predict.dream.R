##' Predict values using dream-model object
##' Predict values using function calibrated using dreamCalibrate, optionally with new data,
##' using various methods of summarising the posterior parameter and output distributions
##' @param object dream-model object
##' @param newdata. new FUN.pars list. If NULL, use object's.
##' @param method CI or a \code{\link{method}} of coef
##' @param level. Requested two-sided level of confidence. For CI method.
##' @param ... arguments to window.dream
##' @return  whatever FUN returns (either numeric, ts or list). For CI, either a matrix with upper and lower bound or list of matrices.
predict.dream_model <- function(object,newdata=NULL,
                                method="sample.ml",level=0.99, ...
                          ){

  ## Check and initialise parameters
  stopifnot(is.null(newdata) || is.list(newdata))
  stopifnot(is.function(method) || !("CI" %in% method && is.null(level)))

###
  ## Fetch function and parameters from dream-model object

  FUN <- object$FUN
  if (is.null(newdata)) newdata <- eval(object$FUN.pars)
  
  par.name <- names(formals(FUN))[1]
  
  wrap <- function(p) {
    newdata[[par.name]] <- p
    do.call(FUN,newdata)
  }
###
  
  ## Predict for desired method(s)

  if (!is.function(method) && method=="CI"){

    sss <- window(object,...)
    
    ff <- apply(as.matrix(sss),1,wrap)

    if (inherits(ff,"matrix")) return(t(apply(ff,1,quantile,c((1-level)/2,1-(1-level)/2))))
    else if (inherits(ff,"numeric")) return(quantile(ff,c((1-level)/2,1-(1-level)/2)))
    else if (inherits(ff,"list")) {
      ## Calculate CI for each series separately
      ## list is of format list[[run.number]][[series.number]]=numeric
      return(
             lapply(1:length(ff[[1]]),function(s){
               ff.s <- sapply(ff,function(x) x[[s]])
               t(apply(ff.s,1,quantile,c((1-level)/2,1-(1-level)/2)))
             })
             )
    } else stop("Unexpected output from application of FUN to matrix of parameters. FUN should return a numeric or list of numerics")
    
  } else {
    return(wrap(coef(object,method=method,...)))
  }

} ##predict.dream
