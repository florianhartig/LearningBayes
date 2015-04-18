##' Extract maximum likelihood parameter values
##' @param dream object
##' @param method. either a function or one of uni.mode,mean,median,sample.ml
##' @param ... arguments to window.dream
##' @return named vector of parameter values
coef.dream <- function(object,method=c("sample.ml","uni.mode","mean","median"),...)
{

  sss <- window(object,...)
  
  if (!is.function(method) && identical(match.arg(method), "sample.ml")) {
    ## TODO: make sure ppp corresponds to sss
    ppp <- window(as.mcmc(object$hist.logp),start=start(sss),thin=thin(sss))
    maxi <- which.max(ppp)
    maxchain <- col(ppp)[maxi]
    maxtime <- row(ppp)[maxi]
    return(sss[[maxchain]][maxtime,])
  }  

  if (!is.function(method)) {
    method <- switch(
                     match.arg(method),
                     "uni.mode"=maxLikCoda,
                     "mean"=function(sss) colMeans(as.matrix(sss)),
                     "median"=function(sss) apply(as.matrix(sss),2,median)
                     )
  }
  
  return(method(sss))
} ##coef.dream
