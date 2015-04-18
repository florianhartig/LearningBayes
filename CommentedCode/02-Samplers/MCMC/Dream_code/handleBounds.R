
##' Checks the bounds of the parameters
##' @param x usually matrix nseq x ndim. vector interpreted as 1 x ndim matrix
##' @param lower vector of length ndim
##' @param upper vector of length ndim
##' @param bound.handling one of reflect, bound, fold, none
##' @return x matrix nseq x ndim now with parameter-wise range [xmin,xmax]
handleBounds <- function(x, lower, upper, bound.handling)
{
  if (is.vector(x)) x<-t(x)
  stopifnot(is.matrix(x))
  
  ## Iterate through parameters
  ## Modify points that are below or above bounds
  for (p in 1:ncol(x)){
    too.low<-which(x[,p]<lower[p])
    too.high<-which(x[,p]>upper[p])
    switch(bound.handling,
           reflect = {
             ## TODO: may still violate bounds if x>2*upper
             x[too.low,p] <- 2*lower[p]-x[too.low,p]
             x[too.high,p] <- 2*upper[p]-x[too.high,p]
           },
           bound = {
             x[too.low,p] <- lower[p]
             x[too.high,p] <- upper[p]
           },
           fold = {
             ## ------- New approach that maintains detailed balance ----------
             ## TODO: may still violate bounds if x>2*upper
             x[too.low,p] <- upper[p]-(lower[p]-x[too.low,p])
             too.high<-which(x[,p]>upper[p])
             x[too.high,p] <- lower[p]+(x[too.high,p]-upper[p])
           },
           none = x,
           rand = {
             x[c(too.low,too.high),p] <- lower[p] +  runif(1)*(upper[p]-lower[p])
           },
           stop("Unrecognised value of 'bound.handling'")
           )#switch
    ##if (bound.handling!="none") stopifnot(all(x[,p]>=lower[p] & x[,p]<=upper[p]))
    if (bound.handling!="none" && !all(x[,p]>=lower[p] & x[,p]<=upper[p])) {
      warning("Bounds violated after correction, using random")
      x <- handleBounds(x,lower,upper,"rand")
    }
  } ##for p
  return(x)
}#handleBounds
