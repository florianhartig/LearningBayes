##' This function calculates the parameters for the exponential power density
##' Equation [20] paper by Thiemann et al. WRR 2001, Vol 37, No 10, 2521-2535
##' @param beta. scalar
##' @return ... Cb, Wb. scalars
## TODO: input valid range?
CalcCbWb <- function(beta){
 
  ## First calculate some dummy variables
  A1 <- gamma(3*(1+beta)/2)
  A2 <- gamma((1+beta)/2)
  ## And use these to derive Cb and Wb
  Cb <- (A1/A2)^(1/(1+beta))
  Wb <- sqrt(A1)/((1+beta)*(A2^1.5))
  return(list(Cb=Cb,Wb=Wb))
}

