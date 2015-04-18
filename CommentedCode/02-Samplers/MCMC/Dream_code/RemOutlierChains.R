##' Finds outlier chains
##' Note: differs from matlab implementation in that it does not remove outliers
##'  removal is done in main dream function
##
##' @param X matrix nseq x ndim+2
##' @param hist.logp matrix [1,ndraw/nseq] x nseq
##' @param control. elements outlierTest,nseq,ndim

##' @return ... list
##'   chain.id. vector/scalar. length [0,nseq]. range [1,nseq]
##'   mean.hist.logp vector. length nseq
##
## depends: ACR.R
##
## TODO: Could pass only x, not X. Could return r.idx, not mean.hist.logp
RemOutlierChains <- function(X,hist.logp,control){
  
  ## Dimensions:
  ## idx.end,idx.start
  ##
  ## q13. vector length 2
  ## iqr,upper.range. scalar
  ## G,t2,Gcrit. scalar
  ## alpha,upper.range,idx,d1 scalar

  
  ## Determine the number of elements of L_density
  idx.end <- nrow(hist.logp)
  idx.start <- floor(0.5*idx.end)
  
  ## Then determine the mean log density of the active chains
  mean.hist.logp <- colMeans(hist.logp[idx.start:idx.end,])
  
  ## Initialize chain_id and Nid
  chain.id <- NULL
  
  ## Check whether any of these active chains are outlier chains
  switch(control$outlierTest,
         'IQR_test'={         
           ## Derive the upper and lower quantile of the data
           q13<-quantile(mean.hist.logp,c(0.75,0.25))
           ## Derive the Inter quartile range
           iqr <- q13[1]-q13[2]
           ## Compute the upper range -- to detect outliers
           ## TODO: shouldn't upper.range be Q3+3*IQR?
           upper.range <- q13[2]-2*iqr
           ## See whether there are any outlier chains
           chain.id <- which(mean.hist.logp<upper.range)
         },
         'Grubbs_test'={
           ## Test whether minimum log_density is outlier
           G <- (mean(mean.hist.logp)-min(mean.hist.logp))/sd(mean.hist.logp)
           
           ## Determine t-value of one-sided interval
           t2 = qt(1 - 0.01/control$nseq,control$nseq-2)^2; ## 95% interval
           
           ## Determine the critical value
           Gcrit <- ((control$nseq-1)/sqrt(control$nseq))*sqrt(t2/(control$nseq-2+t2))
           
           ## Then check this
           if (G > Gcrit) { ## Reject null-hypothesis
             chain.id <- which.min(mean.hist.logp)
           }
         },
         'Mahal_test'={
           ## Use the Mahalanobis distance to find outlier chains
           alpha <- 0.01
           upper.range <- ACR(control$ndim,control$nseq-1,alpha)
           ## Find which chain has minimum log_density
           idx <- which.min(mean.hist.logp)
           ## Then check the Mahalanobis distance
           ## TODO: MATLAB: d1 = mahal(X(idx,1:MCMCPar.n),X(ii,1:MCMCPar.n));
           d1 <- mahalanobis(X[idx,control$ndim],
                             center=colMeans(X[-idx,control$ndim]),
                             cov=cov(X[-idx,control$ndim])
                             )
           ## Then see whether idx is an outlier in X
           if (d1>upper.range) {
             chain.id <- idx
           }
         },
         'None'={
           stop("Outlier detection reached when it should have been turned off")
         },
         stop("Unknown outlierTest specified")
         ) ##switch


  return(list(chain.id=chain.id,mean.hist.logp=mean.hist.logp))

} ##RemOutlierChains
