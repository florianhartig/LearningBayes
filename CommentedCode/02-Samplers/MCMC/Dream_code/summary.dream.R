
summary.dream <- function(object, fraction = 0.5, ...){

  coda.sum <- summary(window(object, fraction = fraction), ...)
  
  cat(sprintf("
Exit message:  %s
Num fun evals: %d
Time (secs):   %.1f
Final R.stats:
",
              object$EXITMSG,
              object$fun.evals,
              object$time
              ))

  R.stat.last <- tail(object$R.stat,1)[,-1]
  if (all(R.stat.last<0)) {
    R.stat.last <- gelman.diag(object$Sequences)$psrf[,1]
  }
  names(R.stat.last) <- colnames(object$R.stat)[-1]
    
  for (i in 1:length(R.stat.last)){
    cat(sprintf("\t%s:\t%f\n",
                names(R.stat.last)[i],
                R.stat.last[i]
                ))
  } ##for    

  cat("
CODA summary for last 50% of MCMC chains:
")

  print(coda.sum)

  cat("
Acceptance Rate
")
  print(summary(object$AR[,2]))

} ##summary.dream
