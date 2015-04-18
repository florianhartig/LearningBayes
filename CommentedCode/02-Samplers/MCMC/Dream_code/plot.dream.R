##' Plot characteristics of dream object
##' @param dream object
##' @param interactive - stop for each plot

##' Uses second half of sequences

plot.dream <- function(x, interactive = TRUE, ...){
  opar <- devAskNewPage(interactive)
  on.exit(devAskNewPage(opar))
  
  ss <- window(x, ...)

  ## Convergence (Gelman plot)
  
  tmp <- try(gelman.plot(ss))

  if (inherits(tmp, "try-error")) {
      plot(x$R.stat[,1],x$R.stat[,2],type="l",ylim=c(0,2))
      for (i in 2:x$control$ndim) lines(x$R.stat[,1],x$R.stat[,i+1],ylim=c(0,2))
      title(main="Evolution of R.stat",sub="Equivalent to gelman.plot")
  }

  ## Trace and parameter density and auto-correlation
  
  print(densityplot(ss))

  print(xyplot(ss))

  print(acfplot(ss))

  ## Acceptance rate
  print(barchart(table(x$AR[,2]), main="Distribution of % acceptance rate"))
  
}##plot.dream
