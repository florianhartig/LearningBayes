print.dream <- function(x,...){
  cat("
Call:
")
  
  print(x$call)
  
  cat("
Control:
")
  for (i in names(x$control)){
    v <- x$control[[i]]
    cat(sprintf("%15s: ",i))
    cat(v,"\n")
  } 

  cat("\nExit condition:",x$EXITMSG,"\n")

}
