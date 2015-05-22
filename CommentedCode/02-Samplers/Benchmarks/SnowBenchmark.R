library(snow)


testFunction<-function(x) {
  for (i in 1:20){
    temp <- rnorm(10000)
    a<- mean(temp)
    b<-median(temp)
    result = a * temp + b
    fit <- lm(result~temp)
  }
  return(1)
}


testparallel <- function(cores=1){
  cl <- makeCluster(cores,type="SOCK")
  elapsed <- system.time(clusterApply(cl, 1:12, testFunction))
  stopCluster(cl)
  return(elapsed[3])
}

cores <- c(1,2,3,4,6,8,12)
cores <- c(1,2,3,4,6)

time <- sapply(cores, testparallel)

par(mfrow = c(2,1))

plot(cores, time, type = "b", main = "total runtime")
plot(cores, time*cores, type = "b", main = "resource time")

result <- data.frame(cores, time, time*cores)
result







