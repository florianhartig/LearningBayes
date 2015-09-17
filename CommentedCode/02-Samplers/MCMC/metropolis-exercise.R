########################################
# Create our data 

trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

x <- (-(sampleSize-1)/2):((sampleSize-1)/2) 
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

plot(x,y)

#######################################
# Likelihood function


# param = c(a,b,sd)
likelihood <- function(param){  
  a <- param[1]
  b <- param[2]
  sd <- param[3]
  
  ymean = a * x + b
  residuals = ymean - y
  singlePointLikelihoods <- dnorm(residuals, mean = 0, sd = sd, log = T)
  return(sum(singlePointLikelihoods))
}

# Example: plot the likelihood profile of the slope a
slopevalues <- function(x){return(likelihood(c(x, trueB, trueSd)))}
slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")

########################################
# Prior 


prior <- function(param){
  a <- param[1] # Slope
  b <- param[2] # Intercept
  sd <- param[3] # Standard deviation
  
  aprior <- dnorm(a, sd = 100, log = T)
  bprior <- dnorm(b, sd = 100, log = T)  
  
  precicion <- 1/sd^2
  
  sdprior <- dunif(precicion, 0,30, log = T)
  
  return(aprior + bprior + sdprior)
}
prior(c(3,7,10))

#########################################
# Posterior

posterior <- function(param){
  if (param[3] <= 0) return(-Inf)
  else return (likelihood(param) + prior(param))
}

posterior(c(3,7,-10))

#########################################
# MCMC



run_metropolis_MCMC <- function(startvalue, iterations){
  
  chain <- matrix(nrow = iterations + 1, ncol = length(startvalue))
  chain[1,] <- startvalue
  
  for (i in 1:iterations){
    proposal <- rnorm(3, mean = chain[i,], sd = c(0.1,0.3,0.3))
    
    jumpProb <- exp(posterior(proposal) - posterior(chain[i,]))
    
    if (runif(1,0,1) <= jumpProb) chain[i+1,] = proposal
    else chain[i+1,] = chain[i,]
  }
  
  return(chain)
}

startvalue <- c(3,7,10)

chain <- run_metropolis_MCMC(startvalue, 5000)

library(coda)

codachain <- mcmc(chain)
plot(codachain)
summary(codachain)

startvalue2 <- c(6,-2,5)
chain2 <- run_metropolis_MCMC(startvalue2, 5000)

codalist <- mcmc.list(mcmc(chain[1000:5000,]), mcmc(chain2[1000:5000,]))

plot(codalist)

gelman.diag(codalist)

summary(lm(y~x))

#####################################
# Correlation plots 
# install.packages("IDPmisc")
library(IDPmisc)

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}

betterPairs(chain)


