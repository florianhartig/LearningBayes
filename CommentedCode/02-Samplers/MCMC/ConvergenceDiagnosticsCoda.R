##############################################################
# Code to demonstrate convergence diagnostics with the CODA packaage in R
#
# see blog post at http://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
# for further explanations
#
# by Florian Hartig http://florianhartig.wordpress.com/
# licence: http://creativecommons.org/licenses/by-nc-sa/3.0/
##############################################################

rm(list=ls(all=TRUE))

# This will load a function to have a more sophisticated pair plot
# Can be ommitted if there is not internet connection
library(devtools)
source_url('https://raw.githubusercontent.com/florianhartig/Cookbook/master/Plotting/Correlation/CorrelationDensityPlotWithIpanel.r')

library(coda)

trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)  
}

prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(mcmc(chain))
}

startvalue = c(4,2,8)
chain = run_metropolis_MCMC(startvalue, 10000)

summary(chain)
plot(chain)
betterPairs(data.frame(chain))
# if this doesn't work, rep

################################
# correlation in posterior

x <- (-(sampleSize-1)/2):((sampleSize-1)/2) + 20
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

chain = run_metropolis_MCMC(startvalue, 10000)
summary(chain)
plot(chain)
betterPairs(data.frame(chain))



######################################
# convergence diagnostics

x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

# first chain
startvalue = c(4,2,8)
chain = run_metropolis_MCMC(startvalue, 10000)

# create second chain
startvalue2 = c(6,-2,12)
chain2 = run_metropolis_MCMC(startvalue, 10000)
combinedchains = mcmc.list(chain, chain2)

# check convergennce
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)










