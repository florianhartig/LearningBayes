##############################################################
# Code to demonstrate adaptive MCMC via covariance adjustment of the proposal
#
# Reference Rosenthal, J. S. (2011) Optimal proposal distributions and adaptive MCMC. Handbook of Markov Chain Monte Carlo, CRC Press.
#
# by Florian Hartig http://florianhartig.wordpress.com/
# licence: http://creativecommons.org/licenses/by-nc-sa/3.0/
##############################################################

rm(list=ls(all=TRUE))

library(devtools)
source_url('https://raw.githubusercontent.com/florianhartig/Cookbook/master/Plotting/Correlation/CorrelationDensityPlotWithIpanel.r')
library(coda)
library(MASS)

trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# creating uncentered predictor variable, will lead to correlation
# between intercept and slope in the posterior
x <- (-(sampleSize-1)/2):((sampleSize-1)/2) + 20
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)


# the following definitions are standard for a linear regression
# if you don't understand why, read 

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

# Now comes the mcmc to sample from the posterior
# We are using a Metropolis-Hastings MCMC. I you are uncertain how this works, first read 
# http://theoreticalecology.wordpress.com/
# In this previous version, the proposal function was simlply creating independent normal draws
# 
# proposalfunction <- function(param){
#   return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
# }
#
# we change now to a multivariate normal version that will be adapted later

sig = diag(x = c(0.1,0.5,0.3), nrow=3, ncol=3)

proposalfunction <- function(param){
  return(mvrnorm(1,mu = param, Sigma= sig))
}

# Here is the mcmc. I adapted the code slightly to allow for running the MCMC
# for a while, stopping, and continuing with the MCMC, which will be useful 
# for the adaptation

run_metropolis_MCMC <- function(iterations){
  startindex = nrow(chain)
  chain = rbind(chain, array(dim = c(iterations,3)))
  for (i in startindex:(startindex+iterations-1)){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

# but first of all, we will just run the MCMC like before and check the convergence as explained 
# in http://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/

chain = array( c(4,2,8), dim = c(1,3))

# Running a non-adapted analysis with deliberately bad proposal function 

# settings for the proposal covariance matrix
sig = diag(x = c(1,1,1), nrow=3, ncol=3)

chain1 = run_metropolis_MCMC(10000)
chain2 = run_metropolis_MCMC(10000)
combinedchains = mcmc.list(mcmc(chain1), mcmc(chain2))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

# adaptation steps 
# we use the samples obtained already to adjust the proposal
# for exaplanations why the scaling factor of 2.38^2/d is optional see 
# Rosenthal, J. S. (2011) Optimal proposal distributions and adaptive MCMC. Handbook of Markov Chain Monte Carlo, CRC Press
chain = array( c(4,2,8), dim = c(1,3))
chain = run_metropolis_MCMC(2000)
sig = 2.38^2 / 3 * cov(chain) 
chain = run_metropolis_MCMC(2000)
sig = 2.38^2 / 3 * cov(chain) 
chain = run_metropolis_MCMC(2000)
sig = 2.38^2 / 3 * cov(chain) 
chain = run_metropolis_MCMC(2000)
sig = 2.38^2 / 3 * cov(chain) 
chain = run_metropolis_MCMC(2000)
sig = 2.38^2 / 3 * cov(chain) 


chain = array( c(4,2,8), dim = c(1,3))
chain1 = run_metropolis_MCMC(10000)
chain2 = run_metropolis_MCMC(10000)
combinedchains = mcmc.list(mcmc(chain1), mcmc(chain2))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)


