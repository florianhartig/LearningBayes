set.seed(123)
rm(list=ls(all=TRUE))
library(rjags)
library(R2jags)


## Creation of the data

#Assume we observe data from an ecological system that creates an exponential size distribution (e.g. tree sizes, see [Taubert, F.; Hartig, F.; Dobner, H.-J. & Huth, A. (2013) On the Challenge of Fitting Tree Size Distributions in Ecology. PLoS ONE, 8, e58036-](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0058036)), but our measurments are performed with a substantial lognormal observation error


meanSize <- 10
trueLogSd <- 1
sampleSize <- 500
truevalues = rexp(rate = 1/meanSize, n = sampleSize)
observations = rlnorm(n = length(truevalues), mean = log(truevalues), sd = trueLogSd)


#Plotting true and observed data

maxV <- ceiling(max(observations,truevalues))
counts <- rbind(
  obs = hist(observations, breaks = 0:maxV, plot = F)$counts,
  true = hist(truevalues, breaks = 0:maxV, plot = F)$counts
)
barplot(log(t(counts)+1), beside=T)


## Fitting a non-hierarchical model leads to bias


#Model specification of a non-hierarchical model in JAGS that does not account for the observation error 

normalModel = textConnection('
                             model {
                             # Priors

                             
                             # Likelihood

                             }
                             ')

# Bundle data
positiveObservations <- observations[observations>0]
data = list(true = positiveObservations, nObs=length(positiveObservations))

# Parameters to be monitored (= to estimate)



## Fitting a hierarchical model removes the bias

#Model specification if hierarchical model that accounts for the observation error in Jags

hierarchicalModel = textConnection('
                                   model {
                                   # Priors

                                   
                                   # Likelihood

                                   }
                                   ')
# Bundle data




