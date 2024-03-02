
# dataset airquality

plot(Ozone ~ Temp, data = airquality)


# Frequentist analysis

fit <- lm(Ozone ~ Temp, data = airquality)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))
par(mfrow = c(2,2))
plot(fit)


# Bayesian analysis

library(rjags)

dat = airquality[complete.cases(airquality),] 
Data = list(y = as.vector(scale(dat$Ozone)), x = as.vector(scale(dat$Temp)), i.max = nrow(dat))

modelCode = "model{

  # Likelihood
  for(i in 1:i.max){
    mu[i] <- a*x[i]+ b
    y[i] ~ dnorm(mu[i],tau)
  }

  # Prior distributions
  
  # For location parameters, normal choice is wide normal
  a ~ dnorm(0,0.0001)
  b ~ dnorm(0,0.0001)

  # For scale parameters, normal choice is decaying
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)

}
"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)

para.names <- c("a","b","sigma")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)
summary(Samples)

###########################################
# Model selection 1 - no model selection 

# is a significant? Calculate the proportion of a > 0

library(BayesianTools)
mean(getSample(Samples, which = 1) > 0)

# 100% of the posterior > 0, i.e. we are quite certain the effect of Temp on Ozone is positive

###########################################
# Model selection 2 - regularization 

# in this case, I put an adaptive shrinkage prior 
# on the regression parameters. Alternative is to
# fix the level of shrinkage ad hoc

modelCode = "model{

# Likelihood
for(i in 1:i.max){
  mu[i] <- a*x[i]+ b
  y[i] ~ dnorm(mu[i],tau)
}

# Prior distributions

# normal with shrinkage
a ~ dnorm(0,tauShrinkage)
b ~ dnorm(0,tauShrinkage)

tauShrinkage ~ dgamma(0.001, 0.001)
sdShrinkage <- 1/sqrt(tauShrinkage)

# For scale parameters, normal choice is decaying
tau ~ dgamma(0.001, 0.001)
sigma <- 1/sqrt(tau)

}
"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)

para.names <- c("a","b","sigma", "sdShrinkage")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)
summary(Samples)


###########################################
# Model selection 3 - DIC


dic = dic.samples(jagsModel, n.iter = 5000)
# penalty is the effective number of parameters of this model


###########################################
# Model selection 4 - Bayes factor and RJ-MCMC

# Not so easy to do with JAGS - see example in appendix of 
# Dormann et al. 2018 (in the readings) or in the vignette of 
# the BayesianTools package 

