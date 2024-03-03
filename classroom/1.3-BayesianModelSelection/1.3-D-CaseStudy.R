# By Florian Hartig. An extended commented version of this code as well as possible updates are available at http://florianhartig.github.io/LearningBayes/. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License


# Run a model selection on the airquality data 

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
# Add a quadratic effect to the model 


###########################################
# Model selection 1 - no model selection 

# look at posteriors
# can get an idea about the strength of the evidence for an effect
# by using the following code

library(BayesianTools)
mean(getSample(Samples, which = 1) > 0)

# index in which = 1 needs to match the parameter you want 

###########################################
# Model selection 2 - regularization 



###########################################
# Model selection 3 - DIC

dic = dic.samples(jagsModel, n.iter = 5000)
# penalty is the effective number of parameters of this model


###########################################
# Model selection 4 - Bayes factor and RJ-MCMC

# Could use BayesianTools package if you wanted 


