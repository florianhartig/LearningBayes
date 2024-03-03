library(rjags)

dat = airquality[complete.cases(airquality),] 
dat = dat[order(dat$Temp),] # order so that we can later make more convenient line plots

Data = list(y = dat$Ozone, 
            x = dat$Temp, 
            i.max = nrow(dat))

modelCode = "
model{
  # Likelihood
  
  for(i in 1:i.max){
    mu[i] <- Temp*x[i]+ intercept
    y[i] ~ dnorm(mu[i],tau)
  }
  sigma <- 1/sqrt(tau)

  # Priors
  intercept ~ dnorm(0,0.0001)
  Temp ~ dnorm(0,0.0001)
  tau ~ dgamma(0.001, 0.001)
  
  # Prior predictive - sample new parameters from prior and make predictions
  
  intercept2 ~ dnorm(0,0.0001)
  Temp2 ~ dnorm(0,0.0001)
  
  for(i in 1:i.max){
    muPrior[i] <- Temp2*x[i]+ intercept2
  }
  
  # Posterior predictive - same as likelihood, but remove data
  
  for(i in 1:i.max){
    yPosterior[i] ~ dnorm(mu[i],tau)
  }
}
"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)
update(jagsModel, n.iter = 1000)
Samples <- coda.samples(jagsModel, variable.names = c("intercept","Temp","sigma"), n.iter = 5000)


# Prior predictive analysis

Samples <- coda.samples(jagsModel, variable.names = c("muPrior"), n.iter = 5000)

pred <- getSample(Samples, start = 300)

# plotting the distributions of predictions 
plot(Ozone ~ Temp, data = dat, ylim = c(-200, 300))
for(i in 1:nrow(x)) lines(dat$Temp, pred[i,])

# what we see here: a priori many regression lines are possible
# you can change your prior to be more narrow and see how this
# would push prior space in a certain area 

# Posterior predictive analysis

# Here we can choose to observe the posterior mean predictions 
# or the observations. We will do both in this case because
# both are inputs to the DHARMa plots

Samples <- coda.samples(jagsModel, variable.names = c("mu", "yPosterior"), n.iter = 5000)

library(BayesianTools)
library(DHARMa)

x <- getSample(Samples, start = 300)

dim(x)

# note - yesterday, we calculated the predictions from the parameters
# here we observe them direct - this is the normal way to calcualte the 
# posterior predictive distribution
posteriorPredDistr = x[,1:111]
posteriorPredSim = x[,112:222]

sim = createDHARMa(simulatedResponse = t(posteriorPredSim), observedResponse = dat$Ozone, fittedPredictedResponse = apply(posteriorPredDistr, 2, median), integerResponse = F)
plot(sim)

# all additional plots in the DHARMa package are possible

# essentially, we see here the same issue as in the standard residual plots 

fit <- lm(Ozone ~ Temp, data = dat)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))
par(mfrow = c(2,2))
plot(fit) # residuals





