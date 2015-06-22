# Exercise for scrip "Beetles"

set.seed(2)
library(R2jags)
library(runjags)

## Simulating dataset

altitude = rep(seq(0,1,len = 50), each = 20)
dataID = 1:1000
spatialCoordinate = rep(seq(0,30, len = 50), each = 20)

# random effects + zeroinflation
plot = rep(1:50, each = 20)
year = rep(1:20, times = 50)

yearRandom = rnorm(20, 0, 1)
plotRandom = rnorm(50, 0, 1)
overdispersion = rnorm(1000, sd = 0.5)
zeroinflation = rbinom(1000,1,0.6)

beetles <- rpois(1000, exp( 0  + 12*altitude - 12*altitude^2 
                            #  + overdispersion   + plotRandom[plot]
                            + yearRandom[year]) * zeroinflation )

data = data.frame(dataID, beetles, altitude, plot, year, spatialCoordinate)

# Measured beetle counts over 20 years on 50 different plots across an altitudinal gradient

plot(year, altitude, cex = beetles/50, pch =2, main = "Beetle counts across altitudinal gradient, triangle is proportional to counts")


### Preparation

library(R2jags)
modelData=as.list(data)
modelData = append(data, list(nobs=1000, nplots = 50, nyears = 20))

## Basic model


modelstring="
model {

  # Likelihood
  for (i in 1:nobs) {
    beetles[i]~dpois(lambda[i]) 
    lambda[i] <- (exp(intercept + alt * altitude[i] + alt2 * altitude[i] * altitude[i] + Ryear[year[i]]) ) * Zero[i] + 0.00000001
  }

  # Effect priors 
  intercept ~ dnorm(0,0.0001)
  alt ~ dnorm(0,0.0001)
  alt2 ~ dnorm(0,0.0001)


  # Random effects

  for (i in 1:nyears) {
    Ryear[i]~dnorm(0,sigmaYear)
  }
  sigmaYear~dgamma(1,2)


  # Zeroinflation

  for (i in 1:nobs) {
    Zero[i]~dbern(zeroMu)
  }
  zeroMu ~ dunif(0,1)


   # Predictions
  for (i in 1:nobs) {
    beetlesPred[i]~dpois(lambda[i])
  }
  Prediction <- sum(beetlesPred)


}
"

#Running this

model=jags(model.file = textConnection(modelstring), data=modelData, n.iter=10000,  parameters.to.save = c("intercept", "alt", "alt2", "Prediction", "Ryear", "sigmaYear"), DIC = F)

plot(model, display.parallel = T)




# Plot the results

altitude <- seq(0,1,len = 50)
plot(data$altitude + runif(1000,-0.02,0.02), log(data$beetles + 1 ))


combinedChainValues <- as.data.frame(combine.mcmc(as.mcmc(model)))

for(i in seq(5,nrow(combinedChainValues), 5)){
  response <- exp(combinedChainValues$intercept[i] + combinedChainValues$alt[i] * altitude + combinedChainValues$alt2[i] * altitude^2)
  lines(altitude, log(response + 1), col = "#22222202", lwd = 3)
}

lines(altitude, log(exp(12*altitude - 12*altitude^2) + 1), col = "red" )



# Bayesian p-value
hist(combinedChainValues$Prediction, breaks = 100, xlim = c(0, 30000))
abline(v=sum(data$beetles), col = "red")


