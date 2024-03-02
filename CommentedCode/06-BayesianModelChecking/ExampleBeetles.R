############ CREATE ZERO-INFLATED GLMM DATA #################

# This first part creates a dataset with beetles counts across an altitudinal gradient (several plots each observed several years), with a random intercept on year and zero-inflation. 

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

plot(year, altitude, cex = beetles/50, pch =2, main = "Beetle counts across altitudinal gradient\n triangle is proportional to counts")

############## Analysis with JAGS ############################

library(R2jags)
modelData=as.list(data)
modelData = append(data, list(nobs=1000, nplots = 50, nyears = 20))
head(data)

# 1) Fit GLM only 

modelstring="
model {
  
  # Likelihood
  for (i in 1:nobs) {
    lambda[i] <- exp(intercept + alt * altitude[i] + alt2 * altitude[i] * altitude[i]) 
    beetles[i]~dpois(lambda[i]) 
  }
  
  # Fixed effect priors 
  intercept ~ dnorm(0,0.0001)
  alt ~ dnorm(0,0.0001)
  alt2 ~ dnorm(0,0.0001)

  # Posterior predictive simulations 
  
  for (i in 1:nobs) {
    beetlesPred[i]~dpois(lambda[i])
  }
  Prediction <- sum(beetlesPred)
}
"

model=jags(model.file = textConnection(modelstring), data=modelData, n.iter=10000,  parameters.to.save = c("intercept", "alt", "alt2", "beetlesPred", "lambda"), DIC = F)

library(DHARMa)
simulations = model$BUGSoutput$sims.list$beetlesPred
pred = apply(model$BUGSoutput$sims.list$lambda, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = data$beetles, fittedPredictedResponse = pred, integerResponse = T)
plotSimulatedResiduals(sim)


# 2) GLMM with random intercept on year, observation-level RE for overdispersion, and zero-inflation

modelstring="
model {

# Likelihood
for (i in 1:nobs) {
lambda[i] <- exp(intercept + alt * altitude[i] + alt2 * altitude[i] * altitude[i] + Ryear[year[i]] + RID[i] ) * Zero[i] + 0.00000001 

beetles[i]~dpois(lambda[i]) 
}

# Fixed effect priors 
intercept ~ dnorm(0,0.0001)
alt ~ dnorm(0,0.0001)
alt2 ~ dnorm(0,0.0001)

# Random effects 

for (i in 1:nyears) {
Ryear[i]~dnorm(0,sigmaYear)
}

for (i in 1:nobs) {
RID[i]~dnorm(0,sigmaID)
}

# Variance priors 
sigmaYear~dgamma(1,2)
sigmaID~dgamma(0.001,0.001)

# Zeroinflation

for (i in 1:nobs) {
Zero[i]~dbern(zeroMu + altZero * altitude[i])
}
zeroMu ~ dunif(0,1)
altZero ~ dnorm(0,0.0001)

# Posterior predictive simulations 
for (i in 1:nobs) {
  beetlesPred[i]~dpois(lambda[i])
}

  Prediction <- sum(beetlesPred)
}
"

model=jags(model.file = textConnection(modelstring), data=modelData, n.iter=10000,  parameters.to.save = c("intercept", "alt", "alt2", "beetlesPred", "Ryear", "sigmaYear", "lambda", "altZero", "zeroMu"), DIC = F)

library(DHARMa)
simulations = model$BUGSoutput$sims.list$beetlesPred
pred = apply(model$BUGSoutput$sims.list$lambda, 2, median)
dim(simulations)
sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = data$beetles, fittedPredictedResponse = pred, integerResponse = T)
plotSimulatedResiduals(sim)


