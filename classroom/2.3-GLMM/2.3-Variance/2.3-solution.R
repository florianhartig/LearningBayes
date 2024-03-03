library(DHARMa)
library(glmmTMB)
library(rjags)


Data = list(SiblingNegotiation = Owls$SiblingNegotiation, 
            SexParent = as.numeric(Owls$SexParent)-1, # dummy coding!
            FoodTreatment = as.numeric(Owls$FoodTreatment)-1,
            LogBroodSize = log(Owls$BroodSize),
            Nest = as.numeric(Owls$Nest),
            nobs = nrow(Owls),
            nNests = length(levels(Owls$Nest)))



modelCode = "model{

for(i in 1:nobs){
SiblingNegotiation[i] ~ dpois(lambda[i])  # poisson error distribution
lambda[i] = exp(eta[i] + Rnest[Nest[i]] + OLRE[i] ) * Zero[i]
eta[i] <- intercept + EffectSexParent*SexParent[i] + EffektFoodTreatment*FoodTreatment[i] + EffektInterSexFood*SexParent[i]*FoodTreatment[i] + LogBroodSize[i]  
}

intercept ~ dnorm(0,0.0001)
EffectSexParent ~ dnorm(0,0.0001)
EffektFoodTreatment ~ dnorm(0,0.0001)

EffektInterSexFood ~ dnorm(0,0.0001)

# Random effect
for (i in 1:nNests) {
Rnest[i]~dnorm(0,tauNest)
}
tauNest~dgamma(0.001,0.001)

# Observation-level RE for obverdispersion
for (i in 1:nobs) {
OLRE[i]~dnorm(0,tauOD)
}
tauOD~dgamma(0.001,0.001)

# Zeroinflation

for (i in 1:nobs) {
Zero[i]~dbern(zeroMu)
}
zeroMu ~ dunif(0,1)

# Posterior predictive simulations 
for (i in 1:nobs) {
SiblingNegotiationPred[i]~dpois(lambdaSim[i])
lambdaSim[i] = exp(eta[i] + RnestSim[Nest[i]] + OLRESim[i] ) * ZeroSim[i]
}

for (i in 1:nNests) {
RnestSim[i]~dnorm(0,tauNest)
}

for (i in 1:nobs) {
OLRESim[i]~dnorm(0,tauOD)
ZeroSim[i]~dbern(zeroMu)
}

}"

inits.fn <- function() list(Zero = rep(1,nrow(Owls)))
jagsModel <- jags.model(file= textConnection(modelCode), inits = inits.fn, data=Data, n.chains = 3)
para.names <- c("intercept","EffectSexParent", "EffektFoodTreatment", "EffektInterSexFood", "zeroMu", "tauOD", "tauNest")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

gelman.diag(Samples)

summary(Samples)

# with so many variables and random effects, having all observed 
# outputs together becomes somewhat tedious. Unfortunately, the
# package rjags (which I think is the easiest to use) has no 
# alternative. You can, hower, use the package R2jags, which 
# returns a list with observed variables, so if your observed 
# variable is a vector with 100 entries, you get them all together
# moreover, the package also 

para.names <- c("eta", "SiblingNegotiationPred")

Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

library(BayesianTools)
library(DHARMa)
x = getSample(Samples)
# note - yesterday, we calcualted the predictions from the parameters
# here we observe them direct - this is the normal way to calcualte the 
# posterior predictive distribution
posteriorPredDistr = x[,1:599]
posteriorPredSim = x[,600:1198]


sim = createDHARMa(simulatedResponse = t(posteriorPredSim), observedResponse = Owls$SiblingNegotiation, fittedPredictedResponse = apply(posteriorPredDistr, 2, median), integerResponse = T)
plot(sim)



