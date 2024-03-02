library(lme4)
library(glmmTMB)
library(DHARMa)

Owls2 = Owls
Owls2$id = 1:nrow(Owls)

m1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest) , data=Owls2 , family = nbinom1, ziformula = ~ 1)
summary(m1)

res <- simulateResiduals(m1)
plot(res)
testDispersion(res)

# How to do the same residual checks in Jags - just for model 1

library(rjags)


Data = list(SiblingNegotiation = Owls$SiblingNegotiation, 
            SexParent = as.numeric(Owls$SexParent)-1, # dummy coding!
            FoodTreatment = as.numeric(Owls$FoodTreatment)-1,
            LogBroodSize = log(Owls$BroodSize),
            Nest = as.numeric(Owls$Nest),
            numNests = nlevels(Owls$Nest),
            nobs = nrow(Owls))


modelCode = "model{

for(i in 1:nobs){
  SiblingNegotiation[i] ~ dpois(lambda[i] * zeroInf[i])  
  
  lambda[i] <- exp(eta[i]) # inverse link function
  
  eta[i] <- intercept + 
            EffectSexParent*SexParent[i] + 
            EffektFoodTreatment*FoodTreatment[i] + 
            EffektInterSexFood*SexParent[i]*FoodTreatment[i] +                LogBroodSize[i] +
            randomEffect[Nest[i]] + 
            OLRE[i]

  zeroInf[i] ~ dbern(pZero)

}

  pZero ~ dunif(0,1)

# Random Effect specification
for(i in 1:numNests){
  randomEffect[i] ~ dnorm(0, tauRE)
}
tauRE ~ dgamma(0.001,0.001)
sigmaRE <- 1/sqrt(tauRE)

for(i in 1:nobs){
  OLRE[i] ~ dnorm(0, tauOLRE)
}
tauOLRE ~ dgamma(0.001,0.001)
sigmaOLRE <- 1/sqrt(tauOLRE)


intercept ~ dnorm(0,0.0001)
EffectSexParent ~ dnorm(0,0.0001)
EffektFoodTreatment ~ dnorm(0,0.0001)
EffektInterSexFood ~ dnorm(0,0.0001)

# Posterior predictive simulations 
for (i in 1:nobs) {
  SiblingNegotiationPred[i]~dpois(lambda[i] * zeroInf[i])
}

}"

inits.fn <- function() list(zeroInf = rep(1,nrow(Owls)))

jagsModel <- jags.model(file= textConnection(modelCode), 
                        data=Data, 
                        inits = inits.fn,
                        n.chains = 3)
para.names <- c("intercept","EffectSexParent", "EffektFoodTreatment", "EffektInterSexFood", "sigmaRE", "sigmaOLRE" , "pZero")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)



library(BayesianTools)
library(DHARMa)
x = getSample(Samples)
# note - yesterday, we calcualted the predictions from the parameters
# here we observe them direct - this is the normal way to calcualte the 
# posterior predictive distribution
posteriorPredDistr = x[,5:(4+599)]
posteriorPredSim = x[,(5+599):(4+2*599)]


sim = createDHARMa(simulatedResponse = t(posteriorPredSim), observedResponse = Owls$SiblingNegotiation, fittedPredictedResponse = apply(posteriorPredDistr, 2, median), integerResponse = T)
plot(sim)



