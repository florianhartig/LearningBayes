
library(glmmTMB)
library(rjags)

Data = list(SiblingNegotiation = Owls$SiblingNegotiation, 
            SexParent = as.numeric(Owls$SexParent)-1, # dummy coding!
            FoodTreatment = as.numeric(Owls$FoodTreatment)-1,
            LogBroodSize = log(Owls$BroodSize),
            nobs = nrow(Owls))


modelCode = "model{

for(i in 1:nobs){
SiblingNegotiation[i] ~ dpois(lambda[i])  # poisson error distribution
lambda[i] <- exp(eta[i]) # inverse link function
eta[i] <- intercept + EffectSexParent*SexParent[i] + EffektFoodTreatment*FoodTreatment[i] + EffektInterSexFood*SexParent[i]*FoodTreatment[i] + LogBroodSize[i]       # linear predictor
}

intercept ~ dnorm(0,0.0001)
EffectSexParent ~ dnorm(0,0.0001)
EffektFoodTreatment ~ dnorm(0,0.0001)
EffektInterSexFood ~ dnorm(0,0.0001)

}"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)
para.names <- c("intercept","EffectSexParent", "EffektFoodTreatment", "EffektInterSexFood")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

plot(Samples)
summary(Samples)


