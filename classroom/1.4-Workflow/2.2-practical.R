## Residual checkes

# We continue with the owl data set - how would you actually check the residuals of a GLM?
# The best way are simulated quantile residuals, which are implemented in DHARMa
# More on DHARMa residual checks in https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

# The following shows a sequence of models, all checked with DHARMa. The example is discussed in a talk at ISEC 2018, see slides [here](https://www.slideshare.net/florianhartig/mon-c5hartig2493). We will discuss about the exact model structures later. 

library(lme4)
library(glmmTMB)
library(DHARMa)

# Fitting the actual hypothesis
m0 <- glm(SiblingNegotiation ~ FoodTreatment, data=Owls , family = poisson)
res <- simulateResiduals(m0)
plot(res)

# Fitting the actual hypothesis
m1 <- glm(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)), data=Owls , family = poisson)
res <- simulateResiduals(m1)
plot(res)

# Adding random effect
m2 <- glmer(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), data=Owls , family = poisson)
res <- simulateResiduals(m2)
plot(res)

# Switching to nbinom1 to account for overdispersion
m3 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), data=Owls , family = nbinom1)
res <- simulateResiduals(m3)
plot(res)
summary(m3)
plotResiduals(Owls$FoodTreatment, res$scaledResiduals)
testDispersion(res)
testZeroInflation(res)

# Adjusting for zero-inflation
m4 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), ziformula = ~ FoodTreatment + SexParent,  data=Owls , family = nbinom1)
summary(m4)
res <- simulateResiduals(m4)
plot(res)
testDispersion(res)
testZeroInflation(res)
plotResiduals(Owls$FoodTreatment, res$scaledResiduals)

# Residuals look fine. The remaining pattern in res ~ pred is probably due to a 
# problem described in https://github.com/florianhartig/DHARMa/issues/16 and 
# does not show a model specification problem

# Despite residuals looking fine, we can also see if we find support for a dependence
# of the variance on a predictor
m5 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), dispformula = ~ FoodTreatment , ziformula = ~ FoodTreatment + SexParent,  data=Owls , family = nbinom1)
summary(m5)
res <- simulateResiduals(m4)
plot(res)


# How to do the same residual checks in Jags - just for model 1

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

  # Posterior predictive simulations 
  for (i in 1:nobs) {
    SiblingNegotiationPred[i]~dpois(lambda[i])
  }

}"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)
para.names <- c("intercept","EffectSexParent", "EffektFoodTreatment", "EffektInterSexFood","lambda", "SiblingNegotiationPred")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

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



