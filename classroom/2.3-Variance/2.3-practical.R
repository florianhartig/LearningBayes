library(DHARMa)
library(glmmTMB)
library(lme4)
library(rjags)

# Implement the following frequentist models with Jags, starting with our previous Jags model

# GLMM with random effect on test
m2 <- glmer(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), data=Owls , family = poisson)
res <- simulateResiduals(m2)
plot(res)

# This, but instead of the neg binomial, use an observation level RE
m3 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), data=Owls , family = nbinom1)
res <- simulateResiduals(m3)
plot(res)
summary(m3)
plotResiduals(Owls$FoodTreatment, res$scaledResiduals)
testDispersion(res)
testZeroInflation(res)

# Add zero-inflation
m4 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), ziformula = ~ 1,  data=Owls , family = nbinom1)
summary(m4)
res <- simulateResiduals(m4)
plot(res)
testDispersion(res)
testZeroInflation(res)
plotResiduals(Owls$FoodTreatment, res$scaledResiduals)


# Bonus: make ZI and variance dependent on predictors
m5 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)) + (1|Nest), dispformula = ~ FoodTreatment , ziformula = ~ FoodTreatment + SexParent,  data=Owls , family = nbinom1)
summary(m5)
res <- simulateResiduals(m4)
plot(res)


# Then, try to extend 