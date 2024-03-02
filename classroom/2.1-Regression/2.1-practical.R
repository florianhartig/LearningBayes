## Owl example (count data)
# The next examples uses the fairly well known Owl dataset which is provided in glmmTMB (see ?Owls for more info about the data).

library(glmmTMB)
library(effects)

m1 <- glm(SiblingNegotiation ~ SexParent, data=Owls , family = poisson)
summary(m1)
plot(allEffects(m1))


# Exercise 1: fit this GLM in Jags! Hint: start with the code from yesterday, and then change the lm to a glm. 
# Hint: SexParent
# sex = as.numeric(Owls$SexParent) - 1 
# sexE * sex[i]



# Exercise 2: implement the following model in Jags
m2 <- glm(SiblingNegotiation ~ FoodTreatment*SexParent + offset(log(BroodSize)), data=Owls , family = poisson)
summary(m2)
plot(allEffects(m2))

















































