rm(list = ls())
Dat = read.table("https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/Aspis_data.txt", stringsAsFactors = T)

# Inspect relationship between body mass and total body lenght
plot(Dat$TL, Dat$BM,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')


# For the analysis we use log-transformed body masses
# and log-transformed and scaled total body lenght (TL)
# because this seems to have less heteroskedasticity
plot(Dat$log_TL.sc, Dat$log_BM)

# Base model
LM <- lm(log_BM ~ log_TL.sc, data = Dat)
summary(LM)
abline(LM)



# Task 1 - base LM --------------------------------------------------------

# Task: Complete the following jags code, run the model, and check that 
# results are identical to the frequentist lm above 

model ="
model{
  # Likelihood
  for(i in 1:n.dat){
    y[i] ~ dnorm(mu[i],tau)

    }
  
  # Prior distributions

  }
"

Data = list()

jagsModel <- jags.model(file = textConnection(model), 
                        data=Data, 
                        n.chains = 3, 
                        n.adapt= 5000)


Samples <- coda.samples(jagsModel, 
                        variable.names = c(), 
                        n.iter = 5000)

summary(Samples)
plot(Samples)
gelman.diag(Samples)


# Task 2 - add categorical predictor  -----------------------------------------


# Body mass may differ between females and males
plot(Dat$TL, Dat$BM,
     pch = as.numeric(Dat$Sex)+1,
     col = as.numeric(Dat$Sex)+1,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')


# Taks: Modify your Jags model to correspond to the following model

LM <- lm(log_BM ~ log_TL.sc * Sex, data = Dat)
summary(LM)
library(effects)
plot(allEffects(LM, partial.residuals = T))

# Hint: reflect about the options to code contrasts for categorical variables,
# see https://theoreticalecology.github.io/AdvancedRegressionModels/2A-LinearRegression.html#categorical-predictors


# Task 3 - add interactions  -----------------------------------------

# Taks: Modify your Jags model to correspond to the following model
# which fits the main effect of AND an interaction with sex

LM <- lm(log_BM ~ log_TL.sc * Sex, data = Dat)
summary(LM)
library(effects)
plot(allEffects(LM, partial.residuals = T))


# Task 3 - mixed effects  -----------------------------------------

#Body sizes may be different in different populations
plot(Dat$TL, Dat$BM,
     pch = as.numeric(Dat$Pop)+1,
     col = as.numeric(Dat$Pop)+1,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')


# Taks: Modify your Jags model to correspond to the following model
# which the previous model with a random intercept on population

library(lme4)
library(lmerTest)
LM <- lmer(log_BM ~ log_TL.sc * Sex + (1|Pop), data = Dat)
# LM <- lm(log_BM ~ log_TL.sc * Sex + Pop, data = Dat) # equivalent FE model
summary(LM)
library(effects)
plot(allEffects(LM, partial.residuals = T))

# note LM assumes that random intercepts come from a normal distribution
x = ranef(LM)
hist(x$Pop$`(Intercept)`)
shapiro.test(x$Pop$`(Intercept)`)

