
# Assume we flipped a coin 10 times, and want to find out if
# it is biased - what can we learn about the probability to
# obtain heads with this coin?

trials = 10
success = 7

# For all three statistical methods, we use the same statistical model
# which is the binomial model. The probability density is available in R
# through he function dbinom - dbinom(6,10,0.9) gives you the probability
# density of obtaining 6/10 heads when the true probability of heads is 0.9

# Try out the function, and observe how the likelihood changes if you change
# the number of successes and the success probability!

# We will now use this model to calculate the three classical inferential
# outputs of statistics - tests, MLE and Bayes

########### NHST (p-values) #####################

# assume the coin is random (0.5), p-value is p >= observed
barplot(dbinom(0:10,10, 0.5), names.arg = 0:10, col = c(rep("grey", 7), rep("red", 4)))

binom.test(7,10,alternative = "greater")
binom.test(7,10,alternative = "less")
binom.test(7,10)


############ MLE ######################

likelihood = function(x) dbinom(7,10, x)
parameterValues = seq(0,1,length.out = 100)

# assume data is fixed, true probability unknown
plot(parameterValues, likelihood(parameterValues), type = "l")

abline(v = parameterValues[which.max(likelihood(parameterValues))],
       col = "red")


############# Bayes ##################

# posterior

likelihood = function(x) dbinom(7,10, x)
prior = function(x) dnorm(x, mean = 0.2, sd = 0.01)
# prior = function(x) dunif(x, 0,1)


par(mfrow = c(2,2))

plot(parameterValues, prior(parameterValues), type = "l", main = "Prior")

plot(parameterValues, likelihood(parameterValues), type = "l", main = "Likelihood")

plot(parameterValues, prior(parameterValues) * likelihood(parameterValues), type = "l", main = "Posterior")

# AFTERNOON

plot(Ozone ~ Temp, data = airquality)

fit <- lm(Ozone ~ scale(Temp), data = airquality)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))

summary(airquality)

airqualityCleaned = airquality[complete.cases(airquality),]

# MLE FIT

fitMLE <- lm(Ozone ~ Temp, data = airqualityCleaned)
summary(fitMLE)

# Bayesian fit using brms

library(brms)
fitBRMS <- brm(Ozone ~ Temp, 
               prior =  set_prior("normal(2,0.001)", 
                                  class = "b", 
                                  coef = "Temp"),
               data = airqualityCleaned,
               )

plot(fitBRMS)
?summary.brmsfit
summary(fitBRMS, priors = T)

plot(conditional_effects(fitBRMS), ask = FALSE)

# Priors: 
# * default: uninformative, brings similar to MLE for simple models
# * mildly regularizing: slight preference for 0, e.g. dnorm(mean = 0, sd = 5)
# * informative priors - say 2 studies Study 1: mean 2, se 0.2, Study 2: mean 2.2, se 0.3 -> consensus via meta-analysis 

# If setting mildly informative priors, scale response and predictors
plot(scale(Ozone) ~ scale(Temp), data = airquality)
fit <- lm(scale(Ozone) ~ scale(Temp), data = airquality)
summary(fit)


## STAN ##

# This is the code that is actually fitted by STAN 
fitBRMS$model


library(rstan)

stanmodelcode <- "
  data {
    int<lower=0> N;
    vector[N] Temp;
    vector[N] Ozone;
  }
  parameters {
    real intercept;
    real TempEffect;
    real<lower=0> sigma;
  }
  model {
    Ozone ~ normal(intercept + TempEffect * Temp, sigma);
  }
"
dat = list(Ozone = airqualityCleaned$Ozone, 
           Temp = airqualityCleaned$Temp, 
           N = nrow(airqualityCleaned))

fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 2012, chains = 3, verbose = TRUE,
            sample_file = file.path(tempdir(), 'norm.csv')) 

rstan::traceplot(fit)
plot(fit)
summary(fit)


library(rjags)

modelCode = "
model{
  # Likelihood 
  for(i in 1:nobs){
    meanExpectation[i] <- temperatureEffect*Temp[i] + intercept
    Ozone[i] ~ dnorm(meanExpectation[i], tau)
  }
  
  # Priors
  temperatureEffect~dnorm(0,0.0001)
  intercept~dnorm(0,0.0001)  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
}
"

Data = list(Ozone = airqualityCleaned$Ozone, 
            Temp = airqualityCleaned$Temp, 
            nobs = nrow(airqualityCleaned))


jagsModel <- jags.model(file= textConnection(modelCode),
                        data = Data,
                        n.chains = 3)

Samples <- coda.samples(jagsModel,
                        variable.names = c("temperatureEffect", "intercept", "sigma"), 
                        n.iter = 5000)

plot(Samples)
summary(Samples)




library(BayesianTools)

likelihood <- function(par){
  a0 = par[1]
  a1 = par[2]
  sigma <- par[3]  
  expectation = a0 + a1 * scale(airqualityCleaned$Temp)
  
  logLikel = sum(dnorm(expectation - airqualityCleaned$Ozone , 
                        sd = sigma, 
                        log = T))
  return(logLikel)
}

likelihood(c(-140, 2.3, 20))

setup <- createBayesianSetup(likelihood = likelihood, 
                             lower = c(-200,-200,0.01), 
                             upper = c(200,200,30), 
                             names = c("a0", "a1", "sigma"))

out <- runMCMC(setup)
summary(out)
plot(out, start = 1000)
marginalPlot(out)
correlationPlot(out, start = 1000)

x = getSample(out, start = 1000, coda = T)
library(coda)
coda::gelman.diag(x)
coda::autocorr.plot(x)



