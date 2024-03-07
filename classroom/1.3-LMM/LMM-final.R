rm(list = ls())
Dat = read.table("https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/Aspis_data.txt", stringsAsFactors = T)

# Inspect relationship between body mass and total body lenght
plot(Dat$TL, Dat$BM,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')

# For the analysis we use log-transformed body masses
# and log-transformed and scaled total body lenght (TL)
plot(Dat$log_TL.sc, Dat$log_BM)

# Body mass may differ between females and males
plot(Dat$TL, Dat$BM,
     pch = as.numeric(Dat$Sex)+1,
     col = as.numeric(Dat$Sex)+1,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')

# Also, body sizes may be different in different populations
plot(Dat$TL, Dat$BM,
     pch = as.numeric(Dat$Pop)+1,
     col = as.numeric(Dat$Pop)+1,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')


# Frequentist models 
LM <- lm(log_BM ~ log_TL.sc, data = Dat)
summary(LM)
abline(LM)

# Effect and interaction with sex
LM <- lm(log_BM ~ log_TL.sc * Sex, data = Dat)
summary(LM)
library(effects)
plot(allEffects(LM, partial.residuals = T))

# mixed model including also population differences
library(lme4)
library(lmerTest)
LM <- lmer(log_BM ~ log_TL.sc * Sex + (1|Pop), data = Dat)
# LM <- lm(log_BM ~ log_TL.sc * Sex + Pop, data = Dat)
summary(LM)
library(effects)
plot(allEffects(LM, partial.residuals = T))
x = ranef(LM)
hist(x$Pop$`(Intercept)`)
shapiro.test(x$Pop$`(Intercept)`)


# Bayesian Version of the model
library(rjags)
model ="
model{
  # Likelihood
  for(i in 1:n.dat){
    y[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha[Pop[i]] + beta.TL * TL[i] + beta.TLm * TL[i] * Sexm[i] + beta.m * Sexm[i]
  }
  
  for(p in 1:n.pop){
    alpha[p] ~ dnorm(mu.alpha, tau.pop)
  }

  # Prior distributions
  mu.alpha ~ dnorm(0,0.001)
  beta.TL ~ dnorm(0,0.001)
  beta.m ~ dnorm(-0.8,10) 
  beta.TLm ~ dnorm(0,0.001)  
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)
  tau.pop <- 1/(sigma.pop*sigma.pop)
  sigma.pop ~ dunif(0,100)  
}
"

Data = list(y = Dat$log_BM, 
            TL = Dat$log_TL.sc,
            Sexm = ifelse(Dat$Sex == "m", 1, 0),
            n.dat = nrow(Dat),
            Pop = Dat$Pop,
            n.pop = max(Dat$Pop))

jagsModel <- jags.model(file = textConnection(model), data=Data, 
                        n.chains = 3, 
                        n.adapt= 5000)

Samples <- coda.samples(jagsModel, 
                        variable.names = c("mu.alpha", 
                                           "alpha",
                                           "beta.TL", 
                                           "beta.TLm", 
                                           "beta.m", 
                                           "sigma", 
                                           "sigma.pop"), 
                        n.iter = 5000)

summary(Samples)
plot(Samples)
gelman.diag(Samples)
BayesianTools::marginalPlot(Samples)

BayesianTools::correlationPlot(Samples, whichParameters = c(8:12))

