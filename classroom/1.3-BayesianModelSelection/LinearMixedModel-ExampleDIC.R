## ---- echo=F, warning=F, message=F---------------------------------------
set.seed(123)
rm(list=ls(all=TRUE))
library(rjags)
library(runjags)
library(lme4)
library(effects)
library(R2jags)

## ---- fig.width=5, fig.height=5------------------------------------------
a <- 5
b <- 10
sigma <- 10
rsigma = 30
group = rep(1:11, each = 5)
randomEffect = rnorm(11, sd = rsigma)

x <- -27:27
y <- a * x + b + rnorm(55,0,sd = sigma) + randomEffect[group]
plot(x,y, col = group, pch = 3)

## ---- fig.width=5, fig.height=5------------------------------------------
fit <- lm(y ~ x)
summary(fit)
plot(allEffects(fit, partial.residuals = T))

## ---- fig.width=5, fig.height=5------------------------------------------
fit <- lmer(y ~ x + (1|group))
summary(fit)
plot(x,y, col = group,  pch = 3)
for(i in 1:11){
  abline(coef(fit)$group[i,1], coef(fit)$group[i,2], col = i)
}


## ------------------------------------------------------------------------
  # 1) Model definition exactly how we created our data 
  modelCode = "
    model{
      
      # Likelihood
      for(i in 1:i.max){
        y[i] ~ dnorm(mu[i],tau)
        mu[i] <- a*x[i] + b
      }

      # Prior distributions
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)
      tau <- 1/(sigma*sigma)
      sigma ~ dunif(0,100)
    }
  "
  
  # 2) Set up a list that contains all the necessary data (here, including parameters of the prior distribution)
  Data = list(y = y, x = x, i.max = length(y))

  # 3) Specify a function to generate inital values for the parameters
  inits.fn <- function() list(a = rnorm(1), b = rnorm(1), sigma = runif(1,1,100))


## ---- fig.width=7, fig.height=7------------------------------------------
  # Compile the model and run the MCMC for an adaptation (burn-in) phase
  jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3, n.adapt= 1000)

  # Specify parameters for which posterior samples are saved
  para.names <- c("a","b","sigma")

  # Continue the MCMC runs with sampling
  Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
  
  # Plot the mcmc chain and the posterior sample for p
  plot(Samples)

  dic = dic.samples(jagsModel, n.iter = 5000)
  dic
  
## ------------------------------------------------------------------------
gelman.diag(Samples)

## ------------------------------------------------------------------------
summary(Samples)

## ---- fig.width=5, fig.height=5------------------------------------------
plot(x,y)
sampleMatrix <- as.matrix(Samples)
selection <- sample(dim(sampleMatrix)[1], 1000)
for (i in selection) abline(sampleMatrix[i,1], sampleMatrix[i,1], col = "#11111105")


## ------------------------------------------------------------------------
  # 1) Model definition exactly how we created our data 
  modelCode = "
    model{
      
      # Likelihood
      for(i in 1:i.max){
        y[i] ~ dnorm(mu[i],tau)
        mu[i] <- a*x[i] + b + r[group[i]]
      }

      # random effect
      for(i in 1:nGroups){
        r[i] ~ dnorm(0,rTau)
      }

      # Prior distributions
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)

      tau <- 1/(sigma*sigma)
      sigma ~ dunif(0,100)

      rTau <- 1/(rSigma*rSigma)
      rSigma ~ dunif(0,100)
    }
  "
  
  # 2) Set up a list that contains all the necessary data (here, including parameters of the prior distribution)
  Data = list(y = y, x = x, i.max = length(y), group = group, nGroups = 11)

  # 3) Specify a function to generate inital values for the parameters
  inits.fn <- function() list(a = rnorm(1), b = rnorm(1), sigma = runif(1,1,100), rSigma = runif(1,1,100))


## ---- fig.width=7, fig.height=7------------------------------------------
  # Compile the model and run the MCMC for an adaptation (burn-in) phase
  jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3, n.adapt= 1000)

  # Specify parameters for which posterior samples are saved
  para.names <- c("a","b","sigma", "rSigma")

  # Continue the MCMC runs with sampling
  Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
  
  # Plot the mcmc chain and the posterior sample for p
  plot(Samples)

## ----  fig.width=18, fig.height=18---------------------------------------
R2JagsResults <- jags(data=Data, inits=inits.fn, parameters.to.save=c("a","b","sigma", "rSigma", "r"), n.chains=3, n.iter=5000, model.file=textConnection(modelCode))

plot(R2JagsResults)
print(R2JagsResults)


dic = dic.samples(jagsModel, n.iter = 5000)
dic


## ---- fig.width=14, fig.height=14, eval= F-------------------------------
## R2JagsCoda <- as.mcmc(R2JagsResults)
## plot(R2JagsCoda)
## summary(R2JagsCoda)

