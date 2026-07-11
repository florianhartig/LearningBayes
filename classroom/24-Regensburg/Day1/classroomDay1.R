
dgamma 
 
values = seq(0,10,length.out = 100)
density = dgamma(values, shape = 1) 
plot(values, density)


# 10 parameter model 
# 100^10 = 100000000000000000000 = 10^20

library(BayesianTools)

density = function(x) dgamma(x, shape = 1, log = T) 
setup = createBayesianSetup(density, lower = 0, upper = 10)
out = runMCMC(setup, settings = list(iterations = 10000), sampler = "DEzs")
plot(out)


str(airquality)
plot(Ozone ~ Temp, data = airquality)

fit <- lm(Ozone ~ Temp, data = airquality)
summary(fit)
abline(fit)

par(mfrow = c(2,2))
plot(fit)

# Bayesian Analysis 

airqualityCleaned = airquality[complete.cases(airquality),]

library(brms)

bayes1 = brm(Ozone ~ Temp , data = airqualityCleaned)
summary(bayes1)
plot(bayes1)

bayes1$model

bayes1$

library(rjags)
Data = list(y = airqualityCleaned$Ozone, 
            x = airqualityCleaned$Temp,
            nobs = nrow(airqualityCleaned))

modelCode = "
model{

  # Likelihood 
  for(i in 1:nobs){
    expected[i] <- Temp*x[i] + intercept
    y[i] ~ dnorm(expected[i], tau) 
  }
  
  # priors
  Temp ~ dnorm(0, 0.0001)
  intercept  ~ dnorm(0, 0.0001)
  tau ~ dgamma(0.001, 0.001)
  
  sigma <- 1/sqrt(tau)

}
"

jagsModel <- jags.model(file = textConnection(modelCode), data = Data, n.chains = 3)

update(jagsModel, n.iter = 2000)

Samples <- coda.samples(jagsModel, 
                        variable.names = c("Temp", "intercept", "sigma"),
                        n.iter = 5000)

plot(Samples)
summary(Samples)



library(BayesianTools)

likelihood <- function(par){
  Temp = par[1]
  intercept = par[2]
  sigma = par[3]
  logLikelihood = sum(dnorm(intercept + Temp * airqualityCleaned$Temp 
                            - airqualityCleaned$Ozone,
                            sd = sigma,
                            log = T))
  return(logLikelihood)
}

likelihood(c(2.4,-150, 20))

setup <- createBayesianSetup(likelihood = likelihood, 
                             lower = c(-200, -200,0), 
                             upper = c(200, 200, 50))

out = runMCMC(setup)
plot(out, start = 1500)
summary(out, start = 1500)

x = getSample(out, start = 1500)
mean(x[,1]>0)

gelmanDiagnostics(out, start = 1500) # Calculate after burn-in!

# Paper: we ran 10.000 MCMC iterations of a DEzs with 3 chains 
# We removed 4500 iterations for burn-in
# Convergence was checked visually and with Gelman-Rubin, all psrf was < 1.05, indicating convergence 


# brms example: 4 chains each with iter = 2000; burn-in = 1000; convergence was checked visually and all Rhat values were < 1.05

correlationPlot(out, start = 1500)




