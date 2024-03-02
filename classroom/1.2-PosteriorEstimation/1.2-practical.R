# note: we will continue to improve this model in 1.3, see comments in this lecture if you use this as a  template

# dataset airquality

airqualityCleaned = airquality[complete.cases(airquality),]
airqualityCleaned = data.frame(scale(airqualityCleaned))
plot(Ozone ~ Temp, data = airqualityCleaned)


# Frequentist inference

fit <- lm(Ozone ~ Temp, data = airqualityCleaned)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))
par(mfrow = c(2,2))
plot(fit) # residuals


# Bayesian regression with brms

library(brms)
fit <- brm(Ozone ~ Temp, data = airqualityCleaned)
summary(fit)
plot(fit, ask = FALSE)
plot(conditional_effects(fit), ask = FALSE)
fit$model # model that is actually fit via 
pp_check(fit) # residual check

# Bayesian analysis with jags

library(rjags)


# 1) Set up a list that contains all the necessary data (here, including parameters of the prior distribution)

Data = list(y = airqualityCleaned$Ozone, x = airqualityCleaned$Temp, nobs = nrow(airqualityCleaned))

# 2) set up the model for JAGS
modelCode = "
model{

  # Likelihood
  for(i in 1:nobs){
    mu[i] <- a*x[i]+ b
    y[i] ~ dnorm(mu[i],tau) # dnorm in jags parameterizes via precision = 1/sd^2
  }

  # Prior distributions
  
  # For location parameters, normal choice is wide normal
  a ~ dnorm(0,0.0001)
  b ~ dnorm(0,0.0001)

  # For scale parameters, normal choice is decaying
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)

}
"

# 3) Specify a function to generate inital values for the parameters 
# (optional, if not provided, will start with the mean of the prior )
inits.fn <- function() list(a = rnorm(1), b = rnorm(1), tau = 1/runif(1,1,100))

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3)

# Specify parameters for which posterior samples are saved
para.names <- c("a","b","sigma")

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

# Plot the mcmc chain and the posterior sample 
plot(Samples)

summary(Samples)


# Bayesian analysis with STAN

library(rstan)

stanmodelcode <- "
  data {
    int<lower=0> N;
    vector[N] x;
    vector[N] y;
  }
  parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
  }
  model {
    y ~ normal(alpha + beta * x, sigma);
  }
"

dat = list(y = airqualityCleaned$Ozone, x = airqualityCleaned$Temp, N = nrow(airqualityCleaned))

fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 2012, chains = 3, verbose = TRUE,
            sample_file = file.path(tempdir(), 'norm.csv')) 
print(fit)
plot(fit)
rstan::traceplot(fit)


# Bayesian analysis via BayesianTools 

library(BayesianTools)

likelihood <- function(par){
  a0 = par[1]
  a1 = par[2]
  sigma <- par[3]  
  logLikel = sum(dnorm(a0 + a1 * airqualityCleaned$Temp  - airqualityCleaned$Ozone , sd = sigma, log = T))
  return(logLikel)
}

setup <- createBayesianSetup(likelihood = likelihood, lower = c(-10,-10,0.01), upper = c(10,10,10), names = c("a0", "a1", "sigma"))

out <- runMCMC(setup)
plot(out)
summary(out, start = 1000)



