# By Florian Hartig. An extended commented version of this code as well as possible updates are available at http://florianhartig.github.io/LearningBayes/. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License


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
fit$model # inspect STAN model code that is produced by brms do the MCMC sampling

# TASK 1: look at brms package help to see what options exist to summarize the outputs
# TASK 2: look at the help to see how you would change the prior for the slope parameter


# Bayesian analysis with jags

library(rjags)

# 1) define the statistical model in JAGS code
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

# 2) Jags or STAN can't look into the R memory. Thus, you have to set up a list that contains all constant data that is used in the model
Data = list(y = airqualityCleaned$Ozone, x = airqualityCleaned$Temp, nobs = nrow(airqualityCleaned))

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

# Task 1: inspect all outputs - are results comparable to the lm and brms results?
# Task 2: how would you change the priors in this model? 
# Note that in Jags and Bugs, normal and many other distributions are parameterized via the precision = 1/sd^2.
# This means that a smaller value is a wider distribution. See also https://people.stat.sc.edu/hansont/stat740/jags_user_manual.pdf
# Task 3: how would you change the assumed residual distribution from a linear regression to a more outlier-robust regression with double-exponential residuals (i.e. exp(abs(x)))?


# Bayesian analysis with STAN

# Stan works practically identical to Jags just that the model structure is specified
# in several blocks, and that there is a block "data" where you have to specifify the
# type of all data that is provided to the model (this doesn't exist in Jags)

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


# Task 1: this STAN model doesn't include priors. Add some priors in the model section using code
# such as 
# alpha ~ normal(20, 5); beta ~ normal(-10, 5);
# You will need to know how the normal is parameterized 
# General Stan references are at https://mc-stan.org/users/documentation/
# For getting the function reference, see https://mc-stan.org/docs/functions-reference/index.html and 
# Look for the normal distribution


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

# Task 1: expect the help for what you can do with the output 
# Task 2: introduce a prior
# Task 3: change the likelihood to a double exponential, i.e. exp(abs(x)), as in Jags

