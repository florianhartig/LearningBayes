#' ---
#' output: html_document
#' editor_options:
#'   chunk_output_type: console
#' ---
#' 
## ---- include=FALSE--------------------------------------------------------------------------------
set.seed(42)

#' 
#' # Posterior estimation
#' 
#' ::: callout-note
#' In this chapter, we will discuss 
#' 
#' :::
#' 
#' ## Fitting a linear regression with different MCMCs
#' 
#' We will use the dataset airquality, just removing NAs and scaling all variables for convenience
#' 
## --------------------------------------------------------------------------------------------------
airqualityCleaned = airquality[complete.cases(airquality),]
airqualityCleaned = data.frame(scale(airqualityCleaned))
plot(Ozone ~ Temp, data = airqualityCleaned)

#' 
#' ### Frequentist inference
#' 
## ---- message = F----------------------------------------------------------------------------------
fit <- lm(Ozone ~ Temp, data = airqualityCleaned)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))
par(mfrow = c(2,2))
plot(fit) # residuals

#' 
#' ### Bayesian analysis with brms
#' 
## ---- message = F, echo=FALSE, results = 'hide'----------------------------------------------------
library(brms)
fit <- brm(Ozone ~ Temp, data = airqualityCleaned)

#' 
## --------------------------------------------------------------------------------------------------
summary(fit)
plot(fit, ask = FALSE)
plot(conditional_effects(fit), ask = FALSE)
fit$model # model that is actually fit via 
pp_check(fit) # residual checks

#' 
#' ### Bayesian analysis with JAGS
#' 
#' The general approach in JAGS is to
#' 
#' 1.  Set up a list that contains all the necessary data
#' 2.  Write the model as a string in the JAGS specific BUGS dialect
#' 3.  Compile the model and run the MCMC for an adaptation (burn-in) phase
#' 
## --------------------------------------------------------------------------------------------------
library(rjags)
Data = list(y = airqualityCleaned$Ozone, x = airqualityCleaned$Temp, nobs = nrow(airqualityCleaned))

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

# Specify a function to generate inital values for the parameters 
# (optional, if not provided, will start with the mean of the prior )
inits.fn <- function() list(a = rnorm(1), b = rnorm(1), tau = 1/runif(1,1,100))

# sets up the model
jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3)

# MCMC sample from model
Samples <- coda.samples(jagsModel, variable.names = c("a","b","sigma"), n.iter = 5000)

# Plot the mcmc chain and the posterior sample 
plot(Samples)
summary(Samples)

#' 
#' ### Bayesian analysis with STAN
#' 
#' Approach is identical to JAGS just that we have to define all variables in the section data
#' 
## ---- message = F, results = 'hide'----------------------------------------------------------------
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

#' 
## --------------------------------------------------------------------------------------------------
print(fit)
plot(fit)
rstan::traceplot(fit)

#' 
#' ### Bayesian analysis via BayesianTools
#' 
#' Here, we don't use a model specification language, but just write out the likelihood as an standard R function. The same can be done for the prior. For simplicity, in this case I just used flat priors using the lower / upper arguments.
#' 
## ----  message = F, results='hide'-----------------------------------------------------------------
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

#' 
## --------------------------------------------------------------------------------------------------
plot(out)
summary(out, start = 1000)

#' 
#' ## Checking and expecting the results
#' 
#' Running the sampler again
#' 
## --------------------------------------------------------------------------------------------------
Samples <- coda.samples(jagsModel, variable.names = c("a","b","sigma"), n.iter = 5000)

#' 
#' ### Convergence checks
#' 
#' Except for details in the syntax, the following is more or less the same for all samplers.
#' 
#' First thing should always be convergence checks. Visual look at the trace plots,
#' 
## --------------------------------------------------------------------------------------------------
plot(Samples)

#' 
#' We want to look at
#' 
#' 1.  Convergence to the right parameter area (seems immediate, else you will see a slow move of the parameters in the traceplot). You should set burn-in after you have converged to the right area
#' 2.  Mixing: low autocorrelation in the chain after convergence to target area (seems excellent in this case)
#' 
#' Further convergence checks should be done AFTER removing burn-in
#' 
## --------------------------------------------------------------------------------------------------
coda::acfplot(Samples)

#' 
#' Formal convergence diagnostics via
#' 
## --------------------------------------------------------------------------------------------------
coda::gelman.diag(Samples)
coda::gelman.plot(Samples)

#' 
#' No fixed rule but typically people require univariate psrf \< 1.05 or \< 1.1 and multivariate psrf \< 1.1 or 1.2
#' 
#' ::: callout-caution
#' Note that the msrf rule was made for estimating the mean / median. If you want to estimate more unstable statistics, e.g. higher quantiles or other values such as the MAP or the DIC (see section on model selection), you may have to run the MCMC chain much longer to get stable outputs.
#' 
## ---- message = F, results='hide'------------------------------------------------------------------

  library(BayesianTools)
  bayesianSetup <- createBayesianSetup(likelihood = testDensityNormal, 
                                       prior = createUniformPrior(lower = -10,
                                                                  upper = 10))
  out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 3000))

#' 
#' The plotDiagnostics function in package BT shows us how statistics develop over time
#' 
## --------------------------------------------------------------------------------------------------
plotDiagnostic(out)

#' :::
#' 
#' ### Summary Table
#' 
## --------------------------------------------------------------------------------------------------
summary(Samples)

#' 
#' Highest Posterior Density intervals
#' 
## --------------------------------------------------------------------------------------------------
HPDinterval(Samples)

#' 
#' ### Plots
#' 
#' Marginal plots show the parameter distribution (these were also created in the standard coda traceplots)
#' 
## --------------------------------------------------------------------------------------------------
BayesianTools::marginalPlot(Samples)

#' 
#' Pair correlation plots show 2nd order correlations
#' 
## --------------------------------------------------------------------------------------------------
# coda
coda::crosscorr.plot(Samples)
#BayesianTools
correlationPlot(Samples)

#' 
#' ### Posterior predictive distribution
#' 
## --------------------------------------------------------------------------------------------------
dat = as.data.frame(Data)[,1:2]
dat = dat[order(dat$x),]
# raw data
plot(dat[,2], dat[,1])

# extract 1000 parameters from posterior from package BayesianTools
x = getSample(Samples, start = 300)
pred = x[,2] + dat[,2] %o% x[,1] 
lines(dat[,2], apply(pred, 1, median))
lines(dat[,2], apply(pred, 1, quantile, probs = 0.2), 
      lty = 2, col = "red")
lines(dat[,2], apply(pred, 1, quantile, probs = 0.8), 
      lty = 2, col = "red")

# alternative: plot all 1000 predictions in transparent color
plot(dat[,2], dat[,1])
for(i in 1:nrow(x)) lines(dat[,2], pred[,i], col = "#0000EE03")

# important point - so far, we have plotted 
# in frequentist, this is know as the confidence vs. the prediction distribution
# in the second case, we add th


pred = x[,2] + dat[,2] %o% x[,1] 
for(i in 1:nrow(x))  {
  pred[,i] = pred[,i] + rnorm(length(pred[,i]), 0, sd = x[i,3])
}

plot(dat[,2], dat[,1])
lines(dat[,2], apply(pred, 1, median))
lines(dat[,2], apply(pred, 1, quantile, probs = 0.2), lty = 2, col = "red")
lines(dat[,2], apply(pred, 1, quantile, probs = 0.8), lty = 2, col = "red")

#alternative plotting
polygon(x = c(dat[,2], rev(dat[,2])), 
        y = c(apply(pred, 1, quantile, probs = 0.2), 
              rev(apply(pred, 1, quantile, probs = 0.8))), 
        col = "#EE000020")

#' 
#' ## Playing around with the pipeline
#' 
#' ### Prior choice
#' 
#' Priors are not scale-free. What that means: dnorm(0,0.0001) might not be an uninformative prior, if the data scale is extremely small so that you might expect huge effect sizes - scaling all variables makes sure we have a good intuition of what "uninformative means".
#' 
#' **Task:** play with the following minimal script for a linear regression to understand how scaling parameter affects priors and thus posterior shapes. In particular, change
#' 
#' -   Multiply Ozone by 1000000 -\> will push sd estimates high
#' -   Multiply Temp by 0.0000001 -\> will push parameter estimates high
#' 
#' Then compare Bayesian parameter estimates and their uncertainty to Bayesian estimates. How would you have to change the priors to fix this problem and keep them uninformative?
#' 
#' **Task: 2** implement mildly informative priors as well as strong shrinkage priors in the regression. Question to discuss: should you put the shrinkage also in the intercept? Why should you center center variables if you include a shrinkage prior on the intercept?
#' 
## --------------------------------------------------------------------------------------------------
library(rjags)

dat = airquality[complete.cases(airquality),] 
# scaling happens here - change 
dat$Ozone = as.vector(scale(dat$Ozone))
dat$Temp = as.vector(scale(dat$Temp)) 


Data = list(y = dat$Ozone, 
            x = dat$Temp, 
            i.max = nrow(dat))

# Model
modelCode = "
model{

  # Likelihood
  for(i in 1:i.max){
    mu[i] <- Temp*x[i]+ intercept
    y[i] ~ dnorm(mu[i],tau)
  }

  # Prior distributions
  
  # For location parameters, typical choice is wide normal
  intercept ~ dnorm(0,0.0001)
  Temp ~ dnorm(0,0.0001)

  # For scale parameters, typical choice is decaying
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)

}
"

# Specify a function to generate inital values for the parameters (optional, if not provided, will start with the mean of the prior )
inits.fn <- function() list(a = rnorm(1), b = rnorm(1), 
                            tau = 1/runif(1,1,100))

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3)

# Run a bit to have a burn-in
update(jagsModel, n.iter = 1000)


# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = c("intercept","Temp","sigma"), n.iter = 5000)


# Bayesian results
summary(Samples)

# MCMC results
fit <- lm(Ozone ~ Temp, data = dat)
summary(fit)

#' 
#' ### Missing data
#' 
#' In the analysis above, we removed missing data. What happens if you are leaving the missing data in in a Jags model? Try it out and discuss what happens
