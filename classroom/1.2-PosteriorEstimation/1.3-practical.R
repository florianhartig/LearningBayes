
# note: this example improves the previous script 1.2 - start with this!

# dataset airquality

plot(Ozone ~ Temp, data = airquality)


# Frequentist analysis

fit <- lm(Ozone ~ Temp, data = airquality)
summary(fit)
library(effects)
plot(allEffects(fit, partial.residuals = T))
par(mfrow = c(2,2))
plot(fit)


# Bayesian analysis

library(rjags)


# 1) Set up a list that contains all the necessary data 
# (here, including parameters of the prior distribution)

dat = airquality[complete.cases(airquality),] # remove NAs, Jags doesn't remove them automatically

# scale the data, because your priors are not scale-free 
# what that means: dnorm(0,0.0001) might not be an 
# uninformative prior, if the data scale is extremely 
# small so that you might expect huge effect sizes - 
# scaling all variables makes sure we have a good 
# intuition of what "uninformative means"

Data = list(y = as.vector(scale(dat$Ozone)), x = as.vector(scale(dat$Temp)), i.max = nrow(dat))



# 2) Model definition exactly how we created our data 
modelCode = "
model{

  # Likelihood
  for(i in 1:i.max){
    mu[i] <- a*x[i]+ b
    y[i] ~ dnorm(mu[i],tau)
  }

  # Prior distributions
  
  # For location parameters, typical choice is wide normal
  a ~ dnorm(0,0.0001)
  b ~ dnorm(0,0.0001)

  # For scale parameters, typical choice is decaying
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)

}
"


# 3) Specify a function to generate inital values for the parameters 
# (optional, if not provided, will start with the mean of the prior )
inits.fn <- function() list(a = rnorm(1), b = rnorm(1), 
                            tau = 1/runif(1,1,100))

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(modelCode), data=Data, init = inits.fn, n.chains = 3)

# Specify parameters for which posterior samples are saved
para.names <- c("a","b","sigma")

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel, variable.names = para.names, 
                        n.iter = 5000)

# Plot the mcmc chain and the posterior sample
plot(Samples)

summary(Samples)