# By Florian Hartig. An extended commented version of this code as well as possible updates are available at http://florianhartig.github.io/LearningBayes/. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License


# Continue with the previous model

rm(list = ls())
dat = airquality[complete.cases(airquality),]
# dat = data.frame(scale(dat))
dat = dat[order(dat$Temp),] # order so that we can later make more convenient line plots

library(rjags)


modelCode = "
  model{
  
    # Likelihood
    for(i in 1:i.max){
      mu[i] <- Temp*x[i]+ intercept
      y[i] ~ dnorm(mu[i],tau)
    }
  
    # Prior distributions
    
    # For location parameters, typical choice is wide normal
    Temp ~ dnorm(0,0.0001)
    intercept ~ dnorm(0,0.0001)
  
    # For scale parameters, typical choice is decaying
    tau ~ dgamma(0.001, 0.001)
    sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)
  
  }
"

Data = list(y = dat$Ozone, x = dat$Temp, i.max = nrow(dat))

jagsModel <- jags.model(file= textConnection(modelCode), 
                        data=Data,
                        n.chains = 3)

Samples <- coda.samples(jagsModel, 
                        variable.names = c("Temp","intercept","sigma"), 
                        n.iter = 5000)


###### Check convergence ########

gelman.diag(Samples)
gelman.plot(Samples)

autocorr.plot(Samples)

rejectionRate(Samples) 
effectiveSize(Samples)

cumuplot(Samples, probs=c(0.0025,0.5,0.9975))

###### Summarize Outputs ########

# basic summary 
summary(Samples)
plot(Samples)

library(BayesianTools)
marginalPlot(Samples)

# Correlations
coda::crosscorr.plot(Samples)
BayesianTools::correlationPlot(Samples, start = 300)

# Highest Posterior Density intervals
coda::HPDinterval(Samples)

# Compare to Frequentist analysis
fit <- lm(Ozone ~ Temp, data = dat)
summary(fit)

###### Posterior Predictive Distribution ########

# In the following example, we forward parameter uncertainties to predictions
# We do this by drawing from the posterior of the parameters and then 
# calculate predictions for each possible parameter combination
# note that this could also be done directly in jags, by observing 
# the model predictions mu[i] in coda.samples(). We will often use this option 
# later in the book 


# extract 1000 parameters from posterior using getSamples() from package 
# BayesianTools. We discard burn-in using the argument start = 300
# then we calculate the predictions for the parameters 
# result is a matrix with 100 prediction lines 
x = data.frame(getSample(Samples, start = 300))
pred = x$intercept +  x$Temp %o% dat$Temp 

# plotting the distributions and uncertainties 
par(mfrow = c(1,2))
plot(Ozone ~ Temp, data = dat, main = ("80% Credible interval"))
lines(dat$Temp, apply(pred, 2, median))
lines(dat$Temp, apply(pred, 2, quantile, probs = 0.1), 
      lty = 2, col = "red")
lines(dat$Temp, apply(pred, 2, quantile, probs = 0.9), 
      lty = 2, col = "red")

# optional: plot all 1000 predictions in transparent color
for(i in 1:nrow(x)) lines(dat$Temp, pred[i,], col = "#0000EE03")

# So far, we have plotted credible intervals,
# i.e. the interval of the true values. We can also plot prediction 
# intervals, i.e. the interval where we expect next observations to land 
# these can be compared to the data to look for outliers etc. 

pred = x$intercept +  x$Temp %o% dat$Temp 
for(i in 1:nrow(x))  {
  pred[i,] = pred[i,] + rnorm(length(pred[i,]), 0, sd = x$sigma[i])
}

plot(Ozone ~ Temp, data = dat, main = ("80% Prediction interval"))
lines(dat$Temp, apply(pred, 2, median))
lines(dat$Temp, apply(pred, 2, quantile, probs = 0.1), lty = 2, col = "red")
lines(dat$Temp, apply(pred, 2, quantile, probs = 0.9), lty = 2, col = "red")

#alternative plotting
polygon(x = c(dat$Temp, rev(dat$Temp)), 
        y = c(apply(pred, 2, quantile, probs = 0.1), 
              rev(apply(pred, 2, quantile, probs = 0.9))), 
        col = "#EE000020")



# Compare to Frequentist analysis
library(effects)
plot(allEffects(fit, partial.residuals = T))


# Task: after running the entire script, scale / center the data and observe what 
# happens to the paramter estimates and their correlation - can you explain why? 
# use dat = data.frame(scale(dat))








