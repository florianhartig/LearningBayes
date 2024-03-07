Dat <- read.table('https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/LizardData.txt')
plot(Dat$Veg,Dat$Count)


# counts ~ poisson(predictedMean)
# predictedMean <- exp(a0 + a1 * veg + a2 * veg^2)

fit <- glm(Count ~ Veg + I(Veg^2), data = Dat, family = "poisson")
summary(fit)

testDispersion(fit)

testDispersion(fit,type = "PearsonChisq")

library(DHARMa)
res <- simulateResiduals(fit, plot = T)

library(glmmTMB)
fit <- glmmTMB(Count ~ Veg + I(Veg^2), data = Dat, family = "poisson", ziformula = ~ Veg + I(Veg^2) )
summary(fit)


library(rjags)
library(DHARMa)
library(BayesianTools)

curve(exp(3.38526  + 0.16947 * x -0.65658 * x^2), -2,2, add = T )

curve(20*plogis(-2.1304  + -0.9436 * x +  0.7052 * x^2), -2,2, add = T )

model = "
  model{

  # Likelihood 
  for(i in 1:nData){
    Lizzards[i] ~ dpois(lambda[i])
    lambda[i] <- predictedMean[i] * present[i] + 0.0000001
    
    predictedMean[i] <- exp(intercept + veg * vegData[i] + veg2 * vegData[i] * vegData[i])
    present[i] ~ dbin(pPresent, 1)
  }
  
  pPresent ~ dunif(0,1)
  
  intercept ~ dnorm(2.6, 0.0001)
  veg ~ dnorm(0, 0.01)
  veg2 ~ dnorm(0, 0.01)
  
  # Posterior predictive simulations 
  for(i in 1:nData){
    LizzardsSim[i] ~ dpois(predictedMean[i])
  }
 }
"

data = list(Lizzards= Dat$Count, 
            vegData = Dat$Veg, 
            nData = nrow(Dat))


jagsModel = jags.model(file = textConnection(model),
                       data = data,
                       n.chains = 3,
                       n.adapt = 5000)

# burnin 
update(jagsModel, n.iter=500)

Samples = coda.samples(jagsModel, 
                       variable.names = c("intercept", "veg", "veg2", "pPresent") ,
                       n.iter = 5000)
                       
summary(Samples)
plot(Samples)
gelman.diag(Samples)


PosteriorSimulations = coda.samples(jagsModel, 
                       variable.names = c("LizzardsSim") ,
                       n.iter = 5000)

x = getSample(PosteriorSimulations)
dim(x)
res <- createDHARMa(t(x), observedResponse = Dat$Count)
plot(res)
testDispersion(res)




# Model specification
model = "
  model{
  # Likelihood
  for(i in 1:n.dat){
    # poisson model p(y|lambda)
    y[i] ~ dpois(lambda[i])
    # logit link function
    log(lambda[i]) <- mu[i]
    # linear predictor on the log scale
    mu[i] <- alpha + beta.Veg*Veg[i] + beta.Veg2*Veg2[i]
    }
  # Priors
  alpha ~ dnorm(0,0.001)
  beta.Veg  ~ dnorm(0,0.001)
  beta.Veg2 ~ dnorm(0,0.001)
  }
 "

###########################################################
# Setting up the JAGS run:

# Set up a list that contains all the necessary data
Model.Data <- list(y = Dat$Count, n.dat = nrow(Dat),
                   Veg = Dat$Veg, Veg2 = Dat$Veg^2)

# Specify a function to generate inital values for the parameters
inits.fn <- function() list(alpha = rnorm(1,0,1),
                            beta.Veg = rnorm(1,0,1),
                            beta.Veg2 = rnorm(1,0,1))

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(model), data=Model.Data, 
                        inits = inits.fn, n.chains = 3, n.adapt= 5000)
# Specify parameters for which posterior samples are saved
para.names <- c('alpha','beta.Veg','beta.Veg2')

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel , variable.names = para.names, n.iter = 5000)

# Statistical summaries of the posterior distributions
summary(Samples)

# Plot MCMC samples
plot(Samples)

# Check convergence
gelman.diag(Samples)

# Correlation plot
correlationPlot(Samples)