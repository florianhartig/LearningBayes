#' ---
#' output: html_document
#' editor_options:
#'   chunk_output_type: console
#' ---
#' 
## ---- include=FALSE--------------------------------------------------------------------------------
set.seed(42)

#' 
#' # Generalized linear mixed models
#' 
#' ::: callout-note
#' In this chapter, we will discuss 
#' 
#' :::
#' 
#' ## Poisson regression
#' 
## --------------------------------------------------------------------------------------------------
Dat <- read.table('https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/LizardData.txt')
plot(Dat$Veg,Dat$Count)

#' 
#' Standard frequentist GLM
#' 
## --------------------------------------------------------------------------------------------------
fit <- glm(Count ~ Veg, data = Dat, family = "poisson")
summary(fit)

#' 
## --------------------------------------------------------------------------------------------------
library(rjags)
library(DHARMa)
library(BayesianTools)

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

#' 
#' ## Adding posterior predictions and residual checks
#' 
## --------------------------------------------------------------------------------------------------
model ="
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

  # Model predictions
  for(i in 1:n.pred){
    y.pred[i] ~ dpois(lambda.pred[i])
    log(lambda.pred[i]) <- mu.pred[i]
    mu.pred[i] <- alpha + beta.Veg*Veg.pred[i] + beta.Veg2*Veg2.pred[i]
    }
  }
 "

###########################################################
# Setting up the JAGS run:

# Set up a list that contains all the necessary data
# Note that for prediction (later used for model checking) 
# we here use the original predictor variables.
Model.Data <- list(y = Dat$Count, n.dat = nrow(Dat),
                   Veg = Dat$Veg, Veg2 = Dat$Veg^2,
                   Veg.pred = Dat$Veg, Veg2.pred = Dat$Veg^2,
                   n.pred = nrow(Dat))

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
#############################################################
# Sample simulated posterior for lizard counts (y.pred)
Pred.Samples <- coda.samples(jagsModel, 
                             variable.names = "y.pred", 
                             n.iter = 5000)

# Transform mcmc.list object to a matrix
Pred.Mat <- as.matrix(Pred.Samples)

# Plot Model predictions against data
Pred.Q <- apply(Pred.Mat,2,quantile,prob=c(0.05,0.5,0.95))
plot(Dat$Veg, Dat$Count)
ord <- order(Dat$Veg)
lines(Dat$Veg[ord], Pred.Q['50%',ord],col='blue',lwd=2)
lines(Dat$Veg[ord], Pred.Q['5%',ord],col='blue')
lines(Dat$Veg[ord], Pred.Q['95%',ord],col='blue')

###########################################################
# Model checking with DHARMa

# Create model checking plots
res = createDHARMa(simulatedResponse = t(Pred.Mat),
                   observedResponse = Dat$Count, 
                   fittedPredictedResponse = apply(Pred.Mat, 2, median),
                   integer = T)
plot(res)
###########################################################

#' 
#' ## Adding an overdispersion term
#' 
## --------------------------------------------------------------------------------------------------
model = "
  model{
  # Likelihood
  for(i in 1:n.dat){
    # poisson model p(y|lambda)
    y[i] ~ dpois(lambda[i])
    # logit link function
    log(lambda[i]) <- mu[i] + eps[i]
    # linear predictor on the log scale
    mu[i] <- alpha + beta.Veg*Veg[i] + beta.Veg2*Veg2[i]
    # overdispersion error terms
    eps[i] ~ dnorm(0,tau.eps) 
    }
  # Priors
  alpha ~ dnorm(0,0.001)
  beta.Veg  ~ dnorm(0,0.001)
  beta.Veg2 ~ dnorm(0,0.001)
  tau.eps ~ dgamma(0.001,0.001)

  # Model predictions
  for(i in 1:n.pred){
    y.pred[i] ~ dpois(lambda.pred[i])
    log(lambda.pred[i]) <- mu.pred[i] + eps.pred[i]
    mu.pred[i] <- alpha + beta.Veg*Veg.pred[i] + beta.Veg2*Veg2.pred[i]
    eps.pred[i] ~ dnorm(0, tau.eps)
  }
  }

"
###########################################################
# Setting up the JAGS run:

# Set up a list that contains all the necessary data
# Note that for prediction (later used for model checking) 
# we here use the original predictor variables.
Model.Data <- list(y = Dat$Count, n.dat = nrow(Dat),
                   Veg = Dat$Veg, Veg2 = Dat$Veg^2,
                   Veg.pred = Dat$Veg, Veg2.pred = Dat$Veg^2,
                   n.pred = nrow(Dat))

# Specify a function to generate inital values for the parameters
inits.fn <- function() list(alpha = rnorm(1,0,1),
                            beta.Veg = rnorm(1,0,1),
                            beta.Veg2 = rnorm(1,0,1),
                            tau.eps = 1
                            )

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(model), data=Model.Data, 
                        inits = inits.fn, n.chains = 3, n.adapt= 5000)
# Specify parameters for which posterior samples are saved
para.names <- c('alpha','beta.Veg','beta.Veg2',
                'tau.eps')
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
#############################################################
# Sample simulated posterior for lizard counts (y.pred)
Pred.Samples <- coda.samples(jagsModel, 
                             variable.names = "y.pred", 
                             n.iter = 5000)

# Transform mcmc.list object to a matrix
Pred.Mat <- as.matrix(Pred.Samples)

# Plot Model predictions against data
Pred.Q <- apply(Pred.Mat,2,quantile,prob=c(0.05,0.5,0.95))
plot(Dat$Veg, Dat$Count)
ord <- order(Dat$Veg)
lines(Dat$Veg[ord], Pred.Q['50%',ord],col='blue',lwd=2)
lines(Dat$Veg[ord], Pred.Q['5%',ord],col='blue')
lines(Dat$Veg[ord], Pred.Q['95%',ord],col='blue')

###########################################################
# Model checking with DHARMa

# Create model checking plots
res = createDHARMa(simulatedResponse = t(Pred.Mat),
                   observedResponse = Dat$Count, 
                   fittedPredictedResponse = apply(Pred.Mat, 2, median),
                   integer = T)
plot(res)
###########################################################



#' 
#' ## Adding zero-inflation
#' 
## --------------------------------------------------------------------------------------------------
model = "
  model{
  # Likelihood
  for(i in 1:n.dat){
    # poisson model p(y|lambda)
    y[i] ~ dpois(lambda.eff[i])
    # effective mean abundance
    lambda.eff[i] <- lambda[i] * Inc[i]
    # binary variable to indicate occupancy
    Inc[i] ~ dbin(p.Inc,1)
    # logit link function
    log(lambda[i]) <- mu[i] + eps[i]
    # linear predictor on the log scale
    mu[i] <- alpha + beta.Veg*Veg[i] + beta.Veg2*Veg2[i]
    # overdispersion error terms
    eps[i] ~ dnorm(0,tau.eps) 
    }
  # Priors
  alpha ~ dnorm(0,0.001)
  beta.Veg  ~ dnorm(0,0.001)
  beta.Veg2 ~ dnorm(0,0.001)
  tau.eps ~ dgamma(0.001,0.001)
  p.Inc ~ dbeta(1,1)

  # Model predictions
  for(i in 1:n.pred){
    y.pred[i] ~ dpois(lambda.eff.pred[i])
    lambda.eff.pred[i] <- lambda.pred[i] * Inc.pred[i]
    Inc.pred[i] ~ dbin(p.Inc, 1)
    log(lambda.pred[i]) <- mu.pred[i] + eps.pred[i]
    mu.pred[i] <- alpha + beta.Veg*Veg.pred[i] + beta.Veg2*Veg2.pred[i]
    eps.pred[i] ~ dnorm(0, tau.eps)
  }
  }
"

###########################################################
# Setting up the JAGS run:

# Set up a list that contains all the necessary data
# Note that for prediction (later used for model checking) 
# we here use the original predictor variables.
Model.Data <- list(y = Dat$Count, n.dat = nrow(Dat),
                   Veg = Dat$Veg, Veg2 = Dat$Veg^2,
                   Veg.pred = Dat$Veg, Veg2.pred = Dat$Veg^2,
                   n.pred = nrow(Dat))

# Specify a function to generate inital values for the parameters
inits.fn <- function() list(alpha = rnorm(1,0,1),
                            beta.Veg = rnorm(1,0,1),
                            beta.Veg2 = rnorm(1,0,1),
                            tau.eps = 1,
                            p.Inc = rbeta(1,1,1),
                            Inc = rep(1,nrow(Dat))
                            )

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file= textConnection(model), data=Model.Data, 
                        inits = inits.fn, n.chains = 3, n.adapt= 5000)
# Specify parameters for which posterior samples are saved
para.names <- c('alpha','beta.Veg','beta.Veg2',
                'tau.eps','p.Inc')
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
#############################################################
# Sample simulated posterior for lizard counts (y.pred)
Pred.Samples <- coda.samples(jagsModel, 
                             variable.names = "y.pred", 
                             n.iter = 5000)

# Transform mcmc.list object to a matrix
Pred.Mat <- as.matrix(Pred.Samples)

# Plot Model predictions against data
Pred.Q <- apply(Pred.Mat,2,quantile,prob=c(0.05,0.5,0.95))
plot(Dat$Veg, Dat$Count)
ord <- order(Dat$Veg)
lines(Dat$Veg[ord], Pred.Q['50%',ord],col='blue',lwd=2)
lines(Dat$Veg[ord], Pred.Q['5%',ord],col='blue')
lines(Dat$Veg[ord], Pred.Q['95%',ord],col='blue')

###########################################################
# Model checking with DHARMa

# Create model checking plots
res = createDHARMa(simulatedResponse = t(Pred.Mat),
                   observedResponse = Dat$Count, 
                   fittedPredictedResponse = apply(Pred.Mat, 2, median),
                   integer = T)
plot(res)

