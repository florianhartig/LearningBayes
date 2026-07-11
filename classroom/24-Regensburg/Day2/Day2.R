Dat = read.table("https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/Aspis_data.txt", stringsAsFactors = T)

# Inspect relationship between body mass and total body lenght
plot(Dat$TL, Dat$BM,
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')
# For the analysis we use log-transformed body masses
# and log-transformed and scaled total body lenght (TL)
plot(Dat$log_TL.sc, Dat$log_BM)

point.symbols <- c(f = 1, m = 4)
plot(Dat$TL, Dat$BM,
     pch = point.symbols[Dat$Sex],
     xlab = 'Total length [mm]',
     ylab = 'Body mass [g]')

library(rjags)

model ="
model{
  # Likelihood
  for(i in 1:n.dat){
    y[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha[Pop[i]] + beta.TL * TL[i] + beta.m * Sexm[i] + beta.mTL * TL[i] * Sexm[i]
    }
  
  for(p in 1:n.pop){
    alpha[p] ~ dnorm(alphaMean, tau.pop)
  }
  
  # Prior distributions
  alphaMean ~ dnorm(0,0.001) # tau = 0.001 corresponds to a sd = 31.62278
  beta.TL ~ dnorm(0,0.001)
  beta.m ~ dnorm(0,0.001)
  beta.mTL ~ dnorm(0,0.001)
  
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)
  tau.pop <- 1/(sigma.pop*sigma.pop)
  sigma.pop ~ dunif(0,100)
  }
"

Data = list(y = Dat$log_BM, 
            TL = Dat$log_TL.sc,
            Sexm = ifelse(Dat$Sex == "m", 1,0),
            Pop = Dat$Pop,
            n.pop = max(Dat$Pop),
            n.dat = nrow(Dat))

jagsModel <- jags.model(file = textConnection(model), data=Data, 
                        n.chains = 3, 
                        n.adapt= 5000)

para.names <- c("alpha","beta.TL", "beta.m", "beta.mTL", "sigma", "alphaMean", " sigma.pop" )

Samples <- coda.samples(jagsModel, variable.names = para.names, 
                        n.iter = 5000)

summary(Samples)
plot(Samples)

library(lme4)
LM <- lmer(log_BM ~ log_TL.sc * Sex + (1|Pop), data = Dat)
summary(LM)

brmsModel <- brm(log_BM ~ log_TL.sc * Sex + (1|Pop), data = Dat)
summary(brmsModel)


# Model selection 


set.seed(1)
dat = data.frame(matrix(runif(20000, -0.5,0.5), ncol = 100))
dat$y = rnorm(200)
dat$y = dat$y + rowSums(dat[,1:10]) 
# Preparing data list for Jags 
Data = list(y = dat$y, x = as.matrix(dat)[,1:100], i.max = nrow(dat))


fullModel = lm(y ~ . , data = dat)
summary(fullModel)

true = c(rep(1,10), rep(0,90))
estimated = coef(fullModel)[-1]
MSE = var(true - estimated)

plotEstimates <- function(estimates, ...){
  MSE = round(var(true - estimates), digits = 4)
  out <- barplot(estimates, las = 2, ylim = c(-0.5, 1.5), ...)
  text(60, 1, paste("MSE", MSE))
  lines(x = c(0,12), y = c(1,1), lwd = 4)
  lines(x = c(12,120), y = c(0,0), lwd = 4)
}
plotEstimates(estimated)

library(rjags)

modelCode0 = "model{
  # Likelihood
  for(i in 1:i.max){
    mu[i] <- inprod(a , x[i,]) + b
    y[i] ~ dnorm(mu[i],tau)
  }
  
  # Prior distributions
  for(i in 1:100){
    a[i] ~ dnorm(0,0.0001)
  }
  b ~ dnorm(0,0.0001) # usually no need and safer not to regularize intercept

  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
}
"

jagsModel0 <- jags.model(file= textConnection(modelCode0), 
                         data=Data, 
                         n.chains = 3)

Samples0 <- coda.samples(jagsModel0, 
                         variable.names = c("a","b","sigma"), 
                         n.iter = 5000)


x0<- summary(Samples0)
est0 <- x0$quantiles[1:100,3] 
plotEstimates(est0)



modelCode1 = "model{
  # Likelihood
  for(i in 1:i.max){
    mu[i] <- inprod(a , x[i,]) + b
    y[i] ~ dnorm(mu[i],tau)
  }
  
  # Prior distributions
  for(i in 1:100){
    a[i] ~ dnorm(0,0.5)
  }
  b ~ dnorm(0,0.0001) # usually no need and safer not to regularize intercept

  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
}
"


jagsModel1 <- jags.model(file= textConnection(modelCode1), 
                         data=Data, 
                         n.chains = 3)


Samples1 <- coda.samples(jagsModel1, 
                         variable.names = c("a","b","sigma"), 
                         n.iter = 5000)

#gelman.diag(Samples)
#summary(Samples)

x1<- summary(Samples1)
est1 <- x1$quantiles[1:100,3] 
plotEstimates(est1)


modelCode2 = "model{

  # Likelihood
  for(i in 1:i.max){
    mu[i] <- inprod(a , x[i,]) + b
    y[i] ~ dnorm(mu[i],tau)
  }
  
  # Prior distributions
  for(i in 1:100){
    a[i] ~ dnorm(0,tauShrinkage)
  }
  b ~ dnorm(0,0.001)

  tauShrinkage ~ dgamma(0.001, 0.001)
  sdShrinkage <- 1/sqrt(tauShrinkage)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
}
"

jagsModel2 <- jags.model(file= textConnection(modelCode2), data=Data, n.chains = 3)

Samples2 <- coda.samples(jagsModel2, 
                           variable.names = c("a","b","sigma", "sdShrinkage"), 
                           n.iter = 5000)

#gelman.diag(Samples)
summary(Samples2)

x2<- summary(Samples2)
est2 <- x2$quantiles[1:100,3] 
plotEstimates(est2)


modelCode3 = "model{
  # Likelihood
  for(i in 1:i.max){
    mu[i] <- inprod(a , x[i,]) + b
    y[i] ~ dnorm(mu[i],tau)
  }

  # Prior distributions
  pind ~ dbeta(5,5)
  for(j in 1:100){
    a_raw[j] ~ dnorm(0,0.01)
    ind[j] ~ dbern(pind)
    a[j] = ind[j] * a_raw[j]
  }
  b ~ dnorm(0,0.01)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
  }
"

jagsModel3 <- jags.model(file= textConnection(modelCode3), data=Data, n.chains = 3)

Samples3 <- coda.samples(jagsModel3, 
                         variable.names = c("a_raw", "ind","sigma"), 
                         n.iter = 5000)

x3<- summary(Samples3)

incl <- x3$quantiles[101:200,3] 
barplot(incl)

condEst <- x3$quantiles[1:100,3] 
plotEstimates(condEst)

est3 <- condEst * incl
plotEstimates(est3)

par(mfrow = c(2,2))
plotEstimates(est0, main = "Uninformative")
plotEstimates(est1, main = "Mildly regularizing")
plotEstimates(est2, main = "Adaptive shrinkage")
plotEstimates(est3, main = "Spike and slap")


dic0 = dic.samples(jagsModel0, n.iter = 5000)
dic0

dic1 = dic.samples(jagsModel1, n.iter = 5000)
dic1

dic2 = dic.samples(jagsModel2, n.iter = 5000)
dic2

dic3 = dic.samples(jagsModel3, n.iter = 5000)
dic3

# brms: I think they also have DIC now, but the programmers prefer WAIC and loocv

# Bayes factor, which corresponds to a LRT, is in the lecture notes!


Dat <- read.table('https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/LizardData.txt')
plot(Dat$Veg,Dat$Count)


Dat <- read.table('https://raw.githubusercontent.com/florianhartig/LearningBayes/master/data/LizardData.txt')
plot(Dat$Veg,Dat$Count)

fit <- glm(Count ~ Veg + I(Veg^2) , data = Dat, family = "poisson")
summary(fit)

library(effects)
plot(allEffects(fit, partial.residuals = T))

library(DHARMa)
res <- simulateResiduals(fit, plot = T)


library(rjags)
library(DHARMa)
library(BayesianTools)

# Model specification
model = "
  model{

  for(i in 1:n.dat){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu[i]
    mu[i] <- alpha + beta.Veg*Veg[i] + beta.Veg2*Veg2[i]
    }

  alpha ~ dnorm(0,0.001)
  beta.Veg  ~ dnorm(0,0.001)
  beta.Veg2 ~ dnorm(0,0.001)
  
  # Model predictions
    for(i in 1:n.dat){
      y.pred[i] ~ dpois(lambda.pred[i])
      log(lambda.pred[i]) <- mu.pred[i]
      mu.pred[i] <- alpha + beta.Veg*Veg[i] + beta.Veg2*Veg2[i]
    }
  }
 "

Model.Data <- list(y = Dat$Count, n.dat = nrow(Dat),
                   Veg = Dat$Veg, Veg2 = Dat$Veg^2)


jagsModel <- jags.model(file= textConnection(model), data=Model.Data, 
                         n.chains = 3, n.adapt= 5000)

# Specify parameters for which posterior samples are saved
para.names <- c('alpha','beta.Veg','beta.Veg2')

# Continue the MCMC runs with sampling
Samples <- coda.samples(jagsModel , variable.names = para.names, n.iter = 5000)

# Statistical summaries of the posterior distributions
summary(Samples)


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

# Create model checking plots
res = createDHARMa(simulatedResponse = t(Pred.Mat),
                   observedResponse = Dat$Count, 
                   fittedPredictedResponse = apply(Pred.Mat, 2, median),
                   integer = T)
plot(res)




library(EcoData)
library(rjags)

nobs = nrow(volcanoisland)

# imagine we had a very bad measurement devide for the altitude
volcanoisland$sAltitudeR = volcanoisland$sAltitude + rnorm(nobs)

plot(log(windObs) ~ sAltitude, data = volcanoisland)
fit = lm(log(windObs) ~ sAltitude, data = volcanoisland)
summary(fit)

abline(fit, col = "red")

fit = lm(log(windObs) ~ sAltitudeR, data = volcanoisland)
summary(fit)

abline(fit, col = "blue")

data = list(WindObs = log(volcanoisland$windObs),
            Altitude = volcanoisland$sAltitudeR,
            plot = as.numeric(volcanoisland$plot),
            nobs = nobs,
            nplots = length(unique(volcanoisland$plot)))


modelCode = "model{

  # Likelihood
  for(i in 1:nobs){
    WindObs[i] ~ dnorm(mu[i],tau)
    mu[i] <- AltitudeEffect*TrueAltitude[plot[i]]+ Intercept
    Altitude[i] ~ dnorm(TrueAltitude[plot[i]],tauMeasure)
  }
  
  # Prior distributions
  AltitudeEffect ~ dnorm(0,0.0001)
  Intercept ~ dnorm(0,0.0001)

  for(i in 1:nplots){
     TrueAltitude[i] ~ dnorm(0,0.0001)
  }

  # For scale parameters, normal choice is decaying
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)

  tauMeasure ~ dgamma(0.001, 0.001)
  sdMeasure <- 1/sqrt(tauMeasure)
}
"

jagsModel <- jags.model(file= textConnection(modelCode), data=data, n.chains = 3)

para.names <- c("AltitudeEffect","Intercept","sigma", "sdMeasure")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)


library(EcoData)
library(effects)
plot(lizardsObs ~ earth , data = volcanoisland)

fit<- glm(lizardsObs ~ earth + windObs , data = volcanoisland, family = binomial)
summary(fit)

plot(allEffects(fit))


Data = list(WindObs = log(volcanoisland$windObs),
            Altitude = volcanoisland$sAltitude[seq(1, 999, by = 10)],
            SoilPlot = unique(volcanoisland$earth),
            LizzardObs =volcanoisland$lizardsObs,
            plot = as.numeric(volcanoisland$plot),
            nobs = nrow(volcanoisland),
            nplots = length(unique(volcanoisland$plot)))


modelCode = "
model{

  # Likelihood
  for(i in 1:nobs){

    LizzardObs[i] ~ dbern(ObservationProb[i] *  LizzardTrue[plot[i]])
    logit(ObservationProb[i]) <- intO + windO * WindObs[i] 

  }

  for(i in 1:nplots){
    LizzardTrue[i] ~ dbern(LizzardSuitability[i])
    logit(LizzardSuitability[i]) <- intL+ SoilL*SoilPlot[i] + altL * Altitude[i]
  }

  # Prior distributions
  intO ~ dnorm(0,0.001)
  windO ~ dnorm(0,0.001)
  intL ~ dnorm(0,0.001)
  SoilL ~ dnorm(0,0.001)
  altL ~ dnorm(0,0.001)

  # posterior predictive simulations

  # Likelihood
  for(i in 1:nobs){
    LizzardObsSim[i] ~ dbern(ObservationProb[i] *  LizzardTrueSim[plot[i]])
  }
  for(i in 1:nplots){
    LizzardTrueSim[i] ~ dbern(LizzardSuitability[i])
  }

}
"

inits.fn <- function() list(LizzardTrue = rep(1,100))
jagsModel <- jags.model(file= textConnection(modelCode), inits = inits.fn, data=Data, n.chains = 3)

para.names <- c("intO","windO","intL", "SoilL", "altL")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)

rm(list=ls())
library(R2jags)


model = "
model {
# Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
  mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   r[t] ~ dnorm(mean.r, tau.proc)
   logN[t+1] <- logN[t] + r[t]
   }

# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(logN[t], tau.obs)
   }

# Population sizes on real scale
for (t in 1:T) {
   N[t] <- exp(logN[t])
   }
}
"


# House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm.counts <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
jags.data <- list(y = log(hm.counts), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

# Call JAGS from R (BRT 3 min)
hm.ssm <- jags(jags.data, inits, parameters, textConnection(model), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(hm.ssm, digits = 3)
plot(hm.ssm)


fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm.counts)
for (i in 1:n.years){
  fitted[i] <- mean(hm.ssm$BUGSoutput$sims.list$N[,i])
  lower[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N[,i], 0.025)
  upper[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N[,i], 0.975)}
m1 <- min(c(fitted, hm.counts, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm.counts, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm.counts, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)



# Population counts (from years 1 to 10)
y <- c(45, 48, 44, 59, 62, 62, 55, 51, 46, 42)

# Capture-recapture data (in m-array format, from years 1 to 10)
m <- matrix(c(11,  0,  0,  0,  0,  0,  0,  0,  0,  70,
              0, 12,  0,  1,  0,  0,  0,  0,  0,  52,
              0,  0, 15,  5,  1,  0,  0,  0,  0,  42,
              0,  0,  0,  8,  3,  0,  0,  0,  0,  51,
              0,  0,  0,  0,  4,  3,  0,  0,  0,  61,
              0,  0,  0,  0,  0, 12,  2,  3,  0,  66,
              0,  0,  0,  0,  0,  0, 16,  5,  0,  44,
              0,  0,  0,  0,  0,  0,  0, 12,  0,  46,
              0,  0,  0,  0,  0,  0,  0,  0, 11,  71,
              10,  2,  0,  0,  0,  0,  0,  0,  0,  13,
              0,  7,  0,  1,  0,  0,  0,  0,  0,  27,
              0,  0, 13,  2,  1,  1,  0,  0,  0,  14,
              0,  0,  0, 12,  2,  0,  0,  0,  0,  20,
              0,  0,  0,  0, 10,  2,  0,  0,  0,  21,
              0,  0,  0,  0,  0, 11,  2,  1,  1,  14,
              0,  0,  0,  0,  0,  0, 12,  0,  0,  18,
              0,  0,  0,  0,  0,  0,  0, 11,  1,  21,
              0,  0,  0,  0,  0,  0,  0,  0, 10,  26), ncol = 10, byrow = TRUE)

# Productivity data (from years 1 to 9)
J <- c(64, 132,  86, 154, 156, 134, 116, 106, 110)
R <- c(21, 28, 26, 38, 35, 33, 31, 30, 33) 

#library(rjags)
library(R2jags)

model = "
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#       1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
# Note that this part is different than in BUGS:
# 1. JAGS seems to be very sensitive to the choice of the prior distribution for the initial population sizes and of the choice of the observation model (priors)
# 2. Since the initial population sizes are used in binomial distributions, the numbers must be integers, otherwise JAGS does not run.
# The following specification seems to work very well and produces similar results as BUGS:

n1 ~ dnorm(25, tauy)T(0,)     # 1-year
nad ~ dnorm(25, tauy)T(0,)    # Adults
N1[1] <- round(n1)
Nad[1] <- round(nad)

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
      }
   
   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j  
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t

# 3.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t]*f[t]
   }
}
"

# Bundle data
jags.data <- list(m = m, y = y, J = J, R = R, nyears = dim(m)[2], r = rowSums(m))

# Initial values
initial<- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), sigma.y = runif(1, 0, 1), n1 = rpois(1, 30), nad = rpois(1, 30))}  
inits<-list(initial(),initial(),initial())
# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "sigma2.y", "lambda")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 25000
nc <- 3

# Call JAGS from R (BRT 2 min)
ipm <- jags(jags.data, inits, parameters, textConnection(model), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Produce Fig. 11-4
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:10){
  lower[i] <- quantile(ipm$BUGSoutput$sims.list$Ntot[,i], 0.025)
  upper[i] <- quantile(ipm$BUGSoutput$sims.list$Ntot[,i], 0.975)
}
plot(ipm$BUGSoutput$mean$Ntot, type = "b", ylim = c(35, 65), ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue", frame = F, cex = 1.5)
segments(1:10, lower, 1:10, upper, col = "blue")
points(y, type = "b", col = "black", pch = 16, lty = 2, cex = 1.5)
legend(x = 1, y = 65, legend = c("Counts", "Estimates"), pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")


library(piecewiseSEM) 
data(keeley)

library(lavaan)

library(lavaanPlot)

k_mod <- "
  rich ~ firesev + cover
  cover ~ firesev"

k_fit_lavaan <- sem(model = k_mod, data = keeley)
summary(k_fit_lavaan)

lavaanPlot(model=k_fit_lavaan, coefs = TRUE, sig = .05)

library(blavaan)

k_fit_blavaan = blavaan(model = k_mod, data = keeley,
                        auto.var=TRUE, auto.fix.first=TRUE,
                        auto.cov.lv.x=TRUE)

summary(k_fit_blavaan)
lavaanPlot(model=k_fit_blavaan, coefs = TRUE, sig = .05)


k_fit_psem <- psem(
  lm(rich ~ firesev + cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  data = keeley
)
summary(k_fit_psem)


rich_mod <- bf(rich ~ firesev + cover)
cover_mod <- bf(cover ~ firesev)

k_fit_brms <- brm(rich_mod +
                    cover_mod + 
                    set_rescor(FALSE), 
                  data=keeley,
                  cores=4, chains = 2)





