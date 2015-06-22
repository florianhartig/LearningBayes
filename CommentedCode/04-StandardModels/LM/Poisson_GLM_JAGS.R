#####################################
# Bayesian analysis of a Poisson GLM 
#####################################
#
# Simplified example following: 
# Kery (2010) Introduction to WinBUGS for Ecologists
# Chapter 15, pages 193 - 202
#
# Relationship between wing length and mite infection of dragonflies
#
# Felix May, 05.06. 2015

library(R2jags)

set.seed(12345)

#----------------------------
# 1. SIMULATE DATA
#----------------------------

# choose parameter values
n1 <- 100                                    # sample size
Winglength <- sort(runif(n1,5,8))            # explanatory variable
cWinglength <- Winglength - mean(Winglength) # center explanatory variable for easier model convergence

a.true <- 2.5    # intercept
b.true <- - 1.1  # slope

eta1 <- a.true+b.true*cWinglength     #linear predictor
lambda.true <- exp(eta1)              # inverse log-link function 

Mites <- rpois(n1,lambda=lambda.true) # add Poisson error

plot(Mites~cWinglength)
lines(lambda.true ~ cWinglength,col="blue")

#-------------------------------
#2. Bayesian analysis using JAGS
#-------------------------------

# Model definition including likelihood and prior
modelCode <- "
   model{
      # Likelihood
      for(i in 1:n.max){
         y[i] ~ dpois(lambda[i])  # poisson error distribution
         lambda[i] <- exp(eta[i]) # inverse link function
         eta[i] <- a + b*x[i]     # linear predictor
      }
      # Prior distributions
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)
   }"

# Set up a list that contains all the necessary data
Data = list(y = Mites, x = cWinglength, n.max = length(Mites))

# Specify a function to generate inital values for the parameters
inits.fn <- function() list(a = rnorm(1,mean=0,sd=10), b = rnorm(1,mean=0,sd=10))

# Start the MCMC sampling in JAGS
jags.fit <- jags(data=Data, 
                 inits=inits.fn,
                 parameters.to.save=c("a","b"), 
                 model.file=textConnection(modelCode),
                 n.chains=3, 
                 n.iter=25000,
                 n.burnin=5000,
                 n.thin=20,
                 DIC=T)

# explore the model object
plot(jags.fit)
print(jags.fit)

# convert to coda mcmc object
jags.mcmc <- as.mcmc(jags.fit)

# explore mcmc object
windows()
plot(jags.mcmc)     # check convergence and posterior distributions
acfplot(jags.mcmc)  # check autocorrelation in chains

summary(jags.mcmc)
gelman.diag(jags.mcmc)  # check convergence
HPDinterval(jags.mcmc)  # posterior credible intervals

# strange first values in the chains despite burnin
# no idea why this happens
head(jags.mcmc[[1]])
head(jags.mcmc[[2]])
head(jags.mcmc[[3]])


#lump the three chains and remove first values
jags.mcmc.lumped <- as.mcmc(rbind(jags.mcmc[[1]][-1,],
                                  jags.mcmc[[2]][-1,],
                                  jags.mcmc[[3]][-1,]))

plot(jags.mcmc.lumped)

# check correlation among variables
pairs(as.data.frame(jags.mcmc.lumped))
cor(as.data.frame(jags.mcmc.lumped))


#-------------------------------
#3. Predictions and uncertainty
#-------------------------------

plot(Mites~cWinglength)

# generate predictions with each posterior sample
pred1 <- matrix(NA,nrow=nrow(jags.mcmc.lumped),ncol=length(Mites))
for (i in 1:nrow(pred1))
   pred1[i,] <- exp(jags.mcmc.lumped[i,"a"] + cWinglength*jags.mcmc.lumped[i,"b"])

# predictions with posterior mean
lambdaPred <- exp(mean(jags.mcmc.lumped[,"a"]) + 
                  cWinglength*mean(jags.mcmc.lumped[,"b"]))
lines(lambdaPred~cWinglength,lwd=2,col="red")

# credible interval for mean mite load
lower1 <- apply(pred1,MARGIN=2,quantile,prob=0.025)
upper1 <- apply(pred1,MARGIN=2,quantile,prob=0.975)

lines(cWinglength,lower1,col="red",lwd=1)
lines(cWinglength,upper1,col="red",lwd=1)

# credible interval for mite load of single individuals
# generate predictions with each posterior sample
pred2 <- matrix(NA,nrow=nrow(jags.mcmc.lumped),ncol=length(Mites))
for (i in 1:nrow(pred1)){
   lambda.pred <- exp(jags.mcmc.lumped[i,"a"] + cWinglength*jags.mcmc.lumped[i,"b"])
   pred2[i,] <- rpois(length(cWinglength),lambda=lambda.pred)                      
}

lower2 <- apply(pred2,MARGIN=2,quantile,prob=0.025)
upper2 <- apply(pred2,MARGIN=2,quantile,prob=0.975)

lines(cWinglength,lower2,col="green",lwd=1,lty=2)
lines(cWinglength,upper2,col="green",lwd=1,lty=2)

# add known true relationship
lines(lambda.true ~ cWinglength,col="blue")

#  **Copyright, reuse and updates**: copyright belongs to author(s) (see author statement at the top of the file). Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License

#Sourcecode and potential future updates available at http://florianhartig.github.io/LearningBayes/ (follow the link under code, and then navigate through the topics to find the location of the file)


