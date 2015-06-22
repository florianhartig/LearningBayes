##################################################
# Bayesian analysis of a simple linear regression 
##################################################
#
# Simplified example following: 
# Kery (2010) Introduction to WinBUGS for Ecologists
# Chapter 11, pages 141 - 150
#
# Relationship between body mass and body length of snakes
#
# Felix May, 05.06. 2015

library(R2jags)

set.seed(12345)

#----------------------------
# 1. SIMULATE DATA
#----------------------------

n1 <- 30                         # sample size
Length <- sort(runif(n1,45,70))  # explanatory variable

a.true <- 30                     # intercept
b.true <- 2.5                    # slope
sigma.true <- 10                 # residual variance

mean.true <- a.true+b.true*Length
Mass <- rnorm(n1,mean=mean.true,sd=sigma.true)

plot(Mass~Length)
lines(mean.true~Length,col="blue")

#-------------------------------
#2. Bayesian analysis using JAGS
#-------------------------------

# Model definition including likelihood and prior

modelCode <- "
   model{
      # Likelihood
      for(i in 1:n.max){
         y[i] ~ dnorm(mu[i],tau) # mean and precision (=1/variance)
         mu[i] <- a + b*x[i]
      }
      # Prior distributions
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)
      tau ~ dunif(0,100)
      sigma <- sqrt(1/tau)
   }"

# Set up a list that contains all the necessary data
Data = list(y = Mass, x = Length, n.max = length(Mass))

# Specify a function to generate inital values for the parameters
inits.fn <- function() list(a = rnorm(1), b = rnorm(1),  tau=runif(1,min=0,max=1))

# Start the MCMC sampling in JAGS
jags.fit <- jags(data=Data, 
                 inits=inits.fn, 
                 parameters.to.save=c("a","b","sigma"), 
                 model.file=textConnection(modelCode),
                 n.chains=3, 
                 n.iter=12000,
                 n.burnin=2000,
                 n.thin=10,
                 DIC=F)

# explore the model object
plot(jags.fit)
print(jags.fit)

# convert to coda mcmc object
jags.mcmc <- as.mcmc(jags.fit)

plot(jags.mcmc)    # plot MCMC traces and posterior distributions
acfplot(jags.mcmc) # check for autocorrelation in the Markov chains

summary(jags.mcmc)
gelman.diag(jags.mcmc) # convergence diagnostis
HPDinterval(jags.mcmc) # highes posterior density intervals (= credible intervals)

#lump the three chains for final analysis
jags.mcmc.lumped <- as.mcmc(rbind(jags.mcmc[[1]],
                                  jags.mcmc[[2]],
                                  jags.mcmc[[3]]))

# check correlation among variables
pairs(as.data.frame(jags.mcmc.lumped))
cor(as.data.frame(jags.mcmc.lumped))

# center length to remove the correlation
cLength <- Length-mean(Length)
plot(Mass~cLength)

# adapt data for JAGS
Data2 = list(y = Mass, x = cLength, n.max = length(Mass))

jags.fit2 <- jags(data=Data2, 
                 inits=inits.fn, 
                 parameters.to.save=c("a","b","sigma"), 
                 model.file=textConnection(modelCode),
                 n.chains=3, 
                 n.iter=12000,
                 n.burnin=2000,
                 n.thin=10,
                 DIC=F)

jags.mcmc2 <- as.mcmc(jags.fit2)
plot(jags.mcmc2)
summary(jags.mcmc2)
gelman.diag(jags.mcmc2)
HPDinterval(jags.mcmc2)

#lump the three chains
jags.mcmc.lumped2 <- as.mcmc(rbind(jags.mcmc2[[1]],
                                   jags.mcmc2[[2]],
                                   jags.mcmc2[[3]]))

# check correlation among variables
pairs(as.data.frame(jags.mcmc.lumped2))
cor(as.data.frame(jags.mcmc.lumped2))
# --> much better after centering


#-------------------------------
#3. Predictions and uncertainty
#-------------------------------

plot(Mass~cLength,pch=16,ylim=c(120,240))

# create predictions with all posterior samples
pred1 <- matrix(NA,nrow=nrow(jags.mcmc.lumped2),ncol=length(Mass))
for (i in 1:nrow(pred1)){
   pred1[i,] <- jags.mcmc.lumped2[i,"a"] + cLength*jags.mcmc.lumped2[i,"b"]
   #lines(cLength,pred1[i,],col="gray")
}
points(Mass~cLength,pch=16)
# Prediction for posterior mean
meanPred <- mean(jags.mcmc.lumped2[,"a"]) + 
            cLength*mean(jags.mcmc.lumped2[,"b"])
lines(meanPred~cLength,lwd=2,col="blue")

# 95% predictive credible interval for the MEAN BODY MASS
lower1 <- apply(pred1,MARGIN=2,quantile,prob=0.025)
upper1 <- apply(pred1,MARGIN=2,quantile,prob=0.975)

lines(cLength,lower1,col="red",lwd=2)
lines(cLength,upper1,col="red",lwd=2)

# 95% credible intervals for the BODY MASS OF ONE SNAKE
# version 1
pred2 <- matrix(NA,nrow=nrow(jags.mcmc.lumped2),ncol=length(Mass))
for (i in 1:nrow(pred2)){
   pred2[i,] <- jags.mcmc.lumped2[i,"a"] + 
                cLength*jags.mcmc.lumped2[i,"b"] +
                rnorm(length(Mass),mean=0,sd=jags.mcmc.lumped2[i,"sigma"])
}

lower2 <- apply(pred2,MARGIN=2,quantile,prob=0.025)
upper2 <- apply(pred2,MARGIN=2,quantile,prob=0.975)

lines(cLength,lower2,col="green",lwd=2,lty=1)
lines(cLength,upper2,col="green",lwd=2,lty=1)

# let JAGS calculate the predicted values

modelCode2 <- "
   model{
      # Likelihood
      for(i in 1:n.max){
         y[i] ~ dnorm(mu[i],tau) # mean and precision (=1/variance)
         ypred[i] ~ dnorm(mu[i],tau)
         mu[i] <- a + b*x[i]
      }
      # Prior distributions
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)
      tau ~ dunif(0,100)
      sigma <- sqrt(1/tau)
}"


# Start the MCMC sampling in JAGS
jags.fit3 <- jags(data=Data2, 
                 inits=inits.fn, 
                 parameters.to.save=c("a","b","sigma","ypred"), 
                 model.file=textConnection(modelCode2),
                 n.chains=3, 
                 n.iter=12000,
                 n.burnin=2000,
                 n.thin=10,
                 DIC=F)

plot(jags.fit3)

jags.mcmc3 <- as.mcmc(jags.fit3)

gelman.diag(jags.mcmc3)


#lump the three chains
jags.mcmc.lumped3 <- as.mcmc(rbind(jags.mcmc3[[1]],
                                   jags.mcmc3[[2]],
                                   jags.mcmc3[[3]]))
dim(jags.mcmc.lumped3)
head(jags.mcmc.lumped3)

pred.mat <- jags.mcmc.lumped3[,4:33]
head(pred.mat)

index.order <- order(as.numeric(sort(as.character(1:30))))

lower3 <- apply(pred.mat[,index.order],MARGIN=2,quantile,prob=0.025)
upper3 <- apply(pred.mat[,index.order],MARGIN=2,quantile,prob=0.975)

plot(Mass~cLength,pch=16,ylim=c(120,240))

lines(cLength,lower3,col="green",lwd=2,lty=1)
lines(cLength,upper3,col="green",lwd=2,lty=1)

#  **Copyright, reuse and updates**: copyright belongs to author(s) (see author statement at the top of the file). Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License

#Sourcecode and potential future updates available at http://florianhartig.github.io/LearningBayes/ (follow the link under code, and then navigate through the topics to find the location of the file)

