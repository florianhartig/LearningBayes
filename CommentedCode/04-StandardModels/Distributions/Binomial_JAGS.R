##################################################
# Bayesian analysis of binomially distributed data 
##################################################
#
#
# Estimate the probability parameter (pi) of a binomial distribution
# Example: seed predation from seed tray
# Felix May, 05.06. 2015

library(R2jags)

set.seed(12345)

#----------------------------
# 1. SIMULATE DATA
#----------------------------


# Set parameters
pi.true <- 0.7 # seed predation probability
N <- 10        # number of seeds per seed tray

# Simulate data
dat1 <- rbinom(n=30,size=N,prob=pi.true) # 30 seed trays in total

hist(dat1,breaks=seq(-0.5,10.5,by=1))
summary(dat1)
abline(v=pi.true*N,col="blue")

#-------------------------------
#2. Bayesian analysis using JAGS
#-------------------------------

# Model definition including likelihood and prior
modelCode <- "
   model{
      # Likelihood
      for(i in 1:n.max){
         y[i] ~ dbin(pi,N.total)
      }
      # Prior distributions
      pi ~ dbeta(1,1) # beta distribution
   }
"

# How does the prior look like?
x <- seq(0,1,0.001)
plot(x,dbeta(x,shape1=1,shape2=1),col=2,type="l",ylab="Prior density")


# Set up a list that contains all the necessary data 
Data = list(y = dat1, N.total=N,n.max = length(dat1))

# Set initial values for the Markov Chains
Inits = list(list("pi"=0.1),list("pi"=0.9),list("pi"=0.5))

# Start the MCMC sampling in JAGS
jags.fit <- jags(data=Data, 
                 inits=Inits, 
                 parameters.to.save=c("pi"), 
                 model.file=textConnection(modelCode),
                 n.chains=3, 
                 n.iter=6000,
                 n.burnin=1000,
                 n.thin=5)

# explore the model object
plot(jags.fit)
print(jags.fit)

# convert to coda mcmc object
jags.mcmc <- as.mcmc(jags.fit)

plot(jags.mcmc)  # traceplot and posterior distribution
summary(jags.mcmc) # numeric output
gelman.diag(jags.mcmc) # numeric check of convergence
HPDinterval(jags.mcmc) # highest posterior density intervals (=credible intervals)

#lump the three chains
jags.mcmc.lumped <- as.mcmc(rbind(jags.mcmc[[1]],jags.mcmc[[2]],jags.mcmc[[3]]))

hist(jags.mcmc.lumped[,"pi"],freq=F,ylab="Posterior density")
abline(v=HPDinterval(jags.mcmc.lumped[,"pi"]),col="red")

# Add the prior
lines(x,dbeta(x,shape1=1,shape2=1),col=3,lty=2)

#---------------------------------------------------------------------------
#Exercise: change the prior - check how posterior changes

# alternative prior --> very low seed predation probability
plot(x,dbeta(x,shape1=1,shape2=20),col=2,type="l",ylab="Prior density")

# Model definition including likelihood and prior
modelCode2 = "
   model{
      # Likelihood
      for(i in 1:n.max){
         y[i] ~ dbin(pi,N.total)
      }
      # Prior distributions
      pi ~ dbeta(1,20)
   }
"

jags.fit2 <- jags(data=Data, 
                 inits=Inits, 
                 parameters.to.save=c("pi"), 
                 model.file=textConnection(modelCode2),
                 n.chains=3, 
                 n.iter=6000,
                 n.burnin=1000,
                 n.thin=5)

jags.mcmc2 <- as.mcmc(jags.fit2)

summary(jags.mcmc)
summary(jags.mcmc2)
HPDinterval(jags.mcmc2)

jags.mcmc.lumped2 <- as.mcmc(rbind(jags.mcmc2[[1]],jags.mcmc2[[2]],jags.mcmc2[[3]]))

# plot both priors and posteriors in one plot
plot(density(jags.mcmc.lumped[,"pi"]),ylab="Posterior density",
     xlim=c(0,1),col=3,main="")
lines(x,dbeta(x,shape1=1,shape2=1),col=3,lty=2)

lines(density(jags.mcmc.lumped2[,"pi"]),col="red")
lines(x,dbeta(x,shape1=1,shape2=20),col="red",lty=2)




#  **Copyright, reuse and updates**: copyright belongs to author(s) (see author statement at the top of the file). Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License

#Sourcecode and potential future updates available at http://florianhartig.github.io/LearningBayes/ (follow the link under code, and then navigate through the topics to find the location of the file)

