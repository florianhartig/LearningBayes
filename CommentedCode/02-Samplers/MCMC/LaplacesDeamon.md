# Different Sampling Algorithms using Laplaces Demon
FlorianHartig  
25 Apr 2015  

Laplaces Demon is (or was, I'm not sure) a package that implements various MCMC algorithms in R, as well as a few other convenience functions. 

The (was) refers to the fact that the package was pulled from CRAM, and it was shortly also removed from GitHub. I have a copy on my GH account, which you should be able to install via


```r
#install.packages("devtools")
library(devtools)
install_github("florianhartig/LaplacesDemon")
```

If everything worked, you should be able to load the package


```r
library(LaplacesDemon)
```

## What can you do with the package

Basically, the package provides a large number of general-purpose MCMCs that expect a target function (likelihood / posterior) and than sample from this. A nice side story of this is that the package has a large number of examples for explicitly defining likelihood functions (as opposed to defining models in bugs-style). You get these examples via



```r
vignette("Examples")
```

There is also a nice introduction to Bayesian inference, via


```r
vignette("BayesianInference")
```

The main thing of interest, however, are the samplers that are available. A description is provided here 


```r
vignette("LaplacesDemonTutorial")
```

I am not sure if all of the algorithms / settings work well and are stable. See, e.g.

* http://wiekvoet.blogspot.de/2014/10/tuning-laplacesdemon.html


Here's the example from the package. 



```r
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#########################  Data List Preparation  #########################
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
     beta <- rnorm(Data$J)
     sigma <- runif(1)
     return(c(beta, sigma))
     }
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     ### Log-Prior
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Likelihood
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }
#library(compiler)
#Model <- cmpfun(Model) #Consider byte-compiling for more speed

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

####################  Automated Factor Slice Sampler  #####################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))
```

```
## 
## Laplace's Demon was called on Fri Nov 20 22:05:27 2015
## 
## Performing initial checks...
## Algorithm: Automated Factor Slice Sampler 
## 
## Laplace's Demon is beginning to update...
## 
## Eigendecomposition will occur every 50 iterations.
## 
## Iteration: 100,   Proposal: Multivariate,   LP: -62.1
## Iteration: 200,   Proposal: Multivariate,   LP: -60.7
## Iteration: 300,   Proposal: Multivariate,   LP: -62.3
## Iteration: 400,   Proposal: Multivariate,   LP: -62.8
## Iteration: 500,   Proposal: Multivariate,   LP: -60.3
## Iteration: 600,   Proposal: Multivariate,   LP: -60.1
## Iteration: 700,   Proposal: Multivariate,   LP: -64.5
## Iteration: 800,   Proposal: Multivariate,   LP: -61.7
## Iteration: 900,   Proposal: Multivariate,   LP: -67.6
## Iteration: 1000,   Proposal: Multivariate,   LP: -63.5
## 
## Assessing Stationarity
## Assessing Thinning and ESS
## Creating Summaries
## Estimating Log of the Marginal Likelihood
## Creating Output
## 
## Laplace's Demon has finished.
```

```r
Fit
```

```
## Call:
## LaplacesDemon(Model = Model, Data = MyData, Initial.Values = Initial.Values, 
##     Covar = NULL, Iterations = 1000, Status = 100, Thinning = 1, 
##     Algorithm = "AFSS", Specs = list(A = Inf, B = NULL, m = 100, 
##         n = 0, w = 1))
## 
## Acceptance Rate: 1
## Algorithm: Automated Factor Slice Sampler
## Covariance Matrix: (NOT SHOWN HERE; diagonal shown instead)
##    beta[1]    beta[2]    beta[3]    beta[4]      sigma 
## 0.03170632 0.11218865 0.10047765 0.27938662 1.13324982 
## 
## Covariance (Diagonal) History: (NOT SHOWN HERE)
## Deviance Information Criterion (DIC):
##          All Stationary
## Dbar  83.444     83.444
## pD    69.224     69.224
## DIC  152.669    152.669
## Initial Values:
## [1]  0.75331105  2.01435467 -0.35513446  2.02816784  0.01331584
## 
## Iterations: 1000
## Log(Marginal Likelihood): -42.75201
## Minutes of run-time: 0.04
## Model: (NOT SHOWN HERE)
## Monitor: (NOT SHOWN HERE)
## Parameters (Number of): 5
## Posterior1: (NOT SHOWN HERE)
## Posterior2: (NOT SHOWN HERE)
## Recommended Burn-In of Thinned Samples: 0
## Recommended Burn-In of Un-thinned Samples: 0
## Recommended Thinning: 4
## Specs: (NOT SHOWN HERE)
## Status is displayed every 100 iterations
## Summary1: (SHOWN BELOW)
## Summary2: (SHOWN BELOW)
## Thinned Samples: 1000
## Thinning: 1
## 
## 
## Summary of All Samples
##                 Mean         SD        MCSE       ESS           LB
## beta[1]    5.0431344  0.1761881 0.005928101 1000.0000   4.81532956
## beta[2]    0.6018762  0.3329698 0.021594155  354.8674   0.08953703
## beta[3]    1.1947535  0.3161568 0.017020539  490.5966   0.60108044
## beta[4]    0.8679975  0.5193960 0.026468908  650.4790   0.22890230
## sigma      0.7627413  1.0383075 0.061934415  455.7384   0.55954025
## Deviance  83.4444566 11.7663986 1.218725684  164.4786  78.02817991
## LP       -62.9001446  5.9084103 0.610929049  164.7706 -67.00061593
##               Median          UB
## beta[1]    5.0419212   5.2832715
## beta[2]    0.5864291   1.1804153
## beta[3]    1.1919701   1.7666925
## beta[4]    0.8756034   1.5277576
## sigma      0.7043232   0.9284134
## Deviance  81.9571591  91.6441175
## LP       -62.1551251 -60.1904380
## 
## 
## Summary of Stationary Samples
##                 Mean         SD        MCSE       ESS           LB
## beta[1]    5.0431344  0.1761881 0.005928101 1000.0000   4.81532956
## beta[2]    0.6018762  0.3329698 0.021594155  354.8674   0.08953703
## beta[3]    1.1947535  0.3161568 0.017020539  490.5966   0.60108044
## beta[4]    0.8679975  0.5193960 0.026468908  650.4790   0.22890230
## sigma      0.7627413  1.0383075 0.061934415  455.7384   0.55954025
## Deviance  83.4444566 11.7663986 1.218725684  164.4786  78.02817991
## LP       -62.9001446  5.9084103 0.610929049  164.7706 -67.00061593
##               Median          UB
## beta[1]    5.0419212   5.2832715
## beta[2]    0.5864291   1.1804153
## beta[3]    1.1919701   1.7666925
## beta[4]    0.8756034   1.5277576
## sigma      0.7043232   0.9284134
## Deviance  81.9571591  91.6441175
## LP       -62.1551251 -60.1904380
```

```r
print(Fit)
```

```
## Call:
## LaplacesDemon(Model = Model, Data = MyData, Initial.Values = Initial.Values, 
##     Covar = NULL, Iterations = 1000, Status = 100, Thinning = 1, 
##     Algorithm = "AFSS", Specs = list(A = Inf, B = NULL, m = 100, 
##         n = 0, w = 1))
## 
## Acceptance Rate: 1
## Algorithm: Automated Factor Slice Sampler
## Covariance Matrix: (NOT SHOWN HERE; diagonal shown instead)
##    beta[1]    beta[2]    beta[3]    beta[4]      sigma 
## 0.03170632 0.11218865 0.10047765 0.27938662 1.13324982 
## 
## Covariance (Diagonal) History: (NOT SHOWN HERE)
## Deviance Information Criterion (DIC):
##          All Stationary
## Dbar  83.444     83.444
## pD    69.224     69.224
## DIC  152.669    152.669
## Initial Values:
## [1]  0.75331105  2.01435467 -0.35513446  2.02816784  0.01331584
## 
## Iterations: 1000
## Log(Marginal Likelihood): -42.75201
## Minutes of run-time: 0.04
## Model: (NOT SHOWN HERE)
## Monitor: (NOT SHOWN HERE)
## Parameters (Number of): 5
## Posterior1: (NOT SHOWN HERE)
## Posterior2: (NOT SHOWN HERE)
## Recommended Burn-In of Thinned Samples: 0
## Recommended Burn-In of Un-thinned Samples: 0
## Recommended Thinning: 4
## Specs: (NOT SHOWN HERE)
## Status is displayed every 100 iterations
## Summary1: (SHOWN BELOW)
## Summary2: (SHOWN BELOW)
## Thinned Samples: 1000
## Thinning: 1
## 
## 
## Summary of All Samples
##                 Mean         SD        MCSE       ESS           LB
## beta[1]    5.0431344  0.1761881 0.005928101 1000.0000   4.81532956
## beta[2]    0.6018762  0.3329698 0.021594155  354.8674   0.08953703
## beta[3]    1.1947535  0.3161568 0.017020539  490.5966   0.60108044
## beta[4]    0.8679975  0.5193960 0.026468908  650.4790   0.22890230
## sigma      0.7627413  1.0383075 0.061934415  455.7384   0.55954025
## Deviance  83.4444566 11.7663986 1.218725684  164.4786  78.02817991
## LP       -62.9001446  5.9084103 0.610929049  164.7706 -67.00061593
##               Median          UB
## beta[1]    5.0419212   5.2832715
## beta[2]    0.5864291   1.1804153
## beta[3]    1.1919701   1.7666925
## beta[4]    0.8756034   1.5277576
## sigma      0.7043232   0.9284134
## Deviance  81.9571591  91.6441175
## LP       -62.1551251 -60.1904380
## 
## 
## Summary of Stationary Samples
##                 Mean         SD        MCSE       ESS           LB
## beta[1]    5.0431344  0.1761881 0.005928101 1000.0000   4.81532956
## beta[2]    0.6018762  0.3329698 0.021594155  354.8674   0.08953703
## beta[3]    1.1947535  0.3161568 0.017020539  490.5966   0.60108044
## beta[4]    0.8679975  0.5193960 0.026468908  650.4790   0.22890230
## sigma      0.7627413  1.0383075 0.061934415  455.7384   0.55954025
## Deviance  83.4444566 11.7663986 1.218725684  164.4786  78.02817991
## LP       -62.9001446  5.9084103 0.610929049  164.7706 -67.00061593
##               Median          UB
## beta[1]    5.0419212   5.2832715
## beta[2]    0.5864291   1.1804153
## beta[3]    1.1919701   1.7666925
## beta[4]    0.8756034   1.5277576
## sigma      0.7043232   0.9284134
## Deviance  81.9571591  91.6441175
## LP       -62.1551251 -60.1904380
```

```r
#Consort(Fit)
#plot(BMK.Diagnostic(Fit))
#PosteriorChecks(Fit)
#caterpillar.plot(Fit, Parms="beta")
#BurnIn <- Fit$Rec.BurnIn.Thinned
#plot(Fit, BurnIn, MyData, PDF=FALSE)
#Pred <- predict(Fit, Model, MyData, CPUs=1)
#summary(Pred, Discrep="Chi-Square")
#plot(Pred, Style="Covariates", Data=MyData)
#plot(Pred, Style="Density", Rows=1:9)
#plot(Pred, Style="ECDF")
#plot(Pred, Style="Fitted")
#plot(Pred, Style="Jarque-Bera")
#plot(Pred, Style="Predictive Quantiles")
#plot(Pred, Style="Residual Density")
#plot(Pred, Style="Residuals")
#Levene.Test(Pred)
#Importance(Fit, Model, MyData, Discrep="Chi-Square")

#############  Adaptive Directional Metropolis-within-Gibbs  ##############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="ADMG", Specs=list(n=0, Periodicity=50))

########################  Adaptive Griddy-Gibbs  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AGG", Specs=list(Grid=GaussHermiteQuadRule(3)$nodes,
#     dparm=NULL, smax=Inf, CPUs=1, Packages=NULL, Dyn.libs=NULL))

##################  Adaptive Hamiltonian Monte Carlo  #####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AHMC", Specs=list(epsilon=0.02, L=2, m=NULL,
#     Periodicity=10))

##########################  Adaptive Metropolis  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

###################  Adaptive Metropolis-within-Gibbs  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMWG", Specs=list(B=NULL, n=0, Periodicity=50))

######################  Adaptive-Mixture Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMM", Specs=list(Adaptive=500, B=NULL, n=0,
#     Periodicity=10, w=0.05))

###################  Affine-Invariant Ensemble Sampler  ###################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AIES", Specs=list(Nc=2*length(Initial.Values), Z=NULL,
#     beta=2, CPUs=1, Packages=NULL, Dyn.libs=NULL))

#################  Componentwise Hit-And-Run Metropolis  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="CHARM", Specs=NULL)

###########  Componentwise Hit-And-Run (Adaptive) Metropolis  #############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="CHARM", Specs=list(alpha.star=0.44))

#################  Delayed Rejection Adaptive Metropolis  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRAM", Specs=list(Adaptive=500, Periodicity=10))

#####################  Delayed Rejection Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRM", Specs=NULL)

##################  Differential Evolution Markov Chain  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DEMC", Specs=list(Nc=3, Z=NULL, gamma=NULL, w=0.1))

#######################  Elliptical Slice Sampler  ########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="ESS", Specs=list(B=NULL))

#############################  Gibbs Sampler  #############################
### NOTE: Unlike the other samplers, Gibbs requires specifying a
### function (FC) that draws from full conditionals.
#FC <- function(parm, Data)
#     {
#     ### Parameters
#     beta <- parm[Data$pos.beta]
#     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
#     sigma2 <- sigma*sigma
#     ### Hyperparameters
#     betamu <- rep(0,length(beta))
#     betaprec <- diag(length(beta))/1000
#     ### Update beta
#     XX <- crossprod(Data$X)
#     Xy <- crossprod(Data$X, Data$y)
#     IR <- backsolve(chol(XX/sigma2 + betaprec), diag(length(beta)))
#     btilde <- crossprod(t(IR)) %*% (Xy/sigma2 + betaprec %*% betamu)
#     beta <- btilde + IR %*% rnorm(length(beta))
#     return(c(beta,sigma))
#     }
##library(compiler)
##FC <- cmpfun(FC) #Consider byte-compiling for more speed
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="Gibbs", Specs=list(FC=FC, MWG=pos.sigma))

#############################  Griddy-Gibbs  ##############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="GG", Specs=list(Grid=seq(from=-0.1, to=0.1, len=5),
#     dparm=NULL, CPUs=1, Packages=NULL, Dyn.libs=NULL))

#######################  Hamiltonian Monte Carlo  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HMC", Specs=list(epsilon=0.001, L=2, m=NULL))

#############  Hamiltonian Monte Carlo with Dual-Averaging  ###############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="HMCDA", Specs=list(A=500, delta=0.65, epsilon=NULL,
#     Lmax=1000, lambda=0.1))

#######################  Hit-And-Run Metropolis  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HARM", Specs=NULL)

##################  Hit-And-Run (Adaptive) Metropolis  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HARM", Specs=list(alpha.star=0.234, B=NULL))

########################  Independence Metropolis  ########################
### Note: the mu and Covar arguments are populated from a previous Laplace
### Approximation.
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=Fit$Covar, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="IM",
#     Specs=list(mu=Fit$Summary1[1:length(Initial.Values),1]))

#########################  Interchain Adaptation  #########################
#Initial.Values <- rbind(Initial.Values, GIV(Model, MyData, PGF=TRUE))
#Fit <- LaplacesDemon.hpc(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="INCA", Specs=list(Adaptive=500, Periodicity=10),
#     LogFile="MyLog", Chains=2, CPUs=2, Type="PSOCK", Packages=NULL,
#     Dyn.libs=NULL)

################  Metropolis-Adjusted Langevin Algorithm  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=1,
#          delta=1, epsilon=c(1e-6,1e-7)))

#############  Metropolis-Coupled Markov Chain Monte Carlo  ###############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MCMCMC", Specs=list(lambda=1, CPUs=2, Packages=NULL,
#     Dyn.libs=NULL))

#######################  Metropolis-within-Gibbs  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MWG", Specs=list(B=NULL))

########################  Multiple-Try Metropolis  ########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MTM", Specs=list(K=4, CPUs=1, Packages=NULL, Dyn.libs=NULL))

##########################  No-U-Turn Sampler  ############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="NUTS", Specs=list(A=500, delta=0.6, epsilon=NULL,
#     Lmax=Inf))

#################  Oblique Hyperrectangle Slice Sampler  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="OHSS", Specs=list(A=Inf, n=0))

#####################  Preconditioned Crank-Nicolson  #####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="pCN", Specs=list(beta=0.1))

######################  Robust Adaptive Metropolis  #######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RAM", Specs=list(alpha.star=0.234, B=NULL, Dist="N",
#     gamma=0.66, n=0))

###################  Random Dive Metropolis-Hastings  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RDMH", Specs=NULL)

##########################  Refractive Sampler  ###########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="Refractive", Specs=list(Adaptive=1, m=2, w=0.1, r=1.3))

###########################  Reversible-Jump  #############################
#bin.n <- J-1
#bin.p <- 0.2
#parm.p <- c(1, rep(1/(J-1),(J-1)), 1)
#selectable <- c(0, rep(1,J-1), 0)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RJ", Specs=list(bin.n=bin.n, bin.p=bin.p,
#          parm.p=parm.p, selectable=selectable,
#          selected=c(0,rep(1,J-1),0)))

########################  Random-Walk Metropolis  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RWM", Specs=NULL)

########################  Reflective Slice Sampler  #######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RSS", Specs=list(m=5, w=1e-5))

##############  Sequential Adaptive Metropolis-within-Gibbs  ##############
#NOTE: The SAMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SAMWG", Specs=list(Dyn=Dyn, Periodicity=50))

##################  Sequential Metropolis-within-Gibbs  ###################
#NOTE: The SMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SMWG", Specs=list(Dyn=Dyn))

#############################  Slice Sampler  #############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=1, Thinning=1,
#     Algorithm="Slice", Specs=list(B=NULL, Bounds=c(-Inf,Inf), m=100,
#     Type="Continuous", w=1))

#################  Stochastic Gradient Langevin Dynamics  #################
#NOTE: The Data and Model functions must be coded differently for SGLD.
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=10, Thinning=10,
#     Algorithm="SGLD", Specs=list(epsilon=1e-4, file="X.csv", Nr=1e4,
#     Nc=6, size=10))

###################  Tempered Hamiltonian Monte Carlo  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="THMC", Specs=list(epsilon=0.001, L=2, m=NULL,
#     Temperature=2))

###############################  t-walk  #################################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5))

#################  Univariate Eigenvector Slice Sampler  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="UESS", Specs=list(A=Inf, B=NULL, m=100, n=0))

##########  Updating Sequential Adaptive Metropolis-within-Gibbs  #########
#NOTE: The USAMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USAMWG", Specs=list(Dyn=Dyn, Periodicity=50, Fit=Fit,
#     Begin=T.m))

##############  Updating Sequential Metropolis-within-Gibbs  ##############
#NOTE: The USMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USMWG", Specs=list(Dyn=Dyn, Fit=Fit, Begin=T.m))

#End
```



---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License





