#' ---
#' output: html_document
#' editor_options:
#'   chunk_output_type: console
#' ---
#' 
## ---- include=FALSE-----------------------------------------------------------------------
set.seed(42)

#' 
#' # Bayesian causal models
#' 
#' This chapter assumes you are aware of causal inference, causal graphs and structural equation models. For a general intro on these topics, see my lecture notes on advanced regression models, section causal inference [here](https://theoreticalecology.github.io/AdvancedRegressionModels/3B-CausalInference.html).
#' 
#' In this chapter, I will show you a few options to fit SEMs with Bayesian methods. We will use the keeley data provided in the package piecewiseSEM
#' 
## ---- message=F---------------------------------------------------------------------------
library(piecewiseSEM) 
data(keeley)

#' 
#' ## Lavaan and blavaan
#' 
#' Laavan is the most popular package to fit multivariate normal SEMs. The advantage of lavaan (over piecewiseSEM which we will use later) is that it can easily include latent factors in the analysis (which we will not use in this example). The downside of lavaan is that it can't handle non-Gaussian responses.
#' 
#' Here the standard lavaan model:
#' 
## -----------------------------------------------------------------------------------------
library(lavaan)
library(lavaanPlot)

k_mod <- "
  rich ~ firesev + cover
  cover ~ firesev"

k_fit_lavaan <- sem(model = k_mod, data = keeley)

lavaanPlot(model=k_fit_lavaan, coefs = TRUE, sig = .05)


#' 
#' Same as Bayesian fit using blavaan, which is a laavan style interface to STAN
#' 
## -----------------------------------------------------------------------------------------
library(blavaan)

k_fit_blavaan = blavaan(model = k_mod, data = keeley,
                auto.var=TRUE, auto.fix.first=TRUE,
                auto.cov.lv.x=TRUE)


summary(k_fit_blavaan)

lavaanPlot(model=k_fit_blavaan, coefs = TRUE, sig = .05)


#' 
#' ## Piecewise SEMs
#' 
#' The piecewiseSEM allows you to create a causal model from a set of standard GLMMs. Advantage is that it can handle standard GLMMs. Disadvantage is that it can't handle latent variables
#' 
## -----------------------------------------------------------------------------------------
k_fit_psem <- psem(
  lm(rich ~ firesev + cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  data = keeley
)
summary(k_fit_psem)


#' 
#' Piecewise SEMs can pretty much in the same way be code in brm:
#' 
## -----------------------------------------------------------------------------------------
library(brms)

rich_mod <- bf(rich ~ firesev + cover)
cover_mod <- bf(cover ~ firesev)

k_fit_brms <- brm(rich_mod +
                  cover_mod + 
                  set_rescor(FALSE), 
                data=keeley,
                cores=4, chains = 2)

plot(k_fit_brms)


#' 
#' ## Freestyle JAGS model
#' 
#' If you want latent variables and non-Gaussian responses, you will have to code by hand. Here an example, which is effectively a piecewise SEM bout could be extended to include latent variables.
#' 
#' TODO
#' 
## ---- eval = F----------------------------------------------------------------------------
## library(EcoData)
## library(effects)
## 
## 
## islandPsem <- psem(
##   lm(windObs ~ sAltitude, data = volcanoisland),
##   glm(lizardsObs ~ sAltitude , family = binomial, data = volcanoisland),
##   glm(beetles ~ windObs + lizardsObs,  family = poisson, data = volcanoisland)
## )
## summary(islandPsem)
## 
## 
## 
## library(rjags)
## library(R2jags)
## 
## # 1) Model definition exactly how we created our data
## modelCode = "
## model{
## 
##   # Likelihood
##   for(i in 1:i.max){
## 
##     Wind[i] ~ dnorm(windPred[i], windPrec)
##     windPred[i] <- intW + altW*Alt[i]
## 
##     Ducks[i] ~ dbern(lambdaD[i])
##     logit(lambdaD[i]) <- intD + altD*Alt[i] + habitatD * Habitat[i] + soilD * SoilTexture[i]
## 
##     Beetles[i] ~ dpois(lambda[i])
##     lambda[i] <- exp(mu[i] )
##     mu[i] <- intB + altB*Alt[i] + alt2B*Alt[i]*Alt[i] + windB * Wind[i] + OLREB[i]
##   }
## 
##   # Prior distributions
##   intB ~ dnorm(0,0.001)
##   altB ~ dnorm(0,0.001)
##   alt2B ~ dnorm(0,0.001)
##   windB ~ dnorm(0,0.001)
##   for(i in 1:i.max){
##     OLREB[i] ~ dnorm(0,precOLREB)
##   }
##   precOLREB ~ dgamma(0.001,0.001)
## 
##   intW ~ dnorm(0,0.001)
##   altW ~ dnorm(0,0.001)
##   windPrec <- 1/(windSD * windSD)
##   windSD ~ dunif(0,100)
## 
##   intD ~ dnorm(0,0.001)
##   altD ~ dnorm(0,0.001)
##   habitatD ~ dnorm(0,0.001)
##   soilD ~ dnorm(0,0.001)
## 
##   # posterior predictive simulations
##   for(i in 1:i.max){
##     yPred[i] ~ dpois(lambda[i])
##   }
## }
## "
## 
## windPartiallObs <- islandData$windObs
## sel = sample.int(1000, 500)
## windPartiallObs[sel] = NA
## 
## # 2) Set up a list that contains all the necessary data (here, including parameters of the prior distribution)
## Data = list(Beetles = islandData$beetles, Alt = islandData$sAltitude, i.max = length(islandData$sAltitude), Wind = windPartiallObs, Ducks = islandData$ducks, Habitat = islandData$habitatQuality, SoilTexture = islandData$earth)
## 
## # 3) Specify a function to generate inital values for the parameters
## 
## # Out of laziness, we don't provide inits for the other parameters. For a real study, provide overdispersed sampling functions for all parameters
## inits.fn <- function() list(intB = rnorm(1), altB = rnorm(1), alt2B = rnorm(1))
## 
## 
## library(R2jags)
## 
## R2JagsResults <- jags(data=Data, inits=inits.fn, parameters.to.save=c("intB","altB","alt2B", "intW", "altW", "windB", "intD", "altD", "habitatD", "soilD"), n.chains=3, n.iter=10000, model.file=textConnection(modelCode), DIC = F)
## 
## plot(R2JagsResults)
## print(R2JagsResults)
## 
## 
## 
## library(DHARMa)
## simulations = R2JagsResults$BUGSoutput$sims.list$yPred
## pred = apply(simulations, 2, median)
## dim(simulations)
## sim = createDHARMa(simulatedResponse = t(simulations), observedResponse = islandData$beetles, fittedPredictedResponse = pred, integerResponse = T)
## plot(sim)
## 
## 
## plotResiduals(islandData$year, sim$scaledResiduals, asFactor = T)
## 
## testSpatialAutocorrelation(sim)
## 

