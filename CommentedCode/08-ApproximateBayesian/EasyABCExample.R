library(EasyABC)

# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
data =  rnorm(10, mean =5.3, sd = 2.7)

# we want to use ABC to infer the parameters that were used.
# we sample from the same model and use mean and variance
# as summary statstitics for the model and the data.

# observed summary statistics
summarydata = c(mean(data), sd(data))

# stochastic model generates a sample for given par and returns summary statistics
model <- function(par){ 
  samples <- rnorm(10, mean =par[1], sd = par[2])
  return(c(mean(samples), sd(samples)))
}

# call to EasyABC with the ABC-MCMC algorithm Marjoram, P.; Molitor, J.; Plagnol, V. & 
# Tavare, S. (2003) Markov chain Monte Carlo without likelihoods. Proc. Natl. Acad. Sci. USA, 100, 15324-15328.
# with some automatic adjustment options 
ABC_Marjoram_original<-ABC_mcmc(method="Marjoram", model=model, 
  prior=list(c("unif",0,10),c("unif",1,5)), 
  summary_stat_target=summarydata, n_rec = 10000)


str(ABC_Marjoram_original)
par(mfrow=c(2,1))
hist(ABC_Marjoram_original$param[5000:10000,1], main = "Posterior for mean")
hist(ABC_Marjoram_original$param[5000:10000,2], main = "Posterior for standard deviation")
