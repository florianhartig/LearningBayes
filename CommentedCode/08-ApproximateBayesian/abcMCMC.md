A simple Approximate Bayesian Computation MCMC (ABC-MCMC) in R
====





Approximate Bayesian Computing and similar techniques, which are based on calculating approximate likelihood values based on samples from a stochastic simulation model, have attracted a lot of attention in the last years, owing to their promise to provide a general statistical technique for stochastic processes of any complexity, without the limitations that apply to “traditional” statistical models due to the problem of maintaining “tractable” likelihood functions. 

If you want to have more background on this algorithm, read the excellent paper by Marjoram et al. (2003) who proposed this algorithm for the first time. 


```r
library(coda)
```

```
## Warning: package 'coda' was built under R version 3.1.3
```

```r
# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
data =  rnorm(10, mean =5.3, sd = 2.7)
 
 
# we want to use ABC to infer the parameters that were used. 
# we sample from the same model and use mean and variance
# as summary statstitics. We return true for ABC acceptance when
# the difference to the data is smaller than a certain threshold
 
meandata <- mean(data)
standarddeviationdata <- sd(data)
 
ABC_acceptance <- function(par){
   
  # prior to avoid negative standard deviation
  if (par[2] <= 0) return(F) 
   
  # stochastic model generates a sample for given par
  samples <- rnorm(10, mean =par[1], sd = par[2])
 
  # comparison with the observed summary statistics
  diffmean <- abs(mean(samples) - meandata)
  diffsd <- abs(sd(samples) - standarddeviationdata)
  if((diffmean < 0.2) & (diffsd < 0.3)) return(T) else return(F)
}
 
 
# we plug this in in a standard metropolis hastings MCMC, 
# with the metropolis acceptance exchanged for the ABC acceptance
 
run_MCMC_ABC <- function(startvalue, iterations){
 
    chain = array(dim = c(iterations+1,2))
    chain[1,] = startvalue
 
    for (i in 1:iterations){
         
        # proposalfunction
        proposal = rnorm(2,mean = chain[i,], sd= c(0.7,0.7))
         
        if(ABC_acceptance(proposal)){
            chain[i+1,] = proposal
        }else{
            chain[i+1,] = chain[i,]
        }
    }
    return(mcmc(chain))
}
 
posterior <- run_MCMC_ABC(c(4,2.3),300000)
```

The result should look something like that:

```r
plot(posterior)
```

![](abcMCMC_files/figure-html/plot_posterior-1.png) 


## References for further reading

Marjoram, P.; Molitor, J.; Plagnol, V. & Tavaré, S. (2003) Markov chain Monte Carlo without likelihoods


---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
