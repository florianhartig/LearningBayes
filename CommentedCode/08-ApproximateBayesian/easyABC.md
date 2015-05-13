The EasyABC package for Approximate Bayesian Computation in R
====






Approximate Bayesian Computation (ABC) is a relatively new method that allows treating any stochastic model (IBM, stochastic population model, …) in a statistical framework by generating “approximate” likelihood values through simulating from the model. We provide a gentle introduction to ABC and some alternative approaches in our recent Ecology Letters review on “statisitical inference for stochastic simulation models”.

ABC has a huge potential as a solution for many typical ecological problems, but to make this more widely known is currently hindered by the fact that you have to code everything by hand, which excludes a large number of users.

The EasyABC package, available from CRAN (developed by Franck Jabot, Thierry Faure, Nicolas Dumoulin and maintained by Nicolas Dumoulin.), implements a number of algorithms for the three main sampling strategies used in ABC, namely Rejection Sampling, Sequential Monte Carlo (SMC) and Markov Chain Monte Carlo (MCMC). All those are also discussed in our review. The use of the package is relatively straightforward. 



```r
library(EasyABC)
```

```
## Warning: package 'EasyABC' was built under R version 3.1.3
```

```
## Loading required package: abc
```

```
## Warning: package 'abc' was built under R version 3.1.3
```

```
## Loading required package: abc.data
```

```
## Warning: package 'abc.data' was built under R version 3.1.3
```

```
## Loading required package: nnet
## Loading required package: quantreg
## Loading required package: SparseM
## 
## Attaching package: 'SparseM'
## 
## The following object is masked from 'package:base':
## 
##     backsolve
## 
## Loading required package: MASS
## Loading required package: locfit
## locfit 1.5-9.1 	 2013-03-22
```

```r
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
```

```
## List of 8
##  $ param              : num [1:10000, 1:2] 7.49 7.23 7.23 7.23 7.23 ...
##  $ stats              : num [1:10000, 1:2] 6.23 6.03 6.03 6.03 6.03 ...
##  $ dist               : num [1:10000(1d)] 0.006045 0.000747 0.000747 0.000747 0.000747 ...
##  $ stats_normalization: num [1:2] 3.08 1.34
##  $ epsilon            : num 0.0338
##  $ nsim               : num 109991
##  $ n_between_sampling : num 10
##  $ computime          : num 14.5
```

```r
par(mfrow=c(2,1))
hist(ABC_Marjoram_original$param[5000:10000,1], main = "Posterior for mean")
hist(ABC_Marjoram_original$param[5000:10000,2], main = "Posterior for standard deviation")
```

![](easyABC_files/figure-html/main-1.png) 
<br />
<br />

## References for further reading

Hartig, F., Calabrese, J. M., Reineking, B., Wiegand, T. and Huth, A. (2011), Statistical inference for stochastic simulation models – theory and application. Ecology Letters, 14: 816–827. 


---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
