---
layout: page
title: Bayesian statistics 
category: hubs
---

# Introduction


* Hobbs, N. T. & Hilborn, R. (2006) Alternatives to statistical hypothesis testing in ecology: A guide to self teaching Ecol. Appl., 16, 5-19
* Ellison, A. M. (2004) Bayesian inference in ecology Ecol. Lett., 7, 509-520


## Prior choice

* Kass, R. E. & Wasserman, L. (1996) The selection of prior distributions by formal rules. J. Am. Stat. Assoc., 91, 1343-1370


### Foundations of Bayesian statistics, Bayes vs. Frequentists

* Efron, B. (2013) A 250-year argument: Belief, behavior, and the bootstrap Bulletin Of The American Mathematical Society, 50, 129-146
* Gelman, A. & Robert, C. P. (2010) ”Not only defended but also applied”: The perceived absurdity of Bayesian inference ArXiv e-prints
* Fisher, R. A. (1922) On the mathematical foundations of theoretical statistics Philos. T. Roy. Soc. A., 222, 309-368
* Kass, R. (2011) Statistical inference: The big picture Stat. Sci., 26, 1-9
* Jaynes, E. (1976) Confidence intervals vs. Bayesian intervals Foundations of probability theory, statistical inference, and statistical theories of science, 2, 175-257.


### Introductory textbooks

A necessarily subjective and incomplete list of textbooks and papers that are interesting in the context of learning Bayesian statistics, compiled with input from Joe Chipperfield and Jörn Pagel for our summer schools on Bayesian statistics.
Textbooks

#### Basic

* Kéry, M. (2010) Introduction to WinBUGS for Ecologists. Academic Press.
* Kruschke, J. F. (2010) Doing Bayesian Data Analysis: A Tutorial with R and BUGS. Academic Press.

#### Comprehensive

* Lunn D. et al. (2012) The BUGS Book: A Practical Introduction to Bayesian Analysis. Chapman and Hall/CRC.
* Gelman, A.; Carlin, J. B.; Stern, H. S. & Rubin, D. B. (2003) Bayesian Data Analysis. Chapman & Hall, London.

#### Hierarchical

* Kéry, M. and Schaub, M. (2011) Bayesian population analysis using WinBUGS. Academic Press.
* Banerjee, S. et al. (2009) Hierarchical Modeling and Anallysis for Spatial Data. Chapman and Hall/CRC.
* Clark, J. S. and Gelfand, A. E. (2006) Hierarchical Modelling for the Environmental Sciences. Oxford University Press.

# Working with Bayesian Methods

## MCMC sampling



* Andrieu, C.; de Freitas, N.; Doucet, A. & Jordan, M. I. (2003) An introduction to MCMC for machine learning Mach. Learning, 50, 5-43
* Andrieu, C. & Thoms, J. (2008) A tutorial on adaptive MCMC Stat. Comput., 18, 343-373.

### Summaries of the posterior 

http://www.sumsar.net/blog/2014/10/probable-points-and-credible-intervals-part-one/



# Hierarchical Models

## Mixed models

## Observer models 

## State-space models

## Literature

* Wikle, C. K. (2003) Hierarchical Bayesian models for predicting the spread of ecological processes Ecology, 84, 1382-1394
* Clark, J. S. (2003) Uncertainty and variability in demography and population growth: A hierarchical approach Ecology, 84, 1370-1381
* Clark, J. S. & Gelfand, A. E. (2006) A future for models and data in environmental science. Trends in Ecology & Evolution, 21, 375-380
* Cressie, N.; Calder, C. A.; Clark, J. S.; Hoef, J. M. V. & Wikle, C. K. (2009) Accounting for uncertainty in ecological analysis: the strengths and limitations of hierarchical statistical modeling Ecol. Appl., 19, 553-570
* Marion, G.; McInerny, G. J.; Pagel, J.; Catterall, S.; Cook, A. R.; Hartig, F. & O’Hara, R. B. (2012) Parameter and uncertainty estimation for process-oriented population and distribution models: data, statistics and the niche J. Biogeogr., 39, 2225–2239
* Cook, A.; Marion, G.; Butler, A. & Gibson, G. (2007) Bayesian Inference for the Spatio-Temporal Invasion of Alien Species Bull. Math. Biol., 69, 2005-2025
* Pagel, J. & Schurr, F. M. (2011) Forecasting species ranges by statistical estimation of ecological niches and spatial population dynamics Global Ecol. Biogeogr.

# Bayesian Model selection

* DIC
* Bayes Factor
* Reversible Jump MCMC (RJ-MCMC)

While DIC is really a different idea, it should be noted that BF and RJ-MCMC have the same underlying logic, i.e. they are based on the full posterior P(D|M, \phi) p(M) p(\phi) , just that for the BF the parameters are marginalized out. 

## DIC

* Note that there are potentially differences in DIC calculations in bugs and jags. 

## Bayes Factors

* Kass, R. E. & Raftery, A. E. (1995) Bayes Factors J. Am. Stat. Assoc., 90, 773-795

https://radfordneal.wordpress.com/2008/08/17/the-harmonic-mean-of-the-likelihood-worst-monte-carlo-method-ever/

http://www.biomedcentral.com/1471-2105/14/85

http://hedibert.org/bayes-factor-computing-marginal-likelihoods-savage-dickey-ratio-reversible-jump-mcmc-bayesian-model-averaging-and-deviance-information-criterion/


## RJ-MCMC



# Fitting (stochastic) process-based models

See our extra hub on this topic

* Van Oijen, M.; Rougier, J. & Smith, R. (2005) Bayesian calibration of process-based forest models: bridging the gap between models and data Tree Physiol., 25, 915-927
* Beaumont, M. A. (2010) Approximate Bayesian computation in evolution and ecology Annu. Rev. Ecol. Evol. Syst., 41, 379-406
* Csilléry, K.; Blum, M. G. B.; Gaggiotti, O. E. & François, O. (2010) Approximate Bayesian Computation (ABC) in practice Trends in Ecology & Evolution, 25, 410-418
* Hartig, F.; Calabrese, J. M.; Reineking, B.; Wiegand, T. & Huth, A. (2011) Statistical inference for stochastic simulation models – theory and application Ecol. Lett., 14, 816-827
* Jabot, F. & Chave, J. (2009) Inferring the parameters of the neutral theory of biodiversity using phylogenetic information and implications for tropical forests Ecol. Lett., 12, 239-248
* Hartig, F.; Dyke, J.; Hickler, T.; Higgins, S. I.; O’Hara, R. B.; Scheiter, S. & Huth, A. (2012) Connecting dynamic vegetation models to data – an inverse perspective J. Biogeogr., 39, 2240-2252.
