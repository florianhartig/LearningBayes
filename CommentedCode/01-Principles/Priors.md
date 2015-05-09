# Priors in Bayesian analysis
FlorianHartig  
9 May 2015  






The choice of prior (prior eliciation) is key to Bayesian analysis, and it is arguably the most contentious step in the whole procedure. 

The neccessity to specify priors was the main reason why Fisher and others developed MLE and NHST. 

If we speak about priors, we have to first make some disticitons 

* Uninformative priors are priors that are meant to express no prior information on a particular parameter
* Informative priors are the opposite 

## Scaling and scale-invariance of priors 





## Default choices for uniformative priors 

1. For **scale parameters** (something that affects the output linearly, like slope or intercept in a regression), use flat or quasi flat priors such as
  * A bounded uniform distribution
  * A wide normal distribution
  
2. It is possible to put a bit tighter priors around scale parameters to get the Bayesian analogue of Lasso or Ridge regression, see 
  * Park, T. & Casella, G. (2008) 
  * Kyung, M.; Gill, J.; Ghosh, M.; Casella, G. et al. (2010) Penalized regression, standard errors, and Bayesian lassos. Bayesian Analysis, 5, 369-411.
  * http://stats.stackexchange.com/questions/95395/ridge-regression-bayesian-interpretation?rq=1


2. For **variance parameters** (something like the standard deviation in a linear regression), use decaying parameters such as
  * 1/x (standard choice according to Jeffrey's prior)
  * inverse-gamma
  
3. For **variance hyperparameters in hierarchical models**, use
  * inverse-gamma
  * half-t family (suggested by Gelman, 2006)
  


## Readings

### Uninformative priors 

Kass, R. E. & Wasserman, L. (1996) The selection of prior distributions by formal rules. J. Am. Stat. Assoc., American Statistical Association, 91, 1343-1370.

Jeffreys, H. (1946) An Invariant Form for the Prior Probability in Estimation Problems. Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, The Royal Society, 186, 453-461.

Jaynes, E. (1968) Prior probabilities. Systems Science and Cybernetics, IEEE Transactions on, IEEE, 4, 227-241.

Tibshirani, R. (1989) Noninformative priors for one parameter of many. Biometrika, 76, 604-608.

Park, T. & Casella, G. (2008) The Bayesian Lasso. Journal of the American Statistical Association, 103, 681-686.

Irony, T. Z. & Singpurwalla, N. D. (1997) Non-informative priors do not exist -- a dialogue with José M. Bernardo. J. Stat. Plan. Infer., 65, 159-177.

Gelman, A.; Jakulin, A.; Pittau, M. G. & Su, Y.-S. (2008) A weakly informative default prior distribution for logistic and other regression models. The Annals of Applied Statistics, JSTOR, , 1360-1383.

Gelman, A. (2006) Prior distributions for variance parameters in hierarchical models. Bayesian Analysis, Citeseer, 1, 515-533.

Fong, Y.; Rue, H. & Wakefield, J. (2010) Bayesian inference for generalized linear mixed models. Biostatistics, 11, 397-412.

Ferguson, T. (1974) Prior distributions on spaces of probability measures. The Annals of Statistics, JSTOR, 2, 615-629.


### Informative priors 

Choy, S. L.; O'Leary, R. & Mengersen, K. (2009) Elicitation by design in ecology: using expert opinion to inform priors for Bayesian statistical models. Ecology, 90, 265-277




---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License

