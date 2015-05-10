Scripts to demonstrate how an MCMC works 

* [Intro to MH-MCMC](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/Metropolis.md)
* [Convergence diagnostics](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/Convergence.md)
* [Different MCMC algorithms with Laplaces Demon](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/LaplacesDeamon.md)
* [The DREAM algorithm](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/Dream.md)



### Software for MCMC sampling

* For sampling statistical models, most people use Bugs, Jags or Stan that include not only a sampler, but also a language to specify the likelihood and priors. See explanations [here](https://github.com/florianhartig/LearningBayes/tree/master/CommentedCode/03-Software)

* If you want a sampler that is able to sample from a general target distribution, a large range of samplers is in the package Laplaces Demon (see [script](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/LaplacesDeamon.md)). If you want to compare samplers or learn about different ones, I would start there. 

* Packages that implement one or a small number of samplers include 

  * [adaptMCMC](http://cran.r-project.org/web/packages/adaptMCMC/)  (adaptive MH)


### References for further reading 

http://arxiv.org/abs/1504.01896




---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
