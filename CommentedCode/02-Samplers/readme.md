This section provides examples and explanations of numerical algorithms that vary parameters to infer the parameters that fit best, uncertainty around the best fit, and sensitivity of the model to parameter changes. 

The **aim of this section** is to understand how these algorithms work. In many cases, you will not need to apply them directly, because they are included and chosen automatically in more Bayesian such as JAGS or STAN that will be presented in the section 3. However, apart from understanding what's going on, there may be cases where you can't specify a model in JAGS or STAN, and hence the need to fall back to the more hands-on implementations that are provided here. 

The topics that are covered here are

* [Optimization](https://github.com/florianhartig/LearningBayes/tree/master/CommentedCode/02-Samplers/Optimization)
* Sampling algorithms
  * [Rejection Sampling](https://github.com/florianhartig/LearningBayes/tree/master/CommentedCode/02-Samplers/Rejection)
  * [Markov-Chain Monte Carlo](https://github.com/florianhartig/LearningBayes/tree/master/CommentedCode/02-Samplers/MCMC)
  * [Sequential Monte Carlo](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/SMC/SMC.md)
* Sensitivity analysis
* Benchmarking of optimization and sampling algorithms
