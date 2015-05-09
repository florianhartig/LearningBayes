Bayesian inference uses the same models than normal statistics. It simply is a different way to evaluate models. Instead of a significance tests, or a maximum likelihood estimate, Bayesian inference calculates a posterior distribution. The posterior is defined as 

P(phi|D) = c * p(D|phi) * p(phi) 

where phi are the parameters of the model, D is the data, and c is a normalization constant that normalizes the posterior to 1. The following scripts explain this approach in more detail.

* The script [inference methods](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/InferenceMethods.md) compares the approach of Bayesian analysis to the other two main approaches of statistical inference, NHST and MLE.
* The script [priors](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/Priors.md) speaks about choice of prior distributions
* The script [posterior](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/Posterior.md) explains how to interpret the posterior distribution


---
**Copyright, reuse and updates**: By Florian Hartig. Updates will be posted at https://github.com/florianhartig/LearningBayes. Reuse permitted under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
