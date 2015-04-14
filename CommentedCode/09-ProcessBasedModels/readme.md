---
layout: page
title: Fitting process-based models
category: hubs
author: Florian Hartig
synopsis: A tutorial on how to fit complex or process-based (vegetation) models to data
---

Fitting complex or process-based models
===

# Overview

A common problem in ecological modelling is that a simulation model (vegetation model, forest simulator, individual-based model, movement model) has been designed, and how one would like to know

* How well does the model fit to data
* Which are the parameters that deliver the best fit
* How sensitive is the model when parameters are changed
* What is the uncertainty in parameters and predictions

There are a number of methods and different steps to achieve this. Some of the questions are

* How do I define "fit", i.e. how do we measure the difference between models and data
* How can I vary model parameters systematically to find
 * Sensitivity (how much outputs or fit change with the parameters)
 * Parameters with the best fit
* How do I quantify uncertainty on
 * Model structure (connected to validation)
 * Model parameters
 * Model predictions

## Preparing your model to work on these questions

Before we can start working on these questions at all, we need to make sure that model parameters can be systematically varied from some programming environment. So

* Make your model callable from R, Python or whatever environment you prefer, so that you can change the parameters from this environment, e.g. by creating a runModel(parameter) function
* Write a function to get the model output into R, Python, e.g. getData(dataype)

**Hint**: for many existing methods, the easiest way to do this will be via HD I/O, i.e. your run model function would change the ini file and then run the model, and your getData() function would read the model outputs from HD. However, this may create a lot of HD traffic, and can then cause problems if you run in parallel on larger cluster systems. I would recommend to start via HD, but then consider later whether you can interface directly.


# Defining fit

How should we measure how well the model fits to the data?

I will add more info later, for the momemnt read our discussion in Hartig et al. (2012): Connecting dynamic vegetation models to data - an inverse perspective. <i>J. Biogeogr.</i>, 39, 2240-2252. <a href="http://dx.doi.org/10.1111/j.1365-2699.2012.02745.x">[journal]</a> <a href="http://ecologypapers.blog.com/files/2012/10/Hartig-Connectingdynamicvegetation-2012.pdf">[pdf]</a>


# Sensitivity analsis

to be added



#Bayesian statistics


##What to do to get your model fitted with Bayes

<ol>
	<li>Decide on the likelihood function that is used to compare your model results to data (see Hartig et al. 2012 below, and/or a stats textbook on likelihood, e.g. Bolker, B. (2008) Ecological models with R Princeton UP) and write this as a function in R</li>
	<li>Sample likelihood with MCMC, see blog posts below</li>

</ol>


*Blog posts with example code*

<ul>
	<li><a href="http://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/">An annotated MCMC example in R</a></li>
	<li><a href="http://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/">MCMC convergence analysis in R</a></li>
	<li><a href="http://theoreticalecology.wordpress.com/2012/07/15/a-simple-approximate-bayesian-computation-mcmc-abc-mcmc-in-r/">A simple ABC example in R</a></li>
</ul>



##References and further links

See also <a href="http://florianhartig.wordpress.com/teaching/bayes-intro/literature-suggestions-bayesian-statistics/">our general Bayesian literature list</a>


*Normal Bayes with forest models*

Hartig et al. (2012): Connecting dynamic vegetation models to data - an inverse perspective. <i>J. Biogeogr.</i>, 39, 2240-2252. <a href="http://dx.doi.org/10.1111/j.1365-2699.2012.02745.x">[journal]</a> <a href="http://ecologypapers.blog.com/files/2012/10/Hartig-Connectingdynamicvegetation-2012.pdf">[pdf]</a>

Tutorial by Marcel van Oijen on fitting forest models with Bayes  <a href="http://nora.nerc.ac.uk/6087/1/BC%26BMC_Guidance_2008-12-18_Final.pdf">here</a>.


<ul>

	<li>Van Oijen, M.; Rougier, J. &amp; Smith, R. (2005) Bayesian calibration of process-based forest models: bridging the gap between models and data Tree Physiol., 25, 915-927</li>


	<li>van Oijen et al. (2013) Bayesian calibration, comparison and averaging of six forest models, using data from Scots pine stands across Europe For. Ecol. Manage., 289, 255-268</li>

</ul>




*ABC and other simulation-based approximations*

<ul>

	<li>Jabot, F. &amp; Chave, J. (2009) Inferring the parameters of the neutral theory of biodiversity using phylogenetic information and implications for tropical forests Ecol. Lett., 12, 239-248</li>

	<li>Hartig et al. (2013) Technical Note: Approximate Bayesian parameterization of a complex tropical forest model. Biogeosciences. http://www.biogeosciences.net/11/1261/2014/bg-11-1261-2014.html </li>

	<li> Hartig et al. (2011): Statistical inference for stochastic simulation models - theory and application. <i>Ecol. Lett. </i> 14, 816-827. <a href="http://dx.doi.org/10.1111/j.1461-0248.2011.01640.x">[journal]</a> 
</ul>


<strong>General talks</strong>

<ul>
	<li><a href="http://de.slideshare.net/florianhartig/florian-hartig-gfoe13s">GFÖ 2013</a></li>
	<li><a href="http://florianhartig.files.wordpress.com/2013/12/inversemodelling.pdf">Zürich</a></li>
</ul>



