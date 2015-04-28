Fitting complex or process-based models
===

In principle, there is nothing different when comparing 
 

## How to do this 

### Preparing compiled models to 

See [LinkingAModelToR.md](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/09-BayesAndProcessBasedModels/LinkingAModelToR.md)

### Defining fit

How should we measure how well the model fits to the data? There are different options 

* Ad hoc objectives
* Likelihood
* Posterior = Likelihood + Prior

For an overview of this discussion, see 

* Hartig et al. (2012): Connecting dynamic vegetation models to data - an inverse perspective. <i>J. Biogeogr.</i>, 39, 2240-2252. <a href="http://dx.doi.org/10.1111/j.1365-2699.2012.02745.x">[journal]</a> <a href="http://ecologypapers.blog.com/files/2012/10/Hartig-Connectingdynamicvegetation-2012.pdf">[pdf]</a>


### Getting best fit and uncertainties 

In general, all the approaches using optimization, MCMC and sensitivity analysis in this collection of tuturials are valid for forest models as well, and it's worth looking at them. In particular, also ABC approaches may be interesting to look at. Some references particular to forest models are:

* Tutorial by Marcel van Oijen on fitting forest models with Bayes  <a href="http://nora.nerc.ac.uk/6087/1/BC%26BMC_Guidance_2008-12-18_Final.pdf">here</a>.
* Hartig et al. (2012): Connecting dynamic vegetation models to data - an inverse perspective. <i>J. Biogeogr.</i>, 39, 2240-2252. <a href="http://dx.doi.org/10.1111/j.1365-2699.2012.02745.x">[journal]</a> <a href="http://ecologypapers.blog.com/files/2012/10/Hartig-Connectingdynamicvegetation-2012.pdf">[pdf]</a>
* Van Oijen, M.; Rougier, J. &amp; Smith, R. (2005) Bayesian calibration of process-based forest models: bridging the gap between models and data Tree Physiol., 25, 915-927
* Hartig, F.; Dislich, C.; Wiegand, T. & Huth, A. (2014) Technical Note: Approximate Bayesian parameterization of a process-based tropical forest model. Biogeosciences, 11, 1261-1272.
* Talk in Zürich by Florian Hartig <a href="http://florianhartig.files.wordpress.com/2013/12/inversemodelling.pdf">pdf</a>


