# Information for participants of the course "Introduction to Bayesian Statistics" Leipzig '15
Florian Hartig  
3 Jun 2015  

### Lecturers

* [Florian Hartig](http://florianhartig.wordpress.com/)
* [Felix May](https://www.ufz.de/index.php?en=32342)

### Course plan

The course runs from 9.00 to 17.00, with a one-hour lunch break around 12.30

1. Day 1 - morning (FH)
  * Intro Bayes, and difference to "conventional" stats 
  * Priors
  * Posterior interpretation 
  * MCMC sampling 

2. Day 1 - Afternoon (FM)
  * Intro Jags
  * Binomial model in JAGS
  * Regression in JAGS

3. Day 2 - Morning (FH)
  * Mixed and generalised linear mixed models in JAGS 
  * Model checking, Bayesian p-values 

4. Day 2 - Afternoon (FH)
  * Hierarchical models
  * Spatial models 
  * Outlook INLA 
  * Bayesian Model selection / Model averaging 
  * Approximate Bayesian Computation (ABC)


### Preparations and installation of software

* If you want to prepare for the course in advance (we do not expect that you will, but just in case), have a look at the recommended reading material [here](http://florianhartig.github.io/LearningBayes/), in particular the nice and free [Bayes intro by Michael Clark](http://www3.nd.edu/%7Emclark19/learn/IntroBayes.pdf)

* Please make sure the following software is installed on your laptops
  * [R](http://www.r-project.org/)
  * [Rstudio](http://www.rstudio.com/)
  * [JAGS](http://mcmc-jags.sourceforge.net/)
  * The R packages [rjags](http://cran.r-project.org/web/packages/rjags/index.html) and [R2jags](http://cran.r-project.org/web/packages/R2jags/index.html)

* Check that everything runs fine by running [the following code](https://raw.githubusercontent.com/florianhartig/LearningBayes/master/CommentedCode/03-Software/Jags/testR2jags.R)


```r
# Test of the R2jags system
# Modified from the help file of the jags function


# An example model file is given in:
model.file <- system.file(package="R2jags", "model", "schools.txt")
# Let's take a look:
file.show(model.file)
# you can also write BUGS model as a R function, see below:

#=================#
# initialization  #
#=================#

# data
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)


jags.data <- list("y","sd","J")
jags.params <- c("mu","sigma","theta")
jags.inits <- function(){
  list("mu"=rnorm(1),"sigma"=runif(1),"theta"=rnorm(J))
}

## You can input data in 4 ways
## 1) data as list of character
jagsfit <- jags(data=list("y","sd","J"), inits=jags.inits, jags.params,
                n.iter=10, model.file=model.file)

## 2) data as character vector of names
jagsfit <- jags(data=c("y","sd","J"), inits=jags.inits, jags.params,
                n.iter=10, model.file=model.file)

## 3) data as named list
jagsfit <- jags(data=list(y=y,sd=sd,J=J), inits=jags.inits, jags.params,
                n.iter=10, model.file=model.file)

## 4) data as a file
fn <- "tmpbugsdata.txt"
dump(c("y","sd","J"), file=fn)
jagsfit <- jags(data=fn, inits=jags.inits, jags.params,
                n.iter=10, model.file=model.file)
unlink("tmpbugsdata.txt")

## You can write bugs model in R as a function

schoolsmodel <- function() {
  for (j in 1:J){                     # J=8, the number of schools
    y[j] ~ dnorm (theta[j], tau.y[j]) # data model:  the likelihood
    tau.y[j] <- pow(sd[j], -2)        # tau = 1/sigma^2
  }
  for (j in 1:J){
    theta[j] ~ dnorm (mu, tau)        # hierarchical model for theta
  }
  tau <- pow(sigma, -2)               # tau = 1/sigma^2
  mu ~ dnorm (0.0, 1.0E-6)            # noninformative prior on mu
  sigma ~ dunif (0, 1000)             # noninformative prior on sigma
}

jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params,
                n.iter=10, model.file=schoolsmodel)


#===============================#
# RUN jags and postprocessing   #
#===============================#
jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params,
                n.iter=5000, model.file=model.file)

# display the output
print(jagsfit)
plot(jagsfit)

# or to use some plots in coda
# use as.mcmmc to convert rjags object into mcmc.list
plot(as.mcmc(jagsfit))
```

