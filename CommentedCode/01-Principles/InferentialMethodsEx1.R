rm(list=ls(all=TRUE))

#Analysis of likelihood, p-value, and Bayes for binomial model,
#10 trials, 3 success, unknown coin, want to do inference

trials = 10
success = 8

# GETTING THE MLE ESTIMATE

# task: use dbinom to plot the probability of the observed success as a function of the model parameter ... find the MLE

testvalues <- seq(0,1,0.01)
likelihood <- dbinom(success, trials, testvalues)
plot(testvalues, likelihood, type = "l")
MLE <- testvalues[which.max(likelihood)]
abline(v=MLE, col="red")

# GETTING A P-VALUE FOR FAIR COIN
# want to get p-value for a smaller or equal result (1-tailed) given a fair coin p(k<=kobs|H0:p=0.5)

# TASK: Use dbinom or pbinom to caculate the probability of having less or equal sucess than observed, given the assumption (H0) that the true parameter value is 0.5

H0 = 0.5 

pvalue <- sum(dbinom(8:10, 10, prob = H0) )
pvalue
pvalue < 0.05

binom.test(8, 10, p=0.5, alternative = "greater")

# take care to choose the right tail (lower / upper)


# ALPHA DEMONSTRATION

# Task: assume a null hypothesis, and draw 1000 x data from the null hypothesis, calculate the p-value. How often is the p-value < 0.05, although the null hypothesis is true (Type I error)


trials <- 1000
probsuccess <- 0.5
repetitions <- 10000

outcomes <- rbinom(repetitions, trials, probsuccess)
hist(outcomes)
pvalues <- pbinom(outcomes, trials, 0.5)
hist(pvalues)
sum(pvalues < 0.05) / repetitions


# BAYESIAN POSTERIOR

# plot the Bayesian posterior, define as likelihood * prior for

# a) a flat prior
# b) a prior centered around 0.5, e.g. with a normal distribution

plot(testvalues, likelihood, type = "l")
MLE <- testvalues[which.max(likelihood)]
abline(v=MLE, col="red")

flatprior <- 1
posteriorflat <- likelihood * flatprior / sum(likelihood * flatprior) * 101

plot(testvalues, posteriorflat, col = "green")
lines(testvalues, likelihood)

centeredprior <- dnorm(testvalues, 0.5, 0.1) 
posteriorcentered <- likelihood * centeredprior / sum(likelihood * centeredprior) * 101

plot(testvalues, posteriorcentered, col = "red")
lines(testvalues, likelihood)
lines(testvalues, centeredprior, col= "blue")





