# extended commented version of this example at
# https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/InferenceMethods.md

# Assume we flipped a coin 10 times, and want to find out if
# it is biased - what can we learn about the probability to
# obtain heads with this coin?

trials = 10
success = 7

# For all three statistical methods, we use the same statistical model
# which is the binomial model. The probability density is available in R
# throught he function dbinom - dbinom(6,10,0.9) gives you the probability
# density of obtaining 6/10 heads when the true probability of heads is 0.9
dbinom(6,10, 0.5)

# We will now use this model to calculate the three classical inferential
# outputps of statistics - tests, MLE and Bayes

########### NHST #####################

# assume the coin is random (0.5), p-value is p >= observed
barplot(dbinom(0:10,10, 0.5), names.arg = 0:10, col = c(rep("grey", 7), rep("red", 4)))

# have to do 6 because R's diosyncrasies with lower.tail
pbinom(6,10,0.5, lower.tail = F)

binom.test(7,10,alternative = "greater")

############ MLE ######################

likelihood = function(x) dbinom(7,10, x)
parameterValues = seq(0,1,length.out = 100)

# assume data is fixed, true probability unknown
plot(parameterValues, likelihood(parameterValues), type = "l")

############# Bayes ##################

# posterior

prior = function(x) dnorm(x, mean = 0, sd = 0.1)

par(mfrow = c(2,2))

plot(parameterValues, prior(parameterValues), type = "l", main = "Prior")

plot(parameterValues, likelihood(parameterValues), type = "l", main = "Likelihood")


plot(parameterValues, prior(parameterValues) * likelihood(parameterValues), type = "l", main = "Posterior")
