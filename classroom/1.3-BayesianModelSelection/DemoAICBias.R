# This example shows how AIC selection, followed by a conventional regression analysis of the selected model, massively inflates false positives. CC BY-NC-SA 4.0 Florian Hartig 

set.seed(1)
library(MASS)

dat = data.frame(matrix(runif(20000), ncol = 100))
dat$y = rnorm(200)
fullModel = lm(y ~ . , data = dat)
summary(fullModel)
# 2 predictors out of 100 significant (on average, we expect 5 of 100 to be significant)

selection = stepAIC(fullModel)
summary(lm(y ~ X1 + X2 + X3 + X5 + X7 + X13 + X20 + X23 + X30 + 
             X37 + X42 + X45 + X46 + X47 + X48 + X64 + X65 + X66 + X71 + 
             X75 + X80 + X81 + X87 + X88 + X89 + X90 + X94 + X100, data = dat))

# voila, 15 out of 28 (before 100) predictors significant - looks like we could have good fun to discuss / publish these results!

dat2 = dat
for(i in 1:100) dat2[,i] = scale(dat2[,i])

Data = list(y = dat2$y, x = as.matrix(dat2)[,1:100], i.max = nrow(dat2))

modelCode = "model{

# Likelihood
for(i in 1:i.max){
  mu[i] <- inprod(a , x[i,]) + b
  y[i] ~ dnorm(mu[i],tau)
}

# Prior distributions

for(i in 1:100){
  a[i] = dnorm(0,tauShrinkage)
}
b ~ dnorm(0,tauShrinkage)

tauShrinkage ~ dgamma(0.001, 0.001)
sdShrinkage <- 1/sqrt(tauShrinkage)

# For scale parameters, normal choice is decaying
tau ~ dgamma(0.001, 0.001)
sigma <- 1/sqrt(tau)

}
"

jagsModel <- jags.model(file= textConnection(modelCode), data=Data, n.chains = 3)

para.names <- c("a","b","sigma", "sdShrinkage")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)
summary(Samples)