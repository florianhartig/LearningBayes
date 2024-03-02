a <- 5
b <- 10
sigma <- 10
rsigma = 30
group = rep(1:11, each = 5)
randomEffect = rnorm(11, sd = rsigma)

x <- -27:27
y <- a * x + b + rnorm(55,0,sd = sigma) + randomEffect[group]
plot(x,y, col = group, pch = 3)
dat = data.frame(x, y, group)

library(effects)
# neglect group
fit <- lm(y ~ x, data = dat)
summary(fit)
ranef(fit)
plot(allEffects(fit, partial.residuals = T))

# adjust for group via fixed effect
groupF = as.factor(group)
fit <- lm(y ~ x + groupF)
summary(fit)
plot(allEffects(fit, partial.residuals = T ))

# adjust for group via RE
library(lme4)
fit <- lmer(y ~ x + (1|group), data = dat)
summary(fit)
plot(predict(fit), residuals(fit))
