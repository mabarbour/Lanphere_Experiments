
## simulate a simple linear model
n <- 200
beta <- 0.3
x <- rep(c(0,1), n/2)
Intercept <- 1
y <- Intercept + beta*x + rnorm(n, 0, 1)

## analyse the unscaled and scaled versions
lm.y <- lm(y ~ x)
coef(lm.y)

lm.scaled <- lm(scale(y) ~ scale(x))
coef(lm.scaled)

## for the unscaled model, there is no correspondence between the coefficient and standard deviation
round(coef(lm.y),3)[2]
round(sd(as.vector(coef(lm.y)[2] * t(model.matrix(~x))[2, ])),3)

## for the scaled model though, they are equivalent. perhaps this doesn't hold though when there is imbalance in the model
round(coef(lm.scaled),3)[2]
round(sd(as.vector(coef(lm.scaled)[2] * t(model.matrix(~scale(x)))[2, ])),3)

## Based on the rationale laid out in Nakagawa and Shielzeth (ME&E 2012: http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x/full),
## the logical extension is that I could use the SD for both fixed and random effects (for scaled responses and predictors)
## as standardized effect sizes. This, in combination with the common practice in structural equation models of multiplying
## standardized coefficients to obtain indirect effects, should enable me to calculate the strength of indirect effects
## for both fixed and random effects on the same scale.
