# playing with boral
library(boral)
## Taken straight from help file

library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
                       row.eff = "fixed", n.burnin = 10, n.iteration = 100, 
                       n.thin = 1, calc.ics = FALSE, save.model = TRUE)

summary(spider.fit.nb)

plot(spider.fit.nb, ask = FALSE, mfrow = c(2,2)) ## Plots used in residual analysis, 
## Used to check if assumptions such an mean-variance relationship 
## are adequately satisfied.

lvsplot(spider.fit.nb) ## Biplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.

rescors <- get.residual.cor(spider.fit.nb)

library(corrplot)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE)
## Example 2 - model with no latent variables, no site effects, 
## 	and environmental covariates
X <- scale(spider$x)
spider.fit.nb <- boral(y, X = X, family = "negative.binomial", 
                       num.lv = 0, n.burnin = 10, n.iteration = 100, n.thin = 1)

summary(spider.fit.nb) 
## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D

# add latent variables
spider.fit.nb3 <- boral(y, X = X, family = "negative.binomial", 
                       num.lv = 2, 
                       #n.burnin = 10, n.iteration = 100, n.thin = 1, 
                       save.model = TRUE)

summary(spider.fit.nb3) 
envcors <- get.enviro.cor(spider.fit.nb3)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE)
## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D

## Example 3 - Extend example 2 to demonstrate grouped and individual
## covariate selection
spider.fit.nb2 <- boral(y, X = X, family = "negative.binomial", 
                        num.lv = 0, n.burnin = 10, n.iteration = 100, n.thin = 1,
                        calc.ics = FALSE, ssvs.index = c(-1,-1,-1,0,1,2))

summary(spider.fit.nb2) 

## Example 3 - model fitted to presence-absence data, no site effects, and
## two latent variables
data(tikus)
y <- tikus$abun
y[y > 0] <- 1
y <- y[1:20,] ## Consider only years 1981 and 1983
y <- y[,apply(y > 0,2,sum) > 2] ## Consider only spp with more than 2 presences

tikus.fit <- boral(y, family = "binomial", num.lv = 2, 
                   n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE)

lvsplot(tikus.fit, biplot = FALSE) 
## A strong location between the two sampling years 


## Example 4 - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X
## will be estimated separately and not regressed against traits.
which.traits[[3]] <- 0

fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, 
                    family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, 
                    n.thin = 1, calc.ics = FALSE)

summary(fit.traits)

plot(fit.traits)
lvsplot(fit.traits)
