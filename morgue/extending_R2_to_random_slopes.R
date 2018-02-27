############################################################################
# 
# Supporting information for the article "Extension of Nakagawa & 
#   Schielzeth's R^2[GLMM] to random slopes models" 
# by Paul Johnson
#
# This R script illustrates estimation of the random effects variance
# component from random slopes GLMMs, an extension of the 
# GLMM R-squared method of Nakagawa & Schielzeth (2013)
# for random intercepts models. 
#
# Paul Johnson, Institute of BAHCM, University of Glasgow
# 31th March 2014
#
# Main reference:
#   Nakagawa, S. & Schielzeth, H. (2013). A general and simple method for 
#   obtaining R^2 from generalized linear mixed-effects models. 
#   Methods in Ecology and Evolution, 4, 133-142.
#
############################################################################

# Preamble

# The aim this script is to illustrate, using three examples, estimation of 
# marginal and conditional R-squared for random slopes models. Each example
# uses data where random slopes improve the fit of the model. R-squared values
# are also estimated from models fitted to the same data without random slopes
# (the workaround recommended by Nakagawa & Schielzeth in the absence of a 
# method for random slopes models), to illustrate the circumstances under
# which marginal and conditional R-squared can be expected to differ between
# the two approaches. 



# Load packages

  library(lme4)     # for fitting GLMMs
  library(lattice)  # for the xyplot function
  library(MuMIn)    # for the r.squaredGLMM function
  library(glmmADMB) # for the Owls data

# Remove all objects from memory

  rm(list = ls())

# Example 1

# The first example uses the orange tree growth data set that comes 
# installed with R.
# The orange tree data set has circumference measurements at 7 ages
# (in days) for 5 orange trees. For more details try ?Orange.

# First, convert age units to years, otherwise the lmer optimizer 
# complains

  Orange$ageYears <- Orange$age/365.25

# Plot tree growth trajectories

  xyplot(circumference ~ ageYears | Tree, Orange,
    panel=
      function(x,y){
        panel.xyplot(x,y)
        if(length(unique(x))>1) panel.abline(lm(y~x))
      })  


# The trees start at similar circumferences but grow at different rates, 
# so that slopes vary between trees, suggesting that a random slopes
# model might be suitable. A GLMM with Gaussian errors might fit reasonably
# well. (There are a number of reasons why this model is far from ideal, 
# including non-linearity, heteroskedasticity and possibly autocorrelated
# errors. However, the fit is not too bad, and the model illustrates the
# R-squared method nicely.)

# Fit random intercepts and random slopes models:
 
  orangemod.ri <- 
    lmer(circumference ~ ageYears + (1 | Tree), data = Orange)
  orangemod.rs <- 
    lmer(circumference ~ ageYears + (ageYears | Tree), data = Orange)

# Test the null hypothesis that the random slopes variance is zero.
# The correct p-value is (Stram & Lee, 1994, Biometrics, 50, 1171-1177):
   
  chistat <- anova(orangemod.ri, orangemod.rs)[2,"Chisq"] 
  0.5 * pchisq(chistat, 1, lower.tail=FALSE) + 0.5 * pchisq(chistat, 2, lower.tail=FALSE)

# The random slopes model fits much better.

# Extract the variance components require for the R-squared statistics,
# starting with the better-fitting random slopes model.

# First we need the design matrix (X), the number of observations (n)
# and the fixed effect estimates (Beta)

  X <- model.matrix(orangemod.rs)
  n <- nrow(X)
  Beta <- fixef(orangemod.rs)

# First the fixed effects variance (eqn 27 of Nakagawa & Schielzeth): 

  Sf <- var(X %*% Beta)

# Second, the list of covariance matrices of the random effects.
# Here the list has only a single element because there is only
# one level of random effects.

  (Sigma.list <- VarCorr(orangemod.rs))

# Use equation 11 in the paper to estimate the random effects variance 
# component.
# Using sapply ensures that the same code will work when there are multiple
# random effects (i.e. where length(Sigma.list) > 1)

  Sl <- 
    
      sapply(Sigma.list,
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          diags <- diag(Z %*% Sigma %*% t(Z)) # MAB MODIFICATION
          sum(diags)/n
        })


# As this model is an LMM, the additive dispersion variance, Se, is
# equivalent to the residual variance. The residual standard deviation
# is stored as an attribute of Sigma.list:
  
  Se <- attr(Sigma.list, "sc")^2

# Finally, the distribution-specific variance, Sd, is zero for LMMs:

  Sd <- 0

# Use eqns 29 and 30 from Nakagawa & Schielzeth to estimate marginal and
# conditional R-squared:

  total.var <- Sf + Sl + Se + Sd
  (Rsq.m <- Sf / total.var) 
  (Rsq.c <- (Sf + Sl) / total.var) 

# The fixed effects alone explain 81% of the variance. The combined fixed
# and random effects explain 97% of the variance. 

# Remove the R-squared ingredients to prevent their leaking into the
# next example:

  rm(X, n, Beta, Sf, Sigma.list, Sl, Se, Sd, total.var, Rsq.m, Rsq.c)

# Compare this to the R-squareds from the random intercepts model:

  X <- model.matrix(orangemod.ri)
  n <- nrow(X)
  Beta <- fixef(orangemod.ri)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(orangemod.ri)
  Sl <- 
    sum(
      sapply(Sigma.list,
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          sum(diag(Z %*% Sigma %*% t(Z)))/n
        }))
  Se <- attr(Sigma.list, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  (Rsq.m <- Sf / total.var) 
  (Rsq.c <- (Sf + Sl) / total.var) 

# The fixed effects alone now explain 82% of the variance. The combined
# fixed and random effects explain 93% of the variance. The difference
# between the two marginal R-squared values is small, which is expected 
# given that the design is balanced (see article discussion). 
# Allowing random slopes increased the conditional R-squared from 93% to
# 97%, an indication of how much better the random slopes model fits.
# NB: These comparisons are made purely for illustration. Although we'd
# expect the better fitting model to have the higher 
# conditional R-squared (but not necessarily higher marginal R-squared),
# we can't deduce better fit from a higher conditional R-squared, because
# the statistic is contains no penalty for adding parameters.

# This extension has been implemented in Kamil Barton's MuMIn package
# from version 1.10.0. Check the package version is >=1.10.0:

  installed.packages()["MuMIn","Version"] 

# The following lines:

  r.squaredGLMM(orangemod.rs)
  r.squaredGLMM(orangemod.ri)

# should return the same R^2[GLMM] estimates as above. 



############################################################################

# Remove all objects from memory

  rm(list = ls())

# Example 2

# The second example uses the Owls data (Roulin & Bersier, 2007,
# Animal Behaviour 74: 1099-1106) that comes with the glmmADMB package. 
# This data set records numbers of negotiation calls by owl chicks in 
# 27 nests. We are interested in investigating the relationship
# between the arrival time of adult owls at nests (average 22 visits
# per nest) and the total number of calls by the chicks in a 30 s period
# during the preceding 15 min. Each observation represents a single visit to a nest.
# The intention here is not to fit the best model, but to illustrate estimation of 
# the R-statistics from a random slopes model fitted on a real data set.
# For guidance on how to fit a more sensible model, see "Zero Inflated Models and
# Generalized Linear Mixed Models with R", 2012, Zuur, Saveliev & Ieno.

# First, create an observation-level factor for fitting an observation-level
# random effect. This will be used for fitting a random effect for mopping up
# overdispersion. 

  Owls$obs <- factor(1:nrow(Owls))

# Re-centre arrival time to units of hours after midnight.
# (This will the reduce the correlation between the random intercepts and slopes
# in the random slope model, which helps the fitting algorithm to converge.)

  Owls$ArrivalTime <- Owls$ArrivalTime - 24

# Plot the relationship between arrival time and the number of negotiation calls.

  xyplot(SiblingNegotiation ~ ArrivalTime | Nest, Owls,
    panel=
      function(x,y){
        panel.xyplot(x,y)
        if(length(unique(x))>1) panel.abline(lm(y~x))
      })  

# Based on fitting simple linear regression lines in each nest, a random slopes model
# looks plausible.

# Now fit the models.
# The random effects have a simple nested structure: visits are nested within 
# nests. Arrival time is be fitted as a continuous fixed effect. The number of calls
# is modelled as Poisson-distributed, with the inclusion of a lognormal observation-level
# random effect making the model a Poisson-lognormal GLMM, which is an alternative to
# negative binomial for overdispersed count data (Elston et al., 2001, Parasitology, 122, 563-9).


  owlmod.ri  <- 
    glmer(SiblingNegotiation ~ ArrivalTime + (1|Nest) + (1|obs),
      family="poisson", data=Owls, control=glmerControl(optimizer="bobyqa"))
  summary(owlmod.ri)

# There are substantial variances at the observation level (modelling overdispersion) 
# and at the nest level. The arrival time slope estimate is -0.151, which equates to a
  exp(-0.151) - 1 
# 14% per hour drop in call rate.

# Fit the random slopes model. 

  owlmod.rs  <- 
    glmer(SiblingNegotiation ~ ArrivalTime + (ArrivalTime|Nest) + (1|obs),
      family="poisson", data=Owls, control=glmerControl(optimizer="bobyqa"))
  summary(owlmod.rs)

# The overdispersion (obs) variance is slightly lower, suggesting that this model might explain 
# more variation than the random intercepts model. It is difficult to interpret the intercept
# and slope variance estimates, because they can change dramatically in response to centering.
# The global slope estimate is now -0.172, a 
  exp(-0.172) - 1
# 16% per hour reduction in call rate. The fixed effect estimates are quite similar so we would
# not expect marginal R-squared to differ very much between the two model fits.

# We want to know if the random slopes model is significantly better:
# The correct p-value here is (Stram & Lee, 1994, Biometrics, 50, 1171-1177):

  chistat <- anova(owlmod.ri, owlmod.rs)[2,"Chisq"] 
  0.5 * pchisq(chistat, 1, lower.tail=FALSE) + 0.5 * pchisq(chistat, 2, lower.tail=FALSE)

# The random slopes model is a better fit. 

# Calculate the R-squared statistics for the random slopes model.
# (Note that for the Poisson-lognormal GLMM the observation (overdispersion)
# variance belongs to the Se component so must be separated from the 
# other random effect variances, which are summed to give Sl.)

  X <- model.matrix(owlmod.rs)
  n <- nrow(X)
  Beta <- fixef(owlmod.rs)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(owlmod.rs)
  Sl <- 
    sum(
      sapply(Sigma.list["Nest"],
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          sum(diag(Z %*% Sigma %*% t(Z)))/n
        }))
  Se <- Sigma.list$obs
  Sd <- log(1 + 1/exp(mean(X %*% Beta))) ### See footnote ###
  (total.var <- c(Sf + Sl + Se + Sd))
  (Rsq.m <- c(Sf / total.var)) 
  (Rsq.c <- c((Sf + Sl) / total.var)) 

# The fixed effects explain 5.3% of the variance, while the 
# fixed plus random effects explain 26%.

  rm(X, n, Beta, Sf, Sigma.list, Sl, Se, Sd, total.var, Rsq.m, Rsq.c)

# Calculate the R-squared statistics for the random intercepts
# model. Remember that this model was significantly worse than
# the random slopes and intercepts model.

  X <- model.matrix(owlmod.ri)
  n <- nrow(X)
  Beta <- fixef(owlmod.ri)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(owlmod.ri)
  Sl <- 
    sum(
      sapply(Sigma.list["Nest"],
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          sum(diag(Z %*% Sigma %*% t(Z)))/n
        }))
  Se <- Sigma.list$obs
  Sd <- log(1 + 1/exp(mean(X %*% Beta)))
  total.var <- c(Sf + Sl + Se + Sd)
  (Rsq.m <- c(Sf / total.var)) 
  (Rsq.c <- c((Sf + Sl) / total.var)) 
  
# The fixed effects explain 4.3% of the variance, while the 
# fixed plus random effects explain about 17.6%. How should we interpret
# the difference between the two models in marginal R-squared? 
# It is slightly higher in the random slopes model simply because the 
# slope (Beta) estimate is slightly steeper, giving higher var(X %*% Beta),
# the denominator of marginal R-squared. The conditional R-squared is
# substantially higher in the random slopes model, principally because
# allowing each nest to have its own slope greatly increases
# the total random effect variance component, while reducing the residual
# (overdispersion) variance.

# The code above is only for illustration. The r.squaredGLMM function
# in MuMIn (version >= 1.10.0) gives the same results much more easily,
# while also working for packages other than lme4.

  r.squaredGLMM(owlmod.rs)  
  r.squaredGLMM(owlmod.ri)  

# Warnings may appear stating that the mean of the linear predictor is
# low (see footnote).


############################################################################

# Remove all objects from memory

  rm(list = ls())

# Example 3

# This example uses simulated data to illustrate a situation in which 
# misspecification of the random effects has a drastic effect on marginal 
# R-squared. As discussed in the article, for balanced designs we expect marginal 
# R-squared to be similar between random slopes and random intercepts 
# fits, even if the random slopes model is a much better fit.
# This simulated data set is an example of imbalance in numbers
# (but not covariate distributions) across groups. This example
# is contrived to illustrate the consequences of extreme imbalance for
# marginal R-squared. However ecological data are often imbalanced, and 
# less extreme but nevertheless substantial biases are plausible.

# Simulate data, starting with the unbalanced group sizes of 
# seven groups, a-g

  group.size <- 2^(2:8)
  names(group.size) <- letters[1:length(group.size)]
  group.size
  n.groups <- length(group.size)
  n <- sum(group.size)

# Create group ID

  id <- 
    factor(unlist(lapply(1:length(group.size), 
      function(i) rep(letters[i], group.size[i]))))
  table(id)

# Simulate random slopes from a normal distribution with mean=0 and
# SD=1, ordering them so that the largest group has the most positive slope

  set.seed(289688) # for repeatability
  rand.slopes <- sort(rnorm(n.groups, mean=0, sd=1))
  names(rand.slopes) <- levels(id)
  rand.slopes

# Create a covariate, x, randomly distributed across the groups
# to give a balanced distribution (so the imbalance is solely in 
# the group sizes)

  x <- sample(1:n)

# Simulate the response, y, as having random slopes varying around 
# a global slope of zero, with a residual SD of 20

  y <- sapply(1:length(x), function(i)x[i]*rand.slopes[id[i]]) + rnorm(length(x),0,20)

  plot(y ~ x, pch=as.character(id))

# Fit random intercept and random slopes models:

  simmod.ri <- lmer(y ~ x + (1 | id)) 
  summary(simmod.ri)
  simmod.rs <- lmer(y ~ x + (1 + x | id)) 
  summary(simmod.rs)

# The random slopes model clearly fits better:

  chistat <- anova(simmod.ri, simmod.rs)[2,"Chisq"] 
  0.5 * pchisq(chistat, 1, lower.tail=FALSE) + 0.5 * pchisq(chistat, 2, lower.tail=FALSE)

# Plot global and within-group predictions from each model

  lines(sort(x), predict(simmod.ri, newdata=data.frame(x=sort(x)), re.form=NA), lwd=3)
  invisible(lapply(levels(id), function(i) 
    lines(sort(x), predict(simmod.ri, newdata=data.frame(x=sort(x), id=i), re.form=NULL), lty=2)))

  lines(sort(x), predict(simmod.rs, newdata=data.frame(x=sort(x)), re.form=NA), col=2, lwd=3)
  invisible(lapply(levels(id), function(i) 
    lines(sort(x), predict(simmod.rs, newdata=data.frame(x=sort(x), id=i), re.form=NULL), col=2, lty=2)))

  abline(h=0, col=4, lwd=3)

  legend("topleft", 
    c("Random intercepts model global trend", "Random intercepts model local trends", 
      "Random slopes model global trend", "Random slopes model local trends", 
      "True global trend"),
      lty=c(1, 2, 1, 2, 1), lwd=c(3, 1, 3, 1, 3), col=c(1, 1, 2, 2, 4))

# In the misspecified random intercepts model the slope estimate (0.691, SE=0.013)
# is biased upward by the largest group, g. The estimate from the random slopes model 
# is closer to zero (0.104, SE=0.316). Notice also that the standard error from the 
# random slopes model is much higher, reflecting the much greater uncertainty about 
# the global slope due to the wide variation in slopes among groups.

# Estimate the R-squared statistics for the (correct) random slopes model

  X <- model.matrix(simmod.rs)
  n <- nrow(X)
  Beta <- fixef(simmod.rs)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(simmod.rs)
  Sl <- 
    sum(
      sapply(Sigma.list,
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          sum(diag(Z %*% Sigma %*% t(Z)))/n
        }))
  Se <- attr(Sigma.list, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  (Rsq.m <- Sf / total.var) 
  (Rsq.c <- (Sf + Sl) / total.var) 

# We are mainly interested here in the effect of the misspecification of the 
# model on marginal R-squared. The value from the correct, random slopes, model
# is 0.4%, close to zero, correctly reflecting the fact that the true slope is zero.

  rm(X, n, Beta, Sf, Sigma.list, Sl, Se, Sd, total.var, Rsq.m, Rsq.c)

# Estimate the R-squared statistics for the (incorrect) random intercepts model

  X <- model.matrix(simmod.ri)
  n <- nrow(X)
  Beta <- fixef(simmod.ri)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(simmod.ri)
  Sl <- 
    sum(
      sapply(Sigma.list,
        function(Sigma)
        {
          Z <-X[,rownames(Sigma)]
          sum(diag(Z %*% Sigma %*% t(Z)))/n
        }))
  Se <- attr(Sigma.list, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  (Rsq.m <- Sf / total.var) 
  (Rsq.c <- (Sf + Sl) / total.var) 

# As expected, the biased slope estimate from the random intercepts model has led to 
# an inflated marginal R-squared of 18%, much greater than the true value of zero.

# Alternatively:

  r.squaredGLMM(simmod.rs)  
  r.squaredGLMM(simmod.ri)  



############################################################################


############
# Footnote #
############

# Because the owl models are Poisson GLMMs, the distribution-specific variance (Sd)
# has to be estimated. This requires estimation of Beta0 (see eqn A6 in
# the appendix of Nakagawa & Schielzeth), which Nakagawa & Schielzeth
# recommend estimating as the intercept from a model with centred
# covariates. This is a little tricky for models with categorical
# covariates, and is anyway unnecessary as Beta0 is equivalent to the
# mean of the linear predictor, which is straightforward
# to calculate as mean(X %*% Beta), giving
# Sd = log(1 + 1/exp(mean(X %*% Beta)))
# This change is also in r.squaredGLMM in MuMIn version >= 1.10.0.

# However, note that the the approximation to Sd is expected to be
# unreliable for exp(Beta0) (which estimates lambda, the Poisson mean
# and variance) below about 3. This is because the variance of the log of
# a sample of counts is not defined when any count is zero, and zero
# counts will be frequent (> 5%) in any Poisson sample with lambda < 3.

# In Example 2, exp(Beta0) is 3.3 calls/nest, so Sd should be a reasonable approximation.
# However if the chick call rate per nest is divided by brood size (using an offset)
# to give a call rate per chick (as it would be in a proper model), exp(Beta) 
# is about 1 call/chick, and Sd could be highly inaccurate.
# It is unclear (and beyond the scope of this article) how to deal with
# this. One option could be to eliminate Sd altogether from the R-squared
# formulae, because it represents a component of variance that can never
# be explained. 

