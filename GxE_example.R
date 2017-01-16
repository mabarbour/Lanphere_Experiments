
## Simulate data for GxE

set.seed(1)

# Genotype variation in the control (n = 50)
rel.c <- rnorm(50, 1, sd = 0.2)

# Genetic variation of tolerance, i.e. slopes (var = sd^2)
# The population effect is set to zero for simplicity
gen.var <- rnorm(50, 0, sd = 0.30)

# Genotype means for damage treatment
rel.d <- rel.c + gen.var

# Set within genotype:treatment combination sample size
n.in = 20

# Generate within genotype:treatment data with sigma = 0.3
dat.con <- rnorm(50*n.in, rel.c, sd = 0.3)
dat.treat <- rnorm(50*n.in, rel.d, sd = 0.3)

# Give genotype id
gen.id <- rep(1:50, times = n.in)
# Put together data frame
dat <- data.frame( "gen" = factor(c(gen.id, gen.id)), "treat" = rep(c("con","treat"), each = 50*n.in),
                   "rel.fit" = c(dat.con, dat.treat) )

## Random intercept and slope

# Fitting a model with random intercepts and random slope but no co-varaince between intercepts and slope
library(lme4)

mod.lmer <- lme4::lmer(rel.fit ~ treat + (treat| gen), data = dat)
summary(mod.lmer)
var.calc(mod.lmer)
plotREsim(REsim(mod.lmer))

# This model can also be defined as estimating a variance for each group (control, treat) instead of estimating the differences as in the model above
mod2.lmer <- lme4::lmer(rel.fit ~ treat + (1|gen) + (0 + treat| gen), data = dat)
summary(mod2.lmer)
var.calc(mod2.lmer)
plotREsim(REsim(mod2.lmer))

## Correlated intercept and slope
# We simulate some data
set.seed(1)
rel.c <- sort(rnorm(50, 1, sd = 0.2))
gen.var <- rnorm(50, 0, sd = 0.30)
rel.d <- rel.c + sort(gen.var,  decreasing = TRUE) 
n.in = 20
dat.con <- rnorm(50*n.in, rel.c, sd = 0.3)
dat.treat <- rnorm(50*n.in, rel.d, sd = 0.3)
gen.id <- rep(1:50, times = n.in)
dat2 <- data.frame( "gen" = factor(c(gen.id, gen.id)), "treat" = rep(c("con","treat"), each = 50*n.in),
                    "rel.fit" = c(dat.con, dat.treat) )

# plot the data
with(dat2, interaction.plot(treat, gen, rel.fit))

# We fit the same models again and look at the output

# in lme4
mod2.lme4 <- lmer(rel.fit ~ treat + (0+treat|gen), dat2)
summary(mod2.lme4)

mod2.lme4.corr <- lmer(rel.fit ~ treat + (1+treat|gen), dat2)
summary(mod2.lme4.corr)
var.calc(mod2.lme4.corr)

mod2.lme4.null <- lmer(rel.fit ~ treat + (1|gen), dat2)
summary(mod2.lme4.null)

anova(mod2.lme4.corr, mod2.lme4.null)
