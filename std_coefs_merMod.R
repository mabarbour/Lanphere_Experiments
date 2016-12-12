

# from: Ben Bolker at https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q4/004974.html
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

REsdExtract(arth.abund.2013)

library(nlme)
data(Orthodont)
dat <- as.data.frame(Orthodont)
detach(package:nlme)

library(dplyr)
library(lme4)
fm2 <- lmer(distance ~ age + Sex + (age|Subject), data = dat)
lm.beta.lmer(fm2)

fm4 <- lmer(scale(distance) ~ (1|dat3$age.f) + scale(as.numeric(Sex)) + (age|Subject), data = dat)
fixef(fm4) # interesting, if you set "contr.sum" you effectively get the same coefficient, but it still isn't perfect...but neither is the interpretation of 
summary(fm4)
#lm.beta.lmer(fm4)

#For this example (which like yours has Sex, a two-level factor, as its only non-numeric predictor) we can show that we get the same answer (up to numeric fuzz) by scale()ing:
  
pdat <- with(dat,cbind(distance,age,s=as.numeric(Sex)))
pdat <- scale(pdat)
dat2 <- data.frame(pdat,Subject=dat$Subject)

fm3 <- lmer(distance ~ age + s + (age|Subject), data = dat2)
fixef(fm3)
summary(fm3) # Subject Int SD = 0.63; age:subject SD = 0.17

## change age to non-ordered categorical factor. Can I replicate the effect size of ~0.51?
dat3 <- dat2 %>% mutate(Subject = factor(Subject, ordered = FALSE),
                        age.f = factor(age, ordered = FALSE))

fm6 <- lmer(distance ~ age + s + (1|Subject), data = dat3)
fixef(fm6)
summary(fm6)

# contrasts sum gives you the same results as if you set intercept = -1
fm7 <- lmer(distance ~ (1|age.f) + s + (1|Subject), data = dat3)#, contrasts = list(age.f = "contr.sum"))
fixef(fm7)
summary(fm7)
mean(fixef(fm7)[2:4]) # contr.sum int = -0.24, contr.sum no int = 0, reg contr no int = 0, reg contr int = 0.84
sd(fixef(fm7)[2:4]) # sd = 0.42, 0.58. 0.58, 0.50
max(fixef(fm7)[2:4]) # max = 0.21, 0.707, 0.707, 1.33

## Let's just focus on effects of sex
fm5 <- lmer(distance ~ age + Sex + (1|Subject), data = dat)
lm.beta.lmer(fm5) # Sex = -0.39, age = 0.506

fm6 <- lmer(scale(distance) ~ age + Sex + (1|Subject), data = dat)
fixef(fm6) # # Sex = -0.39, age = 0.506

dat$age.f <- factor(dat$age)
fm7 <- lmer(scale(distance) ~ (1|age.f) + scale(as.numeric(Sex)) + (1|Subject), data = dat)
fixef(fm7) # Sex = -0.39
summary(fm7) # SD for age = 0.577

fm8 <- lmer(scale(distance) ~ (1|age.f) + (1|Sex) + (1|Subject), data = dat) 
summary(fm8) # SD for age = 0.577, SD for Sex = 0.53

fm9 <- lmer(distance ~ (1|age.f) + (1|Sex) + (1|Subject), data = dat) 
summary(fm9) # SD for age = 1.69, SD for sex = 1.55

fm10 <- lmer(distance ~ scale(age) + scale(as.numeric(Sex)) + (1|Subject), data = dat) 
summary(fm10) # SD for age = 1.48, SD for sex = 1.15

fm11 <- lmer(scale(distance) ~ age.f + scale(as.numeric(Sex)) + (1|Subject)-1, data = dat) 
fixef(fm11)
mean(fixef(fm11)[-5])
sd(fixef(fm11)[-5]) # SD = 0.58
summary(fm11)

fm12 <- lmer((distance) ~ age.f + scale(as.numeric(Sex)) + (1|Subject)-1, data = dat) 
fixef(fm12)
mean(fixef(fm12)[-5])
sd(fixef(fm12)[-5]) # SD = 1.71
summary(fm12)

model <- fm6
fixed.sd <- list()
for (i in seq_along(fixef(model))) {
  fixed.sd[[i]] <- data.frame(Factors = names(fixef(model))[i],
                               sd = sd(as.vector(fixef(model)[i] %*% t(model@pp$X)[i, ])))
}
plyr::ldply(fixed.sd) # according to this, the standard deviations used to calculate R2 are commensurate with beta-coefficients (at least when response is scaled too), therefore, they are more/less comparable with sd from random effects (if you believe that R2 are meaningful for mixed effect models)
