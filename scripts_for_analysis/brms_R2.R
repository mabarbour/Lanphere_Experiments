## BRMs R2
# non-Bayesian. Inspired by sem.model.fits function

get_fixed.var #<- 
library(broom)
wind.data <- filter(wind.arth.df, Year == "2013")
brm.model <- w.arth.rich.2012.brm

bayes_R2(w.arth.rich.2012.brm)


bayes_R2(update(w.arth.rich.2012.brm, .~. -Wind.Exposure))

tidy.brm <- tidy(brm.model)

int <- filter(tidy.brm, term == "b_Intercept")$estimate

bWind <- filter(tidy.brm, term == "b_Wind.ExposureUnexposed")$estimate
up95.Wind <- filter(tidy.brm, term == "b_Wind.ExposureUnexposed")$upper
low95.Wind <- filter(tidy.brm, term == "b_Wind.ExposureUnexposed")$lower

wind.fixef.design <- tidy(model.matrix(~Wind.Exposure, data = wind.data))

avg.wind <- sd(int*wind.fixef.design$X.Intercept. + bWind*wind.fixef.design$Wind.ExposureUnexposed)
upper.wind <- sd(int*wind.fixef.design$X.Intercept. + up95.Wind*wind.fixef.design$Wind.ExposureUnexposed)
lower.wind <- sd(int*wind.fixef.design$X.Intercept. + low95.Wind*wind.fixef.design$Wind.ExposureUnexposed)

avg.wind
upper.wind
lower.wind

pp_check(brm.model)

library(bayesplot)
ppc_intervals_grouped(y = filter(wind.arth.df, Year == "2013")$total.rich,
                      yrep = posterior_predict(brm.model),
                      x = wind.fixef.design$Wind.ExposureUnexposed,
                      group = filter(wind.arth.df, Year == "2013")$Genotype,
                      prob = 0.5)

scale(wind.fixef.design$Wind.ExposureUnexposed)

cbind.data.frame(tidy.brm$term, round(tidy.brm[ ,-1],2))

ggplot(filter(wind.arth.df, Year == "2012"), aes(x = Genotype, y = total.rich, group = Wind.Exposure, color = Wind.Exposure)) + 
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = .25))

ggplot(filter(wind.arth.df, Year == "2012"), aes(x = Wind.Exposure, y = total.rich, group = Genotype, color = Genotype)) + 
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = .25))

ggplot(filter(wind.arth.df, Year == "2013"), aes(x = Genotype, y = total.rich, group = Wind.Exposure, color = Wind.Exposure)) + 
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = .25))

avg.wind^2/(avg.wind^2 + sum((tidy.brm$estimate[c(3:6,8)])^2))

(avg.wind^2 + sum((tidy.brm$estimate[c(3:6)])^2))/(avg.wind^2 + sum((tidy.brm$estimate[c(3:6,8)])^2))
  