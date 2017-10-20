## BRMs R2
# non-Bayesian. Inspired by sem.model.fits function

get_fixed.var #<- 

brm.model <- w.arth.abund.2012.brm
str(brm.model)
model.matrix(brm.model[[1]]$formula, data = filter(wind.arth.df, Year == "2012"))

as.data.frame(fixef(brm.model))$Estimate

terms.formula(brm.model[[1]]$formula,  data = filter(wind.arth.df, Year == "2012"))
