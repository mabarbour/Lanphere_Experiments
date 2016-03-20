## load required libraries ----
library(dplyr)
library(reshape2)
library(ggplot2)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
library(lmerTest)
#library(RLRsim)
library(car)
library(vegan)
library(mvabund)
source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')

## upload datasets ----

## ant-aphid: aphid growth rates
aa.aphid.GR <- read.csv("~/Documents/Lanphere_Experiments/final_data/ant_aphid_Aphis_popgrowth_df.csv") %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  select(Date_rel:plant_code, Aphis.growth.rate)
glimpse(aa.aphid.GR)

## ant-aphid: arthropod community
aa.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/ant_aphid_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(#Ants_all = ant_F_obscuripes + ant_black,
         Aphids_nonAphis = aphid_Tuberolachnus + aphid_LG,
         fact.Ant.mound.dist = as.factor(Ant.mound.dist),
         Plot_code = interaction(Block, fact.Ant.mound.dist)) %>%
  select(Block, Genotype, Ant.mound.dist, Aphid.treatment,
         plant_code, fact.Ant.mound.dist, Plot_code,
         aphid_Aphis, ant_F_obscuripes, ant_black, 
         Aphids_nonAphis,
         psyllid:sawfly_larva, syrphid_larva:LTF_Caloptilia)

aa.arth.names <- colnames(select(aa.arth.df,
                                 ant_F_obscuripes:LTF_Caloptilia))

# generate new columns for total abundance and richness
aa.total.abund <- rowSums(
  select(aa.arth.df, ant_F_obscuripes:LTF_Caloptilia))
aa.total.rich <- rowSums(
  select(aa.arth.df, ant_F_obscuripes:LTF_Caloptilia) > 0)
aa.herb.abund <- rowSums(
  select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,gall_R_salicisbattatus,grasshopper,leaftier_Tortricid,LTF_Caloptilia)
)
aa.herb.rich <- rowSums(
  select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,gall_R_salicisbattatus,grasshopper,leaftier_Tortricid,LTF_Caloptilia) > 0
)
aa.herb.conceal <- rowSums(
  select(aa.arth.df, gall_R_salicisbattatus,leaftier_Tortricid,LTF_Caloptilia)
)
aa.herb.nonconceal <- rowSums(
  select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,grasshopper)
)
aa.pred.abund <- rowSums(
  select(aa.arth.df, ant_black, spiders, syrphid_larva)
)
aa.total.pred.abund <- rowSums(
  select(aa.arth.df, ant_F_obscuripes, ant_black, spiders, syrphid_larva)
)
aa.pred.rich <- rowSums(
  select(aa.arth.df, ant_black, spiders, syrphid_larva) > 0
)
aa.arth.df <- mutate(aa.arth.df, 
                     total.abund = aa.total.abund,
                     total.rich = aa.total.rich,
                     herb.abund.nonAphis = aa.herb.abund,
                     herb.rich.nonAphis = aa.herb.rich,
                     herb.abund.conceal = aa.herb.conceal,
                     herb.abund.nonconceal = aa.herb.nonconceal,
                     pred.abund.nonFobs = aa.pred.abund,
                     pred.rich.nonFobs = aa.pred.rich,
                     pred.abund.all = aa.total.pred.abund)

# subset of data where plants had at least one arthropod individual
aa.arth.12.pos <- aa.arth.df %>%
  filter(total.abund > 0) 

# for avoiding pseudoreplication while testing for ant mound distance effect
aa.effect.12 <- aa.arth.12.pos %>%
  select(-plant_code, -Genotype, -Aphid.treatment) %>%
  group_by(Block, fact.Ant.mound.dist, Plot_code) %>%
  summarise_each(funs(mean))

hist(log(aa.arth.df$herb.abund.nonAphis+1))
hist(aa.arth.df$herb.rich.nonAphis)
hist(log(aa.arth.df$herb.abund.conceal+1))
hist(log(aa.arth.df$herb.abund.nonconceal+1))
hist(aa.arth.df$pred.abund.nonFobs)
hist(aa.arth.df$pred.rich.nonFobs)

## wind: arthropod community
# dead plants have already been removed
wind.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_arthropod_df.csv') %>%
  mutate(Block = as.factor(Block),
         Year = as.factor(Year),
         Plot_code = interaction(Block, Wind.Exposure)) 

wind.arth.names <- colnames(select(wind.arth.df, Gracilliaridae_miner:Spider)) # for subsetting community data

# generate new columns for total abundance and richness
total.abund <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Spider))
total.rich <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Spider) > 0)
herb.abund <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Aphididae))
herb.rich <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Aphididae) > 0)
pred.abund <- rowSums(
  select(wind.arth.df, Formica_ant, Spider))
pred.rich <- rowSums(
  select(wind.arth.df, Formica_ant, Spider) > 0)
wind.arth.df <- mutate(wind.arth.df, 
                       total.abund = total.abund,
                       total.rich = total.rich,
                       herb.abund = herb.abund,
                       herb.rich = herb.rich,
                       pred.abund = pred.abund,
                       pred.rich = pred.rich,
                       X = factor(seq(1,362,1)))

# 2012 dataset
# focus dataset on aggregated Family/Order arthropod groupings
w.arth.12.full <- wind.arth.df %>%
  filter(Year == "2012") %>%
  select(Block:plant_code, Plot_code,
         Gracilliaridae_miner:pred.rich) 

# subset of data where plants had at least one arthropod individual
w.arth.12.pos <- w.arth.12.full %>%
  filter(total.abund > 0) 

# for avoiding pseudoreplication while testing for wind effect
w.effect.12 <- w.arth.12.pos %>%
  select(-plant_code, -Genotype) %>%
  group_by(Block, Wind.Exposure, Plot_code) %>%
  summarise_each(funs(mean))

# 2013 dataset
# same structure as 2012 dataset
w.arth.13.full <- wind.arth.df %>%
  filter(Year == "2013") %>%
  select(Block:plant_code, Plot_code, 
         Gracilliaridae_miner:pred.rich) 

w.arth.13.pos <- w.arth.13.full %>%
  filter(total.abund > 0)

w.effect.13 <- w.arth.13.pos %>%
  select(-plant_code, -Genotype) %>%
  group_by(Block, Wind.Exposure, Plot_code) %>%
  summarise_each(funs(mean))

## Wind community analyses ----

## Wind: arthropod abundance analysis ----

# GLMM
arth.abund.glmer <- glmer(total.abund ~ Wind.Exposure*Year +
                            (1|Genotype) + 
                             (1|Block) + 
                             (1|Block:Wind.Exposure) +
                             (1|plant_code) +
                             (1|X),
                           data = wind.arth.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            #Genotype = "contr.sum",
                                            Year = "contr.sum"),
                           family = "poisson",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))
summary(arth.abund.glmer)
overdisp_fun(arth.abund.glmer) 
plot(arth.abund.glmer) 

plotFEsim(FEsim(arth.abund.glmer))

var.calc(update(arth.abund.glmer, .~. -(1|X))) # removed individual-level random effect, because this is already calculated in var.calc

# Likelihood ratio tests
w.arth.3 <- drop1(arth.abund.glmer, test = "Chisq") %>% tidy()
w.arth.2 <- drop1(update(arth.abund.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq") %>% tidy()
w.arth.1 <- drop1(update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Genotype + Year), test = "Chisq") %>% tidy()

arth.abund.anova <- bind_rows(w.arth.1, w.arth.2, w.arth.3) %>%
  filter(term != "<none>") %>%
  mutate(Response = "arthropod abundance",
         Response.type = "community",
         Specific.type = "arthropods",
         Test.type = "Chi-square",
         Experiment = "wind",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2 for significant predictors
arth.abund.main <- update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Genotype + Year)

sem.model.fits(list(arth.abund.glmer, arth.abund.main))

arth.abund.R2 <- 
  data.frame(Factor = c("Wind.Exposure","Genotype","Year"),
             Response = "arthropod abundance",
             Response.type = "community",
             Specific.type = "arthropods",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(arth.abund.main, update(arth.abund.main, .~. -Wind.Exposure)),
                                deltaR2(arth.abund.main, update(arth.abund.main, .~. -Genotype)),
                                deltaR2(arth.abund.main, update(arth.abund.main, .~. -Year))),2))

## Wind: herbivore abundance analysis ----

# GLMM
herb.abund.glmer <- lme4::glmer(herb.abund ~ Wind.Exposure*Year +
                            (1|Genotype) +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|X) +
                            (1|plant_code),
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(herb.abund.glmer)
plotFEsim(FEsim(herb.abund.glmer))

var.calc(update(herb.abund.glmer, .~. -(1|X)))

overdisp_fun(herb.abund.glmer) 
plot(herb.abund.glmer) 

# Likelihood ratio tests
#w.herb.3 <- drop1(herb.abund.glmer, test = "Chisq") %>% tidy() # marginal effect, but note that full model had difficulty converging
w.herb.2 <- drop1(herb.abund.glmer, test = "Chisq") %>% tidy() # 2-way model runs without warnings.
w.herb.1 <- drop1(update(herb.abund.glmer, .~. -Wind.Exposure:Year), test = "Chisq") %>% tidy()

herb.abund.anova <- bind_rows(w.herb.1, w.herb.2, w.herb.3) %>%
  filter(term != "<none>") %>%
  mutate(Response = "herbivore abundance",
         Response.type = "community",
         Specific.type = "arthropods",
         Test.type = "Chi-square",
         Experiment = "wind",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2 for significant predictors
herb.abund.main <- update(herb.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Genotype + Year)

sem.model.fits(list(herb.abund.glmer, herb.abund.main))

herb.abund.R2 <- 
  data.frame(Factor = c("Wind.Exposure","Genotype","Year"),
             Response = "herbivore abundance",
             Response.type = "community",
             Specific.type = "arthropods",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(herb.abund.main, update(herb.abund.main, .~. -Wind.Exposure)),
                                deltaR2(herb.abund.main, update(herb.abund.main, .~. -Genotype)),
                                deltaR2(herb.abund.main, update(herb.abund.main, .~. -Year))),2))

## Wind: predator abundance analysis ----

# GLMM
pred.abund.glmer <- glmer(pred.abund ~ Wind.Exposure*Year +
                            (1|Genotype) + # data don't support more complex interactions with Genotype.
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_code),
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           #Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(pred.abund.glmer)
overdisp_fun(pred.abund.glmer) # no overdispersion
plot(pred.abund.glmer) 

deltaR2(pred.abund.glmer, update(pred.abund.glmer, .~. -Wind.Exposure))
sem.model.fits(pred.abund.glmer)[ ,"Marginal"] - sum(test$var_percent[1:4])
test <- var.calc(pred.abund.glmer)
  
plotFEsim(FEsim(pred.abund.glmer))
plotREsim(REsim(pred.abund.glmer))


# Likelihood ratio tests
w.pred.3 <- drop1(pred.abund.glmer, test = "Chisq") %>% tidy() # marginal effect, but note that full model had difficulty converging
w.pred.2 <- drop1(update(pred.abund.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq") %>% tidy() # still convergence issues
w.pred.1 <- drop1(update(pred.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Genotype + Year), test = "Chisq") %>% tidy() # runs without any convergence issues.

pred.abund.anova <- bind_rows(w.pred.1, w.pred.2, w.pred.3) %>%
  filter(term != "<none>") %>%
  mutate(Response = "predator abundance",
         Response.type = "community",
         Specific.type = "arthropods",
         Test.type = "Chi-square",
         Experiment = "wind",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2 for significant predictors
pred.abund.main <- update(pred.abund.glmer, .~. -Wind.Exposure*Year + Wind.Exposure + Year)

sem.model.fits(list(pred.abund.glmer, pred.abund.main))

# Genotype is not significant, but has a higher R2 value than wind exposure...
pred.abund.R2 <- 
  data.frame(Factor = c("Wind.Exposure","Genotype","Year"),
             Response = "predator abundance",
             Response.type = "community",
             Specific.type = "arthropods",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(pred.abund.main, update(pred.abund.main, .~. -Wind.Exposure)),
                                deltaR2(pred.abund.main, update(pred.abund.main, .~. -Genotype)),
                                deltaR2(pred.abund.main, update(pred.abund.main, .~. -Year))),2))

# tests
pred.1 <- glmer(pred.abund ~ Wind.Exposure*Year +
                            (1|Genotype) +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_code),
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(pred.1) 
sem.model.fits(pred.1)
var.rand <- 0.48183+0.09968
overdisp_fun(pred.1) # no evidence of overdispersion
sem.model.fits(pred.1) # about 4% for wind, but now, Genotype is showing up as no effect.

plot(pred.abund ~ Genotype, wind.arth.df)


## Wind: Caloptilia analysis ----

# GLMM
LTF.abund.glmer <- glmer(Gracilliaridae_miner ~ Wind.Exposure*Year +
                            (1|Genotype) + # data don't support more complex interactions with Genotype.
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_code),
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           #Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(LTF.abund.glmer)
overdisp_fun(LTF.abund.glmer) # no overdispersion
plot(LTF.abund.glmer) 

var.calc(LTF.abund.glmer)
anova(LTF.abund.glmer, update(LTF.abund.glmer, .~. -(Wind.Exposure*Year|Genotype) + (Wind.Exposure + Year|Genotype)), test = "Chisq")
anova(update(LTF.abund.glmer, .~. -(Wind.Exposure*Year|Genotype) + (Wind.Exposure + Year|Genotype)), update(LTF.abund.glmer, .~. -(Wind.Exposure*Year|Genotype) + (1|Genotype)), test = "Chisq")
anova(update(LTF.abund.glmer, .~. -(Wind.Exposure*Year|Genotype) + (1|Genotype)), update(LTF.abund.glmer, .~. -(Wind.Exposure*Year|Genotype)), test = "Chisq")

plotREsim(REsim(LTF.abund.glmer))

## Wind: Cecidomyiidae galler analysis ----
with(filter(wind.arth.df, Year == "2013"), interaction.plot(Wind.Exposure, Genotype, Cecidomyiidae_gall))

# GLMM
gall.abund.glmer <- glmer(Cecidomyiidae_gall ~ Wind.Exposure +
                           (1|Genotype) + # data don't support more complex interactions with Genotype.
                           (1|Block) + 
                           (1|Block:Wind.Exposure), #+
                           #(1|X) +
                           #(1|plant_code),
                         data = filter(wind.arth.df, Year == "2013"),
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(gall.abund.glmer)
overdisp_fun(gall.abund.glmer) # no overdispersion
plot(gall.abund.glmer) 

plotREsim(REsim(gall.abund.glmer))

anova(gall.abund.glmer, update(gall.abund.glmer, .~. -(1|Genotype)))
var.calc(gall.abund.glmer)

## 2012 ----
summary(w.arth.12.full$total.abund)
summary(w.arth.12.full$total.rich)
sum(colSums(w.arth.12.full[ ,wind.arth.names]))
round(colSums(w.arth.12.full[ ,wind.arth.names])/sum(colSums(w.arth.12.full[ ,wind.arth.names]))*100,0)

# total abundance
plot(total.rich ~ total.abund, w.arth.12.full)
hist(w.arth.12.full$total.abund)
hist(w.arth.12.full$total.rich)
summary(w.arth.12.full$total.rich)
summary(w.arth.12.full$total.abund)
colSums(wind.arth.df[ ,wind.arth.names])/sum(colSums(wind.arth.df[ ,wind.arth.names]))


## arthropod probability 
hist(wind.arth.df$total.abund)
with(wind.arth.df, interaction.plot(Wind.Exposure, Genotype, total.abund))
library(piecewiseSEM)

### Test random vs. fixed effect specification
w.total.12 <- glmer(total.abund ~ Wind.Exposure*Genotype +
                      (1|Year) +
                      (1|plant_code) +
                      #(1|X) +
                     (1|Block/Wind.Exposure), 
                   wind.arth.df, family = poisson(link = log), contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"),
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e5)))
w.total.12.main <- update(w.total.12, .~. -Wind.Exposure:Genotype)
w.total.12.geno <- update(w.total.12.main, .~. -Wind.Exposure)
w.total.12.wind <- update(w.total.12.main, .~. -Genotype)

# fixed
mods <- list(w.total.12, w.total.12.main, w.total.12.geno, w.total.12.wind)
sem.model.fits(mods)

# random
rw.total.12 <- glmer(total.abund ~ (1|Wind.Exposure) + (1|Genotype) +
                      (1|Year) +
                      (1|plant_code) +
                      #(1|X) +
                      (1|Block/Wind.Exposure), 
                    wind.arth.df, family = poisson(link = log), contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"),
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=2e5)))
#rw.total.12.main <- update(w.total.12, .~. -(Wind.Exposure|Genotype))
rw.total.12.geno <- update(rw.total.12, .~. -(1|Wind.Exposure))
rw.total.12.wind <- update(rw.total.12, .~. -(1|Genotype))

# Generate null model (intercept and random effects only, no fixed effects)
null.model = update(rw.total.12, formula = paste(". ~ ", get.random.formula(rw.total.12, "~1", modelList = NULL)))

# Get the fixed effects of the null model
null.fixef = as.numeric(fixef(null.model))

varDist = log(1 + 1/exp(null.fixef))
obs = names(unlist(lapply(ranef(rw.total.12), nrow))[unlist(lapply(ranef(rw.total.12), nrow)) == nrow(rw.total.12@pp$X)]) ##????
varDisp =  sum(
  
  sapply(VarCorr(rw.total.12)[obs], function(Sigma) {
    
    X = model.matrix(rw.total.12)  
    
    Z = X[, rownames(Sigma)]
    
    sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X)
    
  } )

summary(rw.total.12)
0.13491/(1.09009+0.08281+0.06778+0.13491+0.02928+0.05382+varDist) # genotype at 9%
0.05382/(1.09009+0.08281+0.06778+0.13491+0.02928+0.05382+varDist) # wind at 3%
r.mods <- list(rw.total.12, #w.total.12.main, 
               rw.total.12.geno, rw.total.12.wind)
sem.model.fits(r.mods)


#sem.model.fits()
anova(w.total.12,
      update(w.total.12, .~. -Wind.Exposure:Genotype))#,
main <- update(w.total.12, .~. -Wind.Exposure:Genotype)
anova(main,
      update(main, .~. -Wind.Exposure))
anova(main,
      update(main, .~. -Genotype))
sem.model.fits(list(w.total.12, 
                    update(w.total.12, .~. -Wind.Exposure*Genotype + Wind.Exposure),
                    update(w.total.12, .~. -Wind.Exposure*Genotype + Genotype)))#, update(w.total.12, .~.-(1|Genotype) + (1|Wind.Exposure) + Genotype)))
sem.model.fits(list(w.total.12, update(w.total.12, .~. -(1|Genotype)+Genotype)))
overdisp_fun(w.total.12)
summary(w.total.12)
exp(0.327)
exp(-0.25256)
exp(-0.25256) # exposed plants receive 22% less arthropods per year
exp(0.19992)
exp(0.3463)
exp(0.5121) # unexposed plants receive 1.7x more arthropods per year than unexposed plants
predict(w.total.12, newdata = newdata)
newdata <- with(wind.arth.df, expand.grid(Wind.Exposure=unique(Wind.Exposure), Year=unique(Year)))
0.04716/(0.04716+0.05898+0.01119+0.07207+0.66929+0.25) # 4.2% of the variance
visreg(w.total.12)
plot(w.total.12)
exp(0.2428)

anova(w.total.12, update(w.total.12, .~. -(1|Genotype)))
anova(w.total.12, update(w.total.12, .~. -(1|Genotype)))

# test G and GxE
anova(w.total.12,
      update(w.total.12, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), # possible G effect
      update(w.total.12, .~. -(Wind.Exposure|Genotype)))

# test wind effect
anova(w.total.12, update(w.total.12, .~. -Wind.Exposure))

## Herbivore probability
with(w.arth.12.full, sum((herb.abund > 0))/length(herb.abund)) # ~49% of plants with a herbivore
w.herb.12 <- glmer((herb.abund > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                     (1|Block/Wind.Exposure), 
                   w.arth.12.full, 
                   family = "binomial")
summary(w.herb.12)
visreg(w.herb.12)

# test G and GxE
anova(w.herb.12,
      update(w.herb.12, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), # possible G effect
      update(w.herb.12, .~. -(Wind.Exposure|Genotype)))

# test wind effect
anova(w.herb.12, update(w.herb.12, .~. -Wind.Exposure))

# random model
w.herb.12.rand <- glmer((herb.abund > 0) ~ (1|Wind.Exposure) + (1|Genotype) + (1|Block/Wind.Exposure), 
                   w.arth.12.full, 
                   family = "binomial")
anova(w.herb.12.rand, update(w.herb.12.rand, .~. -(1|Wind.Exposure))) # interesting, now wind exposure doesn't apparently have a detectable effect...
anova(w.herb.12.rand, update(w.herb.12.rand, .~. -(1|Genotype)))
summary(w.herb.12.rand)
0.24/(0.1479+0.24+0.3976+(pi^2)/3) # 6% for genotype
0.1479/(0.1479+0.24+0.3976+(pi^2)/3) # 4% for wind exposure

## aphid probability
with(w.arth.12.full, interaction.plot(Wind.Exposure, Genotype, Aphididae))
with(w.arth.12.full, sum((Gracilliaridae_miner > 0))/length(Gracilliaridae_miner)) # ~21% of plants with a herbivore
hist(w.arth.12.full$Gracilliaridae_miner)
w.miner.12 <- lmer(log(Aphididae+1) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                     (1|Block/Wind.Exposure), 
                   w.arth.12.full)
summary(w.miner.12)
visreg(w.miner.12)

# test G and GxE
anova(w.miner.12,
      update(w.miner.12, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), # possible G effect
      update(w.miner.12, .~. -(Wind.Exposure|Genotype)))

# test wind effect
anova(w.miner.12, update(w.miner.12, .~. -Wind.Exposure))

# random model
w.miner.12.rand <- glmer((Gracilliaridae_miner>0) ~ (1|Wind.Exposure) + (1|Genotype) + (1|Block/Wind.Exposure), 
                        w.arth.12.full, 
                        family = "binomial")
anova(w.miner.12.rand, update(w.miner.12.rand, .~. -(1|Wind.Exposure))) # interesting, now wind exposure doesn't apparently have a detectable effect...
anova(w.miner.12.rand, update(w.miner.12.rand, .~. -(1|Genotype)))
summary(w.miner.12.rand)
0.5547/(0+0+0.5547+0.2356+(pi^2)/3) # 14% for genotype
0.2356/(0+0+0.5547+0.2356+(pi^2)/3) # 6% for wind exposure

#w.herbabund.12 <- lmer(log(herb.abund) ~ Wind.Exposure + (1|Genotype) + 
 #                        (1|Block/Wind.Exposure), 
  #                     filter(w.arth.12.full, herb.abund > 0))
#summary(w.herbabund.12)
#plot(w.herbabund.12)

# predator probability. Insufficient cases to likely even test for an effect.
with(w.arth.12.full, sum((pred.abund > 0))/length(pred.abund)) # only ~10% of plants had a predator
w.pred.12 <- glmer((pred.abund > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                     (1|Block/Wind.Exposure), 
                   w.arth.12.full, 
                   family = "binomial")
summary(w.pred.12)
visreg(w.pred.12)

# test G and GxE
anova(w.pred.12,
      update(w.pred.12, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), # possible G effect
      update(w.pred.12, .~. -(Wind.Exposure|Genotype)))

# test wind effect
anova(w.pred.12, update(w.pred.12, .~. -Wind.Exposure))


# dissimilarity matrix, full dataset
w.dis.12 <- vegdist(w.arth.12.pos[ ,wind.arth.names],
                    method = "horn")
w.dis.12.sub <- vegdist(w.effect.12[ ,wind.arth.names],
                        method = "horn")

# Testing interaction: block level as strata
# marginally significant
adonis(w.dis.12 ~ Wind.Exposure*Genotype,
       data = w.arth.12.pos,
       strata = w.arth.12.pos$Block)
# non-significant -> meets assumptions of adonis
anova(betadisper(d = w.dis.12,
                 group = with(w.arth.12.pos, 
                              interaction(Wind.Exposure,
                                          Genotype)),
                 bias.adjust = TRUE),
      strata = w.arth.12.pos$Block)

# Testing Genotype: plot level as strata
# significant effect
adonis(w.dis.12 ~ Genotype,
       data = w.arth.12.pos,
       strata = w.arth.12.pos$Plot_code)
meandist(dist = w.dis.12,
         grouping = w.arth.12.pos$Genotype) 
# non-significant -> meets assumptions of adonis
anova(betadisper(d = w.dis.12,
                 group = w.arth.12.pos$Genotype,
                 bias.adjust = TRUE),
      strata = w.arth.12.pos$Plot_code)

# Testing Wind exposure: Block level as strata on subset
# non-significant effect
adonis(w.dis.12.sub ~ Wind.Exposure,
       data = w.effect.12,
       strata = w.effect.12$Block)
# non-significant: meets assumptions of adonis
anova(betadisper(d = w.dis.12.sub,
                 group = w.effect.12$Wind.Exposure,
                 bias.adjust = TRUE),
      strata = w.effect.12$Block)

## 2013

# Overall, arthropod abundance and richness is very small (median = 1 for both) on each plant. This indicates that the primary source of variation is in the presence or absence of arthropods, rather than their abundance.
summary(w.arth.13.full$total.abund)
summary(w.arth.13.full$total.rich)

# Tortricidae (54%), Cecidomyiids (19%), Spiders (9%), and Gracilliaridae (8%) were the most abundant members of the community. Therefore, I'm going to look at how these different groups responded to genotype and wind exposure (in addition to total arthropods)
round(colSums(w.arth.13.full[ ,wind.arth.names])/sum(colSums(w.arth.13.full[ ,wind.arth.names]))*100,0)

## Total arthropods
# plot the data
with(w.arth.13.full, interaction.plot(Wind.Exposure, Genotype, total.abund)) 

# GLMM with binomial distribution
w.total.13 <- glmer((total.abund > 0) ~ Wind.Exposure + 
                      (1|Genotype) +
                     (1|Block), 
                   w.arth.13.full,
                   family = "binomial")
summary(w.total.13)
exp(0.9010) # the odds of finding a herbivore on an unexposed plant increase 2.5-fold
plot(w.total.13)

# No evidence of a G or GxE effect
anova(w.total.13,
      update(w.total.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), 
      update(w.total.13, .~. -(Wind.Exposure|Genotype)))

# Wind exposure influenced the probabilty of an arthropod colonizing the host-plant.
anova(w.total.13, update(w.total.13, .~. -Wind.Exposure))



# herbivore probability
with(w.arth.13.full, sum((herb.abund > 0))/length(herb.abund)) # ~55% of plants with a herbivore
w.herb.13 <- glmer((herb.abund > 0) ~ Wind.Exposure + Genotype + #+ (Wind.Exposure|Genotype) +
                     (1|Block/Wind.Exposure), 
                   w.arth.13.full, 
                   family = "binomial")
summary(w.herb.13)
visreg(w.herb.13)

# test G and GxE
anova(w.herb.13,
      update(w.herb.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)), 
      update(w.herb.13, .~. -(Wind.Exposure|Genotype)))

anova(w.herb.13, update(w.herb.13, .~. -Wind.Exposure))
anova(w.herb.13, update(w.herb.13, .~. -Genotype))

w.herbabund.13 <- lmer(log(herb.abund) ~ Wind.Exposure + (Wind.Exposure|Genotype) + 
                         (1|Block/Wind.Exposure), 
                       filter(w.arth.13.full, herb.abund > 0))
summary(w.herbabund.13)
plot(w.herbabund.13)

# predator probability. 
with(w.arth.13.full, sum((pred.abund > 0))/length(pred.abund)) # only ~19% of plants had a predator

w.pred.13 <- glmer((pred.abund > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
        (1|Block/Wind.Exposure), 
      w.arth.13.full, 
      family = "binomial")
summary(w.pred.13)

# testing G and GxE
anova(w.pred.13, 
      update(w.pred.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)),
      update(w.pred.13, .~. -(Wind.Exposure|Genotype)))

visreg(w.pred.13)

# leaftier probability. 
with(w.arth.13.full, sum((Tortricidiae_leaftier > 0))/length(Tortricidiae_leaftier)) # only ~35% of plants had a predator

w.leaftier.13 <- glmer((Tortricidiae_leaftier > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                     (1|Block:Wind.Exposure), 
                   w.arth.13.full, 
                   family = "binomial", contrasts = list(Wind.Exposure = "contr.sum"))
summary(w.leaftier.13)
plot(w.leaftier.13)

# testing G and GxE
anova(w.leaftier.13, 
      update(w.leaftier.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)),
      update(w.leaftier.13, .~. -(Wind.Exposure|Genotype)))

visreg(w.pred.13)
plot(Tortricidiae_leaftier ~ Genotype, w.arth.13.full)

# Cecidomyiid gall probability. 
with(w.arth.13.full, sum((Cecidomyiidae_gall > 0))/length(Cecidomyiidae_gall)) # only ~18% of plants had a predator

w.gall.13 <- glmer((Cecidomyiidae_gall > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                         (1|Block/Wind.Exposure), 
                       w.arth.13.full, 
                       family = "binomial", contrasts = list(Wind.Exposure = "contr.sum"))
summary(w.gall.13)
plot(w.gall.13)

# testing G and GxE
anova(w.gall.13, 
      update(w.gall.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)),
      update(w.gall.13, .~. -(Wind.Exposure|Genotype)))

visreg(w.gall.13)
plot(Cecidomyiidae_gall ~ Wind.Exposure, w.arth.13.full)

# Cecidomyiid gall probability. 
with(w.arth.13.full, sum((Gracilliaridae_miner > 0))/length(Gracilliaridae_miner)) # only ~8.5% of plants had a predator

w.miner.13 <- glmer((Gracilliaridae_miner > 0) ~ Wind.Exposure + (Wind.Exposure|Genotype) +
                     (1|Block/Wind.Exposure), 
                   w.arth.13.full, 
                   family = "binomial", contrasts = list(Wind.Exposure = "contr.sum"))
summary(w.miner.13)
plot(w.miner.13)

# testing G and GxE
anova(w.miner.13, 
      update(w.miner.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)),
      update(w.miner.13, .~. -(Wind.Exposure|Genotype)))

visreg(w.miner.13)
plot(Gracilliaridae_miner ~ Genotype, w.arth.13.full)

# predator probability. 
with(w.arth.13.full, sum((pred.abund > 0))/length(pred.abund)) # only ~19% of plants had a predator

w.miner.13 <- glmer((pred.abund > 0) ~ Wind.Exposure + (1|Genotype) +
                      (1|Block/Wind.Exposure), 
                    w.arth.13.full, 
                    family = "binomial", contrasts = list(Wind.Exposure = "contr.sum"))
summary(w.miner.13)
plot(w.miner.13)

# testing G and GxE
anova(w.miner.13, 
      update(w.miner.13, .~. -(Wind.Exposure|Genotype) + (1|Genotype)),
      update(w.miner.13, .~. -(Wind.Exposure|Genotype)))

anova(w.miner.13, 
      update(w.miner.13, .~. -(1|Genotype)))

visreg(w.miner.13)
plot(Gracilliaridae_miner ~ Genotype, w.arth.13.full)


# dissimilarity matrix, full dataset
w.dis.13 <- vegdist(w.arth.13.pos[ ,wind.arth.names],
                    method = "horn")
w.dis.13.sub <- vegdist(w.effect.13[ ,wind.arth.names],
                        method = "horn")

# Testing interaction: block level as strata
# non-significant effect
adonis(w.dis.13 ~ Wind.Exposure*Genotype,
       data = w.arth.13.pos,
       strata = w.arth.13.pos$Block)
# non-significant -> meets assumptions of adonis
anova(betadisper(d = w.dis.13,
                 group = with(w.arth.13.pos, 
                              interaction(Wind.Exposure,
                                          Genotype)),
                 bias.adjust = TRUE),
      strata = w.arth.13.pos$Block)

# Testing Genotype: plot level as strata
# non-significant effect
adonis(w.dis.13 ~ Genotype,
       data = w.arth.13.pos,
       strata = w.arth.13.pos$Plot_code)
# non-significant -> meets assumptions of adonis
anova(betadisper(d = w.dis.13,
                 group = w.arth.13.pos$Genotype,
                 bias.adjust = TRUE),
      strata = w.arth.13.pos$Plot_code)

# Testing Wind exposure: Block level as strata on subset
# significant effect
adonis(w.dis.13.sub ~ Wind.Exposure,
       data = w.effect.13,
       strata = w.effect.13$Block)
meandist(dist = w.dis.13.sub,
         grouping = w.effect.13$Wind.Exposure) # 39% dissimilarity
# non-significant: meets assumptions of adonis
anova(betadisper(d = w.dis.13.sub,
                 group = w.effect.13$Wind.Exposure,
                 bias.adjust = TRUE),
      strata = w.effect.13$Block)

## Ant-aphid community analyses ----
# only 2012
summary(aa.arth.df$total.abund)
summary(aa.arth.df$total.rich)
round(colSums(aa.arth.df[ ,aa.arth.names])/sum(colSums(aa.arth.df[ ,aa.arth.names]))*100,0)
sum(colSums(aa.arth.df[ ,aa.arth.names]))

sum(colSums(aa.arth.df[ ,aa.arth.names[-3]])) # excluding non-Aphis aphid counts

# herbivore probability and abundance
with(aa.arth.df, sum((herb.abund.nonAphis > 0))/length(herb.abund.nonAphis)) # 73% of plants with a herbivore

aa.herb.prob <- glmer((herb.abund.nonAphis > 0) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
      data = aa.arth.df, family = "binomial",
      contrasts = list(Aphid.treatment = "contr.sum")) # convergence issues so I dropped most effects in random effect model
plot(aa.herb.prob)
summary(aa.herb.prob)
visreg(aa.herb.prob, xvar = "Ant.mound.dist", by = "Aphid.treatment")

# testing G effect
anova(aa.herb.prob,
      update(aa.herb.prob,.~. -(1|Genotype))) 

aa.herb.prob.rand <- glmer((herb.abund.nonAphis > 0) ~ (1|Ant.mound.dist) + (1|Aphid.treatment) + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
                           data = aa.arth.df, family = "binomial")
summary(aa.herb.prob.rand)

with(aa.arth.df, sum((herb.abund.nonconceal > 0))/length(herb.abund.nonconceal)) # 61% of plants
aa.nonconceal.prob.glmer <- glmer((herb.abund.nonconceal > 0) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
      data = aa.arth.df, family = "binomial",
      contrasts = list(Aphid.treatment = "contr.sum"))
summary(aa.nonconceal.prob.glmer)
visreg(aa.nonconceal.prob.glmer, xvar = "Ant.mound.dist", by = "Aphid.treatment")

with(aa.arth.df, sum((herb.abund.conceal > 0))/length(herb.abund.conceal)) # 39% of plants
aa.conceal.prob.glmer <- glmer((herb.abund.conceal > 0) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
                                  data = aa.arth.df, family = "binomial",
                                  contrasts = list(Aphid.treatment = "contr.sum"))
summary(aa.conceal.prob.glmer)
visreg(aa.conceal.prob.glmer, xvar = "Ant.mound.dist", by = "Aphid.treatment")

#aa.herb.abund <- lmer(log(herb.abund.nonAphis) ~ Aphid.treatment*Ant.mound.dist + (0+Aphid.treatment|Genotype) + (1|Block/fact.Ant.mound.dist),
 #                     data = filter(aa.arth.df, herb.abund.nonAphis > 0)) #,contrasts = list(Genotype = "contr.sum")
#plot(aa.herb.abund)
#summary(aa.herb.abund)
#anova(aa.herb.abund, ddf = "Kenward-Roger")
#visreg(aa.herb.abund, xvar = "Ant.mound.dist", by = "Aphid.treatment")
#visreg(aa.herb.abund, xvar = "Genotype", by = "Aphid.treatment")

# predator probability and abundance. no F obscuripes
with(aa.arth.df, sum((pred.abund.nonFobs > 0))/length(pred.abund.nonFobs)) # 60% of plants with a predator

aa.pred.prob <- glmer((pred.abund.nonFobs > 0) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
                      data = aa.arth.df, family = "binomial",
                      contrasts = list(Aphid.treatment = "contr.sum"))
plot(aa.pred.prob)
summary(aa.pred.prob)

# test G
anova(aa.pred.prob, update(aa.pred.prob, .~. -(1|Genotype))) # possibly a significant effect.

# F_obscuripes
aa.Fobs.prob <- glmer((ant_F_obscuripes > 0) ~ scale(Ant.mound.dist) + Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
                          data = aa.arth.df, family = "binomial",
                          contrasts = list(Aphid.treatment = "contr.sum"))
plot(aa.Fobs.prob)
summary(aa.Fobs.prob)

## all predators
with(aa.arth.df, sum((ant_F_obscuripes > 0))/length(ant_F_obscuripes)) # only 8% of plants with F_obscuriptes
with(aa.arth.df, sum((pred.abund.all > 0))/length(pred.abund.all)) 

aa.pred.all.prob <- glmer((pred.abund.all > 0) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
                      data = aa.arth.df, family = "binomial",
                      contrasts = list(Aphid.treatment = "contr.sum"))
plot(aa.pred.all.prob)
summary(aa.pred.all.prob)

# test G
anova(aa.pred.all.prob, update(aa.pred.all.prob, .~. -(1|Genotype))) # possibly a significant effect.

visreg(aa.pred.all.prob, xvar = "Ant.mound.dist", by = "Aphid.treatment")

summary(glmer((pred.abund.all > 0) ~ (1|Aphid.treatment) + (1|Ant.mound.dist) + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
      data = aa.arth.df, family = "binomial"))

#aa.pred.abund <- lmer(log(pred.abund.nonFobs) ~ scale(Ant.mound.dist)*Aphid.treatment + (1|Genotype) + (1|Block/fact.Ant.mound.dist),
 #                     data = filter(aa.arth.df, pred.abund.nonFobs > 0),
  #                    contrast = list(Aphid.treatment = "contr.sum"))
#plot(aa.pred.abund)
#summary(aa.pred.abund)
#anova(aa.pred.abund, ddf = "Kenward-Roger")
#visreg(aa.pred.prob, xvar = "Ant.mound.dist")
#0.01407/(0.2206+0.01407)

# dissimilarity matrix, full dataset
aa.dis.12 <- vegdist(aa.arth.12.pos[ ,aa.arth.names],
                    method = "horn")
aa.dis.12.sub <- vegdist(aa.effect.12[ ,aa.arth.names],
                        method = "horn")

# Testing 3-way interaction: block level as strata
# non-significant
adonis(aa.dis.12 ~ Ant.mound.dist*Aphid.treatment*Genotype,
       data = aa.arth.12.pos,
       strata = aa.arth.12.pos$Block)
# non-significant -> meets assumptions of adonis
anova(betadisper(d = aa.dis.12,
                 group = with(aa.arth.12.pos, 
                              interaction(Ant.mound.dist,
                                          Aphid.treatment,
                                          Genotype)),
                 bias.adjust = TRUE),
      strata = aa.arth.12.pos$Block)

# Testing 2-way interactions: block level as strata
# non-significant
adonis(aa.dis.12 ~ (Ant.mound.dist + Aphid.treatment + Genotype)^2,
       data = aa.arth.12.pos,
       strata = aa.arth.12.pos$Block)

# Testing Genotype and Aphid Treatment effects: plot level as strata
# significant effect
adonis(aa.dis.12 ~ Genotype + Aphid.treatment,
       data = aa.arth.12.pos,
       strata = aa.arth.12.pos$Plot_code)

# significant -> does NOT meet assumptions of adonis. Interesting...the graphical plot appears okay...
disper.Aphid <- betadisper(d = aa.dis.12,
                 group = aa.arth.12.pos$Aphid.treatment,
                 bias.adjust = TRUE)
anova(disper.Aphid, strata = aa.arth.12.pos$Plot_code)
boxplot(disper.Aphid) # appears to meet assumptions...

# significant -> does NOT meet assumptions of adonis. Interesting...the graphical plots appears okay...
disper.Geno <- betadisper(d = aa.dis.12,
                 group = aa.arth.12.pos$Genotype,
                 bias.adjust = TRUE)
anova(disper.Geno, strata = aa.arth.12.pos$Plot_code)
boxplot(disper.Geno) # appears to meet assumptions...

# ants and Caloptilia appear to be driving effects
plot(capscale(aa.arth.12.pos[ ,aa.arth.names] ~ 
                Genotype + Aphid.treatment,
              data = aa.arth.12.pos, 
              distance = "horn"),
     display = c("cn","sp"))

# ants
hist(aa.arth.12.pos$Ants_all)
ant.lmer <- glmer(Ants_all/total.abund ~ Genotype*Aphid.treatment +
                   (1|Block) + (1|Block:fact.Ant.mound.dist),
                 data = aa.arth.12.pos,
                 weights = total.abund,
                 family = "binomial")
summary(ant.lmer)
plot(ant.lmer) # better but not great
anova(ant.lmer, ddf = "Kenward-Roger")

# Caloptilia
hist(aa.arth.12.pos$LTF_Caloptilia)
LTF.lmer <- lmer(log(LTF_Caloptilia+1) ~ Genotype*Aphid.treatment +
                   (1|Block) + (1|Block:fact.Ant.mound.dist),
                 aa.arth.12.pos)
summary(LTF.lmer)
plot(LTF.lmer) # better but not great
anova(LTF.lmer, ddf = "Kenward-Roger") # only genotype effect


# Testing Ant mound distance: Block level as strata on subset
# non-significant effect
adonis(aa.dis.12.sub ~ Ant.mound.dist,
       data = aa.effect.12,
       strata = aa.effect.12$Block)
# non-significant: meets assumptions of adonis
anova(betadisper(d = aa.dis.12.sub,
                 group = aa.effect.12$fact.Ant.mound.dist,
                 bias.adjust = TRUE),
      strata = aa.effect.12$Block)


## old ----
adonis(wind.effect.13[ ,wind.arth.names] ~ Wind.Exposure, 
       data = wind.effect,
       strata = wind.effect$Block,
       method = "horn")

inter.13 <- with(w.arth.13.pos, interaction(Wind.Exposure,Block))
adonis(decostand(w.arth.13.pos[ ,wind.arth.names], method = "hellinger") ~ Genotype*Wind.Exposure, 
       data =w.arth.13.pos,
       strata = w.arth.13.pos$Block,
       method = "euclidean")

adonis(decostand(w.arth.13.pos[ ,wind.arth.names], method = "hellinger") ~ Genotype, 
       data =w.arth.13.pos,
       strata = inter.13,
       method = "euclidean")




adonis(decostand(wind.effect[ ,wind.arth.names], method = "hellinger") ~ Wind.Exposure, 
       data = wind.effect,
       strata = wind.effect$Block,
       method = "euclidean")

colSums(w.arth.12.full[ ,wind.arth.names])/
  sum(colSums(w.arth.12.full[ ,wind.arth.names])) # dominant groups: aphids, leaf miners, leaf tiers, and spiders

## Analysis

# Total abundance. Had to run model without interaction term, because otherwise a singularity appeared.
sum(table(w.arth.12.full$total.abund)[-1]) # 97 sites with 1 or more arthropods, whereas the other half have zero.

## First test with binomial model
w.abund.glmer <- glmer((total.abund > 0) ~ Genotype + Wind.Exposure +
                       (1|Block) + (1|Block:Wind.Exposure),
                     data = w.arth.12.full,
                     family = "binomial",
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl=list(maxfun=2e4)))
overdisp_fun(w.abund.glmer) # no evidence of overdispersion
plot(w.abund.glmer) 
tt <- getME(w.abund.glmer,"theta")
ll <- getME(w.abund.glmer,"lower")
min(tt[ll==0]) # singularity appears with interaction

summary(w.abund.glmer) # much greater probability of finding an arthropod on unexposed vs. exposed plants.

# second test with only positive values. Tried glmer variations, but the poisson model and one that included individual-level random effects were horribly overdispersed, and the residuals were terrible, so I decided to just use the log(x+1) transformed lmer
w.abund.lmer <- lmer(log(total.abund+1) ~ Genotype*Wind.Exposure +
                        (1|Block) + (1|Block:Wind.Exposure),
                      data = w.arth.12.pos)
plot(w.abund.lmer) # okay, but not great
anova(w.abund.lmer, ddf = "Kenward-Roger")

# total richness with only positive values
plot((total.rich) ~ Genotype, w.arth.12.pos)
w.rich.lmer <- lmer(log(total.rich+1) ~ Genotype*Wind.Exposure +
                       (1|Block) + (1|Block:Wind.Exposure),
                     data = w.arth.12.pos)
plot(w.rich.lmer) # okay, but not great
anova(w.rich.lmer, ddf = "Kenward-Roger")
visreg::visreg(w.rich.lmer, xvar = "Genotype")

# Thoughts, test for an effect of block. If there is none, remove it from the model to save residual df.
# test manyglm
mv.w.arth.12 <- mvabund(w.arth.12.pos[ ,wind.arth.names])

test <- manyglm(mv.w.arth.12 ~ Block + Genotype*Wind.Exposure,
                data = w.arth.12.pos,
                family = "negative.binomial")
plot(test)
anova(test, p.uni = "unadjusted")
test1 <- manyglm(mv.w.arth.12 ~ Genotype + Wind.Exposure,
                data = w.arth.12.pos,
                family = "negative.binomial")
test2 <- manyglm(mv.w.arth.12 ~ Wind.Exposure,
                 data = w.arth.12.pos,
                 family = "negative.binomial")
plot(test) # looks okay
anova.manyglm(test2, test1) # marginal sig. clear effect of wind, clear effect of Genotype.

w.Grac.12 <- glm((Gracilliaridae_miner > 0) ~ Block + Wind.Exposure + Genotype,
                  w.arth.12.pos,
                  family = "binomial")
summary(w.Grac.12)
anova(w.Grac.12, test = "LR")
plot(w.Grac.12)
anova(w.Grac.12, type = 3, ddf = "Kenward-Roger")


w.Aphid.12 <- glmer(Gracilliaridae_miner/total.abund ~ Genotype + Wind.Exposure + (1|Block/Wind.Exposure),
                    family = "binomial",
                    weights = total.abund,
                  data = w.arth.12.pos, control = glmerControl(optimizer = "bobyqa",
                                                               optCtrl=list(maxfun=2e4)))
summary(w.Aphid.12)
plot(w.Aphid.12)
te <- ranef(w.Aphid.12)$`Block:Wind.Exposure`
te2 <- ranef(w.Aphid.12)$'Wind.Exposure:Block'
ggQQ_ranef(te[,1])
anova(w.Aphid.12, update(w.Aphid.12, .~. -Genotype))
anova(w.Aphid.12, test = "LR")
anova(w.Aphid.12, type = 3, ddf = "Kenward-Roger")
visreg::visreg(w.Aphid.12, xvar = "Wind.Exposure", by = "Genotype", scale = "response")

test.aphid <- glm((Aphididae>0) ~ Wind.Exposure + Genotype,
                           data = w.arth.12.pos, family = "binomial")
summary(test.aphid)
plot(test.aphid)
anova(test.aphid)

# Community dissimilarity using all positive values
w.arth.12.hell <- decostand(w.arth.12.pos[ ,wind.arth.names],
                            method = "hellinger")

plot(log(w.arth.12.pos$Gracilliaridae_miner+1) ~ w.arth.12.pos$Genotype)
plot(log(w.arth.12.pos$Aphididae+1) ~ w.arth.12.pos$Genotype)
plot(log(w.arth.12.pos$Aphididae+1) ~ w.arth.12.pos$Wind.Exposure)

inter <- with(w.arth.12.pos, interaction(Block, Wind.Exposure))
w.rda.12 <- rda(w.arth.12.hell ~ Genotype,
                data = w.arth.12.pos)
summary(w.rda.12)
plot(w.rda.12, display = c("sp","cn"))
anova(w.rda.12, by = "margin", strata = inter)

w.arth.12.bin <- w.arth.12.pos[ ,wind.arth.names]
w.arth.12.bin[w.arth.12.bin > 0] <- 1

wind.2012.adonis <- adonis(vegdist(filter(w.arth.12.pos, Genotype != "J")[ ,wind.arth.names], binary = FALSE, method = "bray") ~ Genotype*Wind.Exposure, data = filter(w.arth.12.pos, Genotype != "J"), strata = filter(w.arth.12.pos, Genotype != "J")$Block)
wind.2012.adonis

#wind.12.cap <- capscale(wind.comm.2012 ~ Genotype + Wind.Exposure + Condition(Block), data = wind.arth.2012, distance = "bray")
#summary(wind.12.cap)
#plot(wind.12.cap, display = c("cn","sp"))

anova(betadisper(vegdist(filter(w.arth.12.pos)[ ,wind.arth.names], binary = FALSE, method = "horn"),
                 group = filter(w.arth.12.pos)$Genotype,
                 bias.adjust = TRUE),
      permutations = how(blocks = inter))
boxplot(betadisper(vegdist(w.arth.12.bin, method = "jaccard"), 
                 group = w.arth.12.pos$Genotype,
                 bias.adjust = TRUE))

wind.2012.rda <- rda(decostand(wind.comm.2012,
                               method = "hellinger") ~ 
                       Wind.Exposure +
                       Condition(Block),
                     data = wind.arth.2012)
summary(wind.2012.rda)
plot(wind.2012.rda, display = c("cn","sp"))

# note that using wind.arth.2012 greatly reduces the size of the dataset.
wind.arth.2012.lmer <- wind.arth.df %>%
  select(Year:plant_code, Gracilliaridae_miner:Spider, wind.arth.nonzero, wind.arth.rich) %>%
  filter(Year == 2012)

rich.12 <- glmer((wind.arth.rich>0) ~ Wind.Exposure*Genotype + 
                  (1|Block/Wind.Exposure),
                family = "binomial",
                wind.arth.2012.lmer)
summary(rich.12)
anova(rich.12, ddf = "Kenward-Roger")
plot(rich.12)
plot((wind.arth.rich) ~ Wind.Exposure + Genotype, wind.arth.2012.lmer)

arth.abund.12 <- lmer(log(wind.arth.nonzero+1) ~ Wind.Exposure*Genotype + 
                        (1|Block/Wind.Exposure),
                      wind.arth.2012.lmer)
anova(arth.abund.12, ddf = "Kenward-Roger")
plot(arth.abund.12)
plot(wind.arth.nonzero ~ Wind.Exposure, wind.arth.2012.lmer)

hist(log(wind.arth.2012$Gracilliaridae_miner+1))
w.Grac.12 <- lmer(log(Gracilliaridae_miner+1) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.Grac.12)
plot(w.Grac.12)
anova(w.Grac.12, type = 3, ddf = "Kenward-Roger")
plot(log(Gracilliaridae_miner+1) ~ Genotype, wind.arth.2012)
table(wind.arth.2012$Genotype)
plot(log(Gracilliaridae_miner+1) ~ Wind.Exposure, wind.arth.2012)

w.Grac.12.glmer <- glmer((Gracilliaridae_miner>0) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure), family = "binomial",
                  wind.arth.2012.lmer)
summary(w.Grac.12.glmer)
plot(w.Grac.12.glmer)
anova(w.Grac.12.glmer, update(w.Grac.12.glmer, .~. -Wind.Exposure))
anova(w.Grac.12, type = 3, ddf = "Kenward-Roger")

w.spid.12 <- lmer(log(Spider+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.spid.12)
anova(w.spid.12, type = 3, ddf = "Kenward-Roger")
plot(w.spid.12)

w.Aphid.12 <- lmer(log(Aphididae+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2012)
summary(w.Aphid.12)
anova(w.Aphid.12, type = 3, ddf = "Kenward-Roger")
plot(w.Aphid.12)
plot(log(Aphididae+1) ~ Wind.Exposure, wind.arth.2012)

w.Aphid.12 <- glmer(Aphididae ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                    family = "poisson",
                   wind.arth.2012.lmer)
summary(w.Aphid.12)
anova(w.Aphid.12, type = 3, ddf = "Kenward-Roger")
plot(w.Aphid.12)

w.Tort.12 <- lmer(log(Tortricidiae_leaftier+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.Tort.12)
anova(w.Tort.12, type = 3, ddf = "Kenward-Roger")
plot(w.Tort.12)

anova(betadisper(vegdist(wind.comm.2012), 
                 group = wind.arth.2012$Wind.Exposure,
                 bias.adjust = TRUE))
anova(betadisper(vegdist(wind.comm.2012), 
                 group = wind.arth.2012$Genotype,
                 bias.adjust = TRUE)) # violates assumption of adonis...

# 2013
wind.arth.2013 <- wind.arth.df %>%
  select(Year:plant_code, Gracilliaridae_miner:Spider) %>%
  filter(Year == 2013, wind.arth.nonzero > 0)

wind.comm.2013 <- wind.arth.2013 %>%
  select(Gracilliaridae_miner:Spider) 
colSums(wind.comm.2013)/sum(colSums(wind.comm.2013)) # dominant groups: leaf miners, leaf tiers, galls, and spiders

wind.2013.adonis <- adonis(wind.comm.2013 ~ Genotype*Wind.Exposure, data = wind.arth.2013, strata = wind.arth.2013$Block)
wind.2013.adonis # interesting... only an effect of wind exposure

w.spid.13 <- lmer(log(Spider+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.spid.13)
anova(w.spid.13, type = 3, ddf = "Kenward-Roger")
plot(w.spid.13)

hist(wind.arth.2013$Cecidomyiidae_gall)
w.Cecid.13 <- lmer(log(Cecidomyiidae_gall+1) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.Cecid.13)
anova(w.Cecid.13, type = 3, ddf = "Kenward-Roger")
plot(w.Cecid.13)
plot(log(Cecidomyiidae_gall+1) ~ Wind.Exposure + Genotype, wind.arth.2013)

w.Cecid.13.glmer <- glmer((Cecidomyiidae_gall>0) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure), family = "binomial",
                   wind.arth.2013)
summary(w.Cecid.13.glmer)

hist(log(wind.arth.2013$Tortricidiae_leaftier+1))
w.Tort.13 <- lmer(log(Tortricidiae_leaftier+1) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.Tort.13)
anova(w.Tort.13, type = 3, ddf = "Kenward-Roger")
plot(w.Tort.13)
plot(log(Tortricidiae_leaftier+1) ~ Wind.Exposure, wind.arth.2013)

w.Tort.13.glmer <- glmer((Tortricidiae_leaftier>0) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure), family = "binomial",
                  wind.arth.2013)
summary(w.Tort.13)
anova(w.Tort.13, type = 3, ddf = "Kenward-Roger")

plot(betadisper(vegdist(wind.comm.2013), 
                 group = wind.arth.2013$Wind.Exposure,
                 bias.adjust = TRUE)) # violates assumption of adonis
plot(betadisper(vegdist(wind.comm.2013), 
                 group = wind.arth.2013$Genotype,
                 bias.adjust = TRUE)) # doesn't violates assumption of adonis...


## ant-aphid growth rates ----
# no relationship with distance to ant mound
ggplot(aa.aphid.GR, aes(x = Genotype, y = Aphis.growth.rate, color = Genotype)) +
  geom_boxplot() +
  geom_point(aes(fill = Genotype), position = position_jitterdodge(jitter.width = 2)) + 
  #stat_smooth(se = FALSE) +
  facet_wrap(~Date_rel, nrow = 2)

# residuals look a bit binomial in their distribution. I think this is because the different durations put different lower limits to the decline in population growth rate. Now, I'm thinking the best thing to do would be to just used Aphid densities over time. It is a lot more of an intuitive response variable, and I could account for zeros with a poisson or maybe zero-inflated poisson.
aphid.GR.lmer <- lmer(Aphis.growth.rate ~ Ant.mound.dist + Genotype*Date_rel + (1|Block/Ant.mound.dist) + (1|plant_code),
                      aa.aphid.GR)
summary(aphid.GR.lmer)
#confint(profile(aphid.GR.lmer))
plot(aphid.GR.lmer) # things that need to be accounted for.
anova(aphid.GR.lmer, ddf = "Kenward-Roger")


#mutate(non.leaftier.herb.abund = LTF_Caloptilia + tentmine_Phyllonorycter + gall_R_rigidae + gall_R_salicisbrassicoides + gall_Pontania + gall_Aculus + leafhopper_C_reductus + leafhopper_green + leafhopper_unk + sawfly_larva + caterpillar_looper + caterpillar_LB + caterpillar_unk + red_scale + psyllid + grasshopper,
#      total.herb.abund = non.leaftier.herb.abund + leaftier_Tortricid,
#     total.omniv.abund = stinkbug + ant_F_obscuripes,
#    total.pred.abund = spider_Theridion + spider_BY + spider_NW + spider_Tetragnathid + spider_CS + spider_Larionoides)

#herbs <- colnames(wind.2013.vis.df.max)[5:21]
#omnivs <- colnames(wind.2013.vis.df.max)[22:23]
#preds <- colnames(wind.2013.vis.df.max)[24:29]

#wind.2013.max.gg <- wind.2013.vis.df.max %>%
# filter(Dead < 1) %>%
#gather(Species, Abundance, LTF_Caloptilia:total.pred.abund)

#library(ggplot2)
#ggplot(filter(wind.2013.max.gg, Species %in% c("leaftier_Tortricid", "non.leaftier.herb.abund", "total.pred.abund")),
#      aes(x = Genotype, y = Abundance, color = Wind.Exposure)) +
#geom_boxplot() + 
#facet_wrap(~Species, nrow = 3)

# graph doesn't look too convincing. Aphid population growth is virtually always negative...
ggplot(filter(aa.2012.vis.aphidgrowth, Aphid.Treatment == "aphid"),
       aes(x = Genotype, y = Aphid.growth.rate, color = Genotype)) +
  geom_boxplot() +
  facet_wrap(~relative.Date, ncol = 2)

ggplot(filter(aa.2012.vis.alive, Aphid.Treatment == "aphid"),
       aes(x = total.aphid.abund, y = Red.Ants)) +
  geom_point() +
  stat_smooth(method = "loess")

library(lme4)
aphid.lmer <- glmer(Red.Ants.pres ~ total.aphid.abund + (1|Genotype) + (1|relative.Date) + (1|Block/Distance.to.Ant.Mound), aa.2012.vis.alive, family = "binomial")
summary(aphid.lmer)
confint(profile(aphid.lmer))
anova(aphid.lmer)
plot(aphid.lmer)
plot(total.aphid.abund ~ relative.Date, aa.2012.vis.alive)


