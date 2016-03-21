## load required libraries ----
library(dplyr)
library(reshape2)
library(ggplot2)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
#library(lmerTest)
#library(RLRsim)
library(lme4)
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

## Wind community analyses

## Wind: arthropod abundance analysis ----

# GLMM
arth.abund.glmer <- glmer(total.abund ~ Wind.Exposure*Year*Genotype +
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
summary(arth.abund.glmer)
overdisp_fun(arth.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(arth.abund.glmer) 

# Likelihood ratio tests
arth.abund.anova <- anova.table(arth.abund.glmer, test = "Chisq", experiment = "wind")

## Calculate R2 for significant predictors
arth.abund.up <- update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure*Year + (1|Genotype) - (1|X))

arth.abund.R2 <- var.table(arth.abund.up, experiment = "wind")

## Wind: arthropod richness analysis ----
plot(total.rich ~ total.abund, wind.arth.df)

# GLMM
arth.rich.glmer <- glmer(total.rich ~ Wind.Exposure*Year*Genotype +
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
print(summary(arth.rich.glmer), correlation = TRUE)
overdisp_fun(arth.rich.glmer) # no overdispersion
plot(arth.rich.glmer) 

arth.rich.anova <- anova.table(arth.rich.glmer, test = "Chisq", experiment = "wind")

arth.rich.up <- update(arth.rich.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

arth.rich.R2 <- var.table(arth.rich.up, experiment = "wind")

## Wind: rarefied arthropod richness analysis ----
wind.arth.df$total.rarerich <- rarefy(wind.arth.df[ ,wind.arth.names], 2) - 1

hist(filter(wind.arth.df, total.abund > 2)$total.rarerich)

# GLMM
arth.rarerich.glmer <- lmer(total.rarerich ~ Wind.Exposure*Year+ Wind.Exposure*Genotype +
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           #(1|X) +
                           (1|plant_code),
                         data = filter(wind.arth.df, total.abund > 2),
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"))
print(summary(arth.rarerich.glmer), correlation = TRUE)
overdisp_fun(arth.rarerich.glmer) # no overdispersion
plot(arth.rarerich.glmer) 

arth.rarerich.anova <- anova.table(arth.rarerich.glmer, test = "F", experiment = "wind")

arth.rarerich.up <- update(arth.rarerich.glmer, .~. -Wind.Exposure*Year -Wind.Exposure*Genotype + Wind.Exposure -(1|plant_code) -(1|Block:Wind.Exposure))

plotFEsim(FEsim(arth.rarerich.up))

arth.rarerich.R2 <- var.table(arth.rarerich.up, experiment = "wind")

## Wind: herbivore abundance analysis ----
hist(wind.arth.df$herb.abund)

# GLMM
herb.abund.glmer <- glmer(herb.abund ~ (Wind.Exposure + Year + Genotype)^2 + # solves the problem with the model being identified
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
print(summary(herb.abund.glmer), correlation = TRUE) 
overdisp_fun(herb.abund.glmer) 
plot(herb.abund.glmer) 

# Wald Chi-square test
herb.abund.anova <- anova.table(herb.abund.glmer, test = "Chisq", experiment = "wind")


## Calculate R2 for significant predictors
herb.abund.up <- update(herb.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + (1|Genotype) - (1|X)) # removing individual-level effect, because the R2 calculation already accounts for overdispersion.

herb.abund.R2 <- var.table(herb.abund.up, "wind")

## Wind: predator abundance analysis ----

# GLMM
pred.abund.glmer <- glmer(pred.abund ~ Wind.Exposure*Year + Genotype +
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
summary(pred.abund.glmer)
overdisp_fun(pred.abund.glmer) # no overdispersion
plot(pred.abund.glmer) 

# Wald Chi-square tests
pred.abund.anova <- anova.table(pred.abund.glmer, test = "Chisq", experiment = "wind")

## Calculate R2 for significant predictors
pred.abund.up <- update(pred.abund.glmer, .~. -Wind.Exposure:Year -Genotype -(1|Block)) # dropped block because 0% of variance explained

pred.abund.R2 <- var.table(pred.abund.up, experiment = "wind")

## Wind: Caloptilia analysis ----
hist(wind.arth.df$Gracilliaridae_miner)
with(wind.arth.df, table(Gracilliaridae_miner, Year))

# GLMM
LTF.abund.glmer <- glmer(Gracilliaridae_miner ~ Wind.Exposure*Year + 
                           Genotype + # only model without convergence
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
summary(LTF.abund.glmer)
overdisp_fun(LTF.abund.glmer) # no overdispersion
plot(LTF.abund.glmer) 

## Wald Chi-square test
LTF.abund.anova <- anova.table(LTF.abund.glmer, test = "Chisq", experiment = "wind")

## calculate variance components
LTF.abund.up <- update(LTF.abund.glmer, .~. -Wind.Exposure:Year
                       -Genotype + (1|Genotype))
LTF.abund.R2 <- var.table(LTF.abund.up, experiment = "wind")

## Wind: Cecidomyiidae galler analysis ----
with(wind.arth.df, table(Cecidomyiidae_gall, Year))

# GLMM
gall.abund.glmer <- glmer(Cecidomyiidae_gall ~ Wind.Exposure+Genotype + Year +
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
summary(gall.abund.glmer)
overdisp_fun(gall.abund.glmer) # no overdispersion
plot(gall.abund.glmer) 

gall.abund.anova <- anova.table(gall.abund.glmer, test = "Chisq", experiment = "wind")

gall.abund.up <- update(gall.abund.glmer, .~. -Genotype + (1|Genotype) -(1|Block))

gall.abund.R2 <- var.table(gall.abund.up, experiment = "wind")

## Wind: Spider analysis ----
spider.abund.glmer <- glmer(Spider ~ Wind.Exposure*Year + Genotype +
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
summary(spider.abund.glmer)
overdisp_fun(spider.abund.glmer) # no overdispersion
plot(spider.abund.glmer) 

spider.abund.anova <- anova.table(spider.abund.glmer, test = "Chisq", experiment = "wind")

spider.abund.up <- update(spider.abund.glmer, .~. -Wind.Exposure:Year -Genotype -(1|Block))

spider.abund.R2 <- var.table(spider.abund.up, experiment = "wind")

## Wind: Aphididae analysis ----
with(wind.arth.df, table(Aphididae, Year)) # no aphids in 2013, therefore, I only analyze 2012

aphid.abund.glmer <- glmer(Aphididae ~ Wind.Exposure + Genotype +
                              (1|Block) + 
                              (1|Block:Wind.Exposure) +
                              #(1|X) +
                              (1|plant_code),
                            data = filter(wind.arth.df, Year == "2012"),
                            contrasts = list(Wind.Exposure = "contr.sum",
                                             Genotype = "contr.sum",
                                             Year = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
summary(aphid.abund.glmer)
overdisp_fun(aphid.abund.glmer) # no overdispersion
plot(aphid.abund.glmer) 

anova.table(aphid.abund.glmer, test = "Chisq", experiment = "wind") # nothing significant

## Wind: Tortricid leaf tier analysis
with(wind.arth.df, table(Tortricidiae_leaftier, Year)) # no leaftiers in 2013, therefore, I only analyze 2012

leaftier.abund.glmer <- glmer(Tortricidiae_leaftier ~ Wind.Exposure*Year + Wind.Exposure*Genotype + 
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
summary(leaftier.abund.glmer)
overdisp_fun(leaftier.abund.glmer) # no overdispersion
plot(leaftier.abund.glmer) 

leaftier.anova <- anova.table(leaftier.abund.glmer, test = "Chisq", experiment = "wind") # nothing significant

leaftier.up <- update(leaftier.abund.glmer, .~. -Wind.Exposure:Year -Wind.Exposure:Genotype -Wind.Exposure -Genotype +(1|Genotype))

plotREsim(REsim(leaftier.up))  
leaftier.R2 <- var.table(leaftier.up, experiment = "wind")

## Wind: community analysis ----
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
colSums(aa.arth.df[ ,aa.arth.names])
round(colSums(aa.arth.df[ ,aa.arth.names])/sum(colSums(aa.arth.df[ ,aa.arth.names]))*100,0)

# look at response from spiders, leafhoper, leaftier, LTF, ant_black, Aphids_nonAphis and aphid growth rates.

## ant-aphid: arthropod abundance analysis ----

# GLMM
aa.arth.abund.glmer <- glmer(total.abund ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                            (1|plant_code) +
                            (1|Block/fact.Ant.mound.dist),
                          data = aa.arth.df,
                          contrasts = list(Aphid.treatment = "contr.sum",
                                           Genotype = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
print(summary(aa.arth.abund.glmer), correlation = TRUE)
overdisp_fun(aa.arth.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(aa.arth.abund.glmer) 

aa.arth.abund.anova <- anova.table(aa.arth.abund.glmer, test = "Chisq", experiment = "ant-aphid")

aa.arth.abund.up <- update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (Aphid.treatment|Genotype) - (1|plant_code))

plotFEsim(FEsim(aa.arth.abund.up))
with(aa.arth.df, interaction.plot(Ant.mound.dist, Aphid.treatment, total.abund))
with(aa.arth.df, interaction.plot(Aphid.treatment, Genotype, total.abund))
plotREsim(REsim(aa.arth.abund.up))

aa.arth.abund.R2 <- var.table(aa.arth.abund.up, experiment = "ant-aphid")

## ant-aphid: arthropod richness analysis ----
plot(total.rich ~ total.abund, aa.arth.df)

# GLMM
aa.arth.rich.glmer <- glmer(total.rich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                           #(1|plant_code) +
                           (1|Block/fact.Ant.mound.dist),
                         data = aa.arth.df,
                         contrasts = list(Aphid.treatment = "contr.sum",
                                          Genotype = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
print(summary(aa.arth.rich.glmer), correlation = TRUE)
overdisp_fun(aa.arth.rich.glmer) # no overdispersion
plot(aa.arth.rich.glmer) 

aa.arth.rich.anova <- anova.table(aa.arth.rich.glmer, test = "Chisq", experiment = "ant-aphid")

aa.arth.rich.up <- update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + (1|Genotype))

plotREsim(REsim(aa.arth.rich.up))

aa.arth.rich.R2 <- var.table(aa.arth.rich.up, experiment = "ant-aphid")

## ant-aphid: rarefied arthropod richness analysis ----
aa.arth.df$total.rarerich <- rarefy(aa.arth.df[ ,aa.arth.names], 2) - 1

hist(filter(aa.arth.df, total.abund > 2)$total.rarerich)

# GLMM. Note that logit-transformation (common for proportion data) did not qualitatively affect the outcome or the residuals
aa.arth.rarerich.glmer <- lmer(total.rarerich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                              #(1|plant_code) +
                              (1|Block/fact.Ant.mound.dist),
                            data = filter(aa.arth.df, total.abund > 2),
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"))
print(summary(aa.arth.rarerich.glmer), correlation = TRUE)
plot(aa.arth.rarerich.glmer) 

aa.arth.rarerich.anova <- anova.table(aa.arth.rarerich.glmer, test = "F", experiment = "ant-aphid")

aa.arth.rarerich.up <- update(aa.arth.rarerich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype)

aa.arth.rarerich.R2 <- var.table(aa.arth.rarerich.up, experiment = "ant-aphid")

## ant-aphid: herbivore abundance analysis ----
hist(aa.arth.df$herb.abund.nonAphis)

# GLMM
aa.herb.nonAphis.glmer <- glmer(herb.abund.nonAphis ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                              (1|plant_code) +
                              (1|Block/fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.herb.nonAphis.glmer), correlation = TRUE)
overdisp_fun(aa.herb.nonAphis.glmer) # no overdispersion
plot(aa.herb.nonAphis.glmer) 

aa.herb.nonAphis.anova <- anova.table(aa.herb.nonAphis.glmer, test = "Chisq", experiment = "ant-aphid")

aa.herb.nonAphis.up <- update(aa.herb.nonAphis.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (Aphid.treatment|Genotype) -(1|plant_code))

plotREsim(REsim(aa.herb.nonAphis.up))
with(aa.arth.df, interaction.plot(Aphid.treatment, Genotype, herb.abund.nonAphis))
with(aa.arth.df, interaction.plot(Ant.mound.dist, Aphid.treatment, herb.abund.nonAphis))

aa.herb.nonAphis.R2 <- var.table(aa.herb.nonAphis.up, experiment = "ant-aphid")

## ant-aphid: predator abundance analysis ----

# GLMM
aa.pred.nonFobs.glmer <- glmer(pred.abund.nonFobs ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                  #(1|plant_code) +
                                  (1|Block/fact.Ant.mound.dist),
                                data = aa.arth.df,
                                contrasts = list(Aphid.treatment = "contr.sum",
                                                 Genotype = "contr.sum"),
                                family = "poisson",
                                control=glmerControl(optimizer="bobyqa",
                                                     optCtrl=list(maxfun=2e5)))
print(summary(aa.pred.nonFobs.glmer), correlation = TRUE)
overdisp_fun(aa.pred.nonFobs.glmer) # no overdispersion
plot(aa.pred.nonFobs.glmer) 

aa.pred.nonFobs.anova <- anova.table(aa.pred.nonFobs.glmer, test = "Chisq", experiment = "ant-aphid")

aa.pred.nonFobs.up <- update(aa.pred.nonFobs.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + Aphid.treatment + (1|Genotype))

plotFEsim(FEsim(aa.pred.nonFobs.up))
plotREsim(REsim(aa.pred.nonFobs.up))
with(aa.arth.df, plot(pred.abund.nonFobs ~ Aphid.treatment + Genotype))
with(aa.arth.df, interaction.plot(Ant.mound.dist, Aphid.treatment, pred.abund.nonFobs))

aa.pred.nonFobs.R2 <- var.table(aa.pred.nonFobs.up, experiment = "ant-aphid")

## ant-aphid: Formica obscuripes abundance analysis ----

# GLMM
aa.Fobs.abund.glmer <- glmer(ant_F_obscuripes ~ scale(Ant.mound.dist) + Aphid.treatment + (1|Genotype) +
                                 (1|plant_code) +
                                 (1|Block/fact.Ant.mound.dist),
                               data = aa.arth.df,
                               contrasts = list(Aphid.treatment = "contr.sum",
                                                Genotype = "contr.sum"),
                               family = "poisson",
                               control=glmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun=2e5)))
print(summary(aa.Fobs.abund.glmer), correlation = TRUE)
overdisp_fun(aa.Fobs.abund.glmer) # no overdispersion
plot(aa.Fobs.abund.glmer) 

aa.Fobs.abund.anova <- anova.table(aa.Fobs.abund.glmer, test = "Chisq", experiment = "ant-aphid")

anova(aa.Fobs.abund.glmer, update(aa.Fobs.abund.glmer, .~. -(1|Genotype))) # no genotype effect

aa.Fobs.abund.up <- update(aa.Fobs.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + Aphid.treatment -(1|Genotype) -(1|Block/fact.Ant.mound.dist) + (1|Block:fact.Ant.mound.dist) -(1|plant_code)) 

plotFEsim(FEsim(aa.Fobs.abund.up))
plotREsim(REsim(aa.Fobs.abund.up))

with(aa.arth.df, plot(ant_F_obscuripes ~ Aphid.treatment))

aa.Fobs.abund.R2 <- var.table(aa.Fobs.abund.up, experiment = "ant-aphid")

## ant-aphid: Aphis growth rates analysis ----

# calculate average aphid growth rate for each plant.
aa.aphid.GR.sum <- aa.aphid.GR %>%
  group_by(Block, Genotype, Ant.mound.dist, Aphid.treatment, plant_code) %>%
  dplyr::summarise(mean.Aphid.GR = mean(Aphis.growth.rate, na.rm = TRUE)) %>%
  mutate(fact.Ant.mound.dist = as.factor(Ant.mound.dist))
hist(aa.aphid.GR.sum$mean.Aphid.GR) # note though that growth rates are virtually all negative.

aphid.GR.lmer <- lmer(mean.Aphid.GR ~ scale(Ant.mound.dist)*Genotype + (1|Block/fact.Ant.mound.dist),
                      aa.aphid.GR.sum,
                      contrasts = list(Genotype = "contr.sum"))
summary(aphid.GR.lmer)
plot(aphid.GR.lmer)

aphid.GR.anova <- anova.table(aphid.GR.lmer, test = "F", type = 3, experiment = "ant-aphid") # effect of genotype on aphid growth rate.

aphid.GR.up <- update(aphid.GR.lmer, .~. -scale(Ant.mound.dist)*Genotype + (1|Genotype))

aphid.GR.R2 <- var.table(aphid.GR.up, experiment = "ant-aphid")



# only 2012 ----
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

## Ant-aphid: community analysis

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


