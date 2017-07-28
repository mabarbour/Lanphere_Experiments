############################################
## Description: This script analyzes arthropod community data for the Lanphere experiments.
## Code author: Matt Barbour
## Email: barbour@zoology.ubc.ca
############################################

## load required libraries ----
#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(psych)
library(lme4)
library(car)
library(vegan)
library(effects)
library(ggplot2)
library(cowplot)

library(RCurl) # for loading github files directly.

script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/autoplot.custom.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/veganCovEllipse.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

## upload datasets ----

## ant-aphid: aphid growth rates
aa.aphid.GR <- read.csv('final_data/ant_aphid_Aphis_popgrowth_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  select(Date_rel:plant_ID, Aphis.growth.rate)
glimpse(aa.aphid.GR)

## ant-aphid: arthropod community
aa.arth.df <- read.csv('final_data/ant_aphid_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(#Ants_all = ant_F_obscuripes + ant_black,
         #Aphids_nonAphis = aphid_Tuberolachnus + aphid_LG,
         X = as.factor(X),
         Block = as.factor(Block),
         fact.Ant.mound.dist = as.factor(Ant.mound.dist),
         ord.Ant.mound.dist = ordered(Ant.mound.dist),
         GxE = C(interaction(Genotype,Aphid.treatment), "contr.sum")) #%>% #,
         #Plot_code = paste(Block, fact.Ant.mound.dist, sep = "_")) %>%

aa.arth.names <- colnames(select(aa.arth.df, Gracilliaridae_miner:Spider))
 # select(Block, Genotype, Ant.mound.dist, Aphid.treatment,
  #       plant_ID, fact.Ant.mound.dist, Plot_code,
   #      aphid_Aphis, ant_F_obscuripes, ant_black, 
    #     Aphids_nonAphis,
     #    psyllid:sawfly_larva, syrphid_larva:LTF_Caloptilia)

#aa.arth.names <- colnames(select(aa.arth.df,
 #                                ant_black:LTF_Caloptilia)) # excluding Aphis farinosa and Formica obscuripes

# generate new columns for total abundance and richness
#aa.total.abund <- rowSums(
 # select(aa.arth.df, ant_black:LTF_Caloptilia))
#aa.total.rich <- rowSums(
 # select(aa.arth.df, ant_black:LTF_Caloptilia) > 0)
#aa.arth.df$total.rarerich <- rarefy(aa.arth.df[ ,aa.arth.names], 2) - 1
#aa.na.rarerich <- which(aa.arth.df$total.abund < 2) # can only test rarefied richness when at least 2 individuals were sampled.
#aa.arth.df$total.rarerich[aa.na.rarerich] <- NA
#aa.herb.abund <- rowSums(
#  select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,gall_R_salicisbattatus,grasshopper,leaftier_Tortricid,LTF_Caloptilia)
#)
#aa.herb.rich <- rowSums(
 # select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,gall_R_salicisbattatus,grasshopper,leaftier_Tortricid,LTF_Caloptilia) > 0
#)
#aa.herb.conceal <- rowSums(
 # select(aa.arth.df, gall_R_salicisbattatus,leaftier_Tortricid,LTF_Caloptilia)
#)
#aa.herb.nonconceal <- rowSums(
 # select(aa.arth.df, Aphids_nonAphis,psyllid,leafhopper,froghopper,sawfly_larva,grasshopper)
#)
#aa.pred.abund <- rowSums(
 # select(aa.arth.df, ant_black, spiders, syrphid_larva)
#)
#aa.total.pred.abund <- rowSums(
 # select(aa.arth.df, ant_F_obscuripes, ant_black, spiders, syrphid_larva)
#)
#aa.pred.rich <- rowSums(
 # select(aa.arth.df, ant_black, spiders, syrphid_larva) > 0
#)
#aa.arth.df <- mutate(aa.arth.df, 
 #                    total.abund = aa.total.abund,
  #                   total.rich = aa.total.rich)#,
                    # herb.abund.nonAphis = aa.herb.abund,
                     #herb.rich.nonAphis = aa.herb.rich,
                     #herb.abund.conceal = aa.herb.conceal,
                     #herb.abund.nonconceal = aa.herb.nonconceal,
                     #pred.abund.nonFobs = aa.pred.abund,
                     #pred.rich.nonFobs = aa.pred.rich,
                     #pred.abund.all = aa.total.pred.abund)
#aa.arth.df$GxE <- with(aa.arth.df, C(interaction(Genotype,Aphid.treatment), "contr.sum"))

# upload ant-aphid trait data
#aa.trait.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/ant_aphid_trait_df.csv') %>% 
 # filter(Year == "2012") %>% 
  #mutate(plant_ID = paste(Block, Plant_Position, sep = "_")) %>%
  #select(plant_ID, Height, all.shoot.count, mature.shoot.avg.length, leaf_WC, leaf_trichome.density)

#aa.arth.df <- left_join(aa.arth.df, aa.trait.df, by = "plant_ID")

# subset of data where plants had at least one arthropod individual
aa.arth.12.pos <- aa.arth.df %>%
  filter(total.abund > 1) 


## wind: arthropod community
# dead plants have already been removed
wind.arth.df <- read.csv('final_data/wind_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(X = as.factor(X),
         Block = as.factor(Block),
         Year = as.factor(Year),
         Plot_code = interaction(Block, Wind.Exposure),
         GxE = C(interaction(Genotype, Wind.Exposure), contr = "contr.sum", how.many = 9))  # create a new variable for the interaction to permit testing of main effects with type 3 sum of squares. 

## upload trait data from wind experiment for determining direct and indirect effects
#w.trait.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_trait_df.csv') %>%
 # tbl_df() %>%
  #mutate(plant_ID = paste(treatment, block, plant.position, sep = "_"),
   #      Year = as.factor(Year)) %>%
  #select(Year, plant_ID, Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC, SLA, leaf_C_N) 
#glimpse(w.trait.df)

#wind.arth.df <- left_join(wind.arth.df, w.trait.df, by = c("Year","plant_ID"))

wind.arth.names <- colnames(select(wind.arth.df, Gracilliaridae_miner:Spider)) # for subsetting community data

# generate new columns for total abundance and richness
#wind.arth.df$total.abund <- rowSums(
 # select(wind.arth.df, Gracilliaridae_miner:Spider))
#wind.arth.df$total.rich <- rowSums(
 # select(wind.arth.df, Gracilliaridae_miner:Spider) > 0)
#herb.abund <- rowSums(
 # select(wind.arth.df, Gracilliaridae_miner:Aphididae))
#herb.rich <- rowSums(
 # select(wind.arth.df, Gracilliaridae_miner:Aphididae) > 0)
#pred.abund <- rowSums(
 # select(wind.arth.df, Formica_ant, Spider))
#pred.rich <- rowSums(
 # select(wind.arth.df, Formica_ant, Spider) > 0)
#wind.arth.df$total.rarerich <- rarefy(wind.arth.df[ ,wind.arth.names], 2) - 1
#na.rarerich <- which(w.arth.df$total.abund < 2) # can only test rarefied richness when at least 2 individuals were sampled.
#w.arth.df$total.rarerich[na.rarerich] <- NA
#wind.arth.df <- mutate(wind.arth.df, 
      #                 total.abund = total.abund,
      #                 total.rich = total.rich)#,
                       #herb.abund = herb.abund,
                       #herb.rich = herb.rich,
                       #pred.abund = pred.abund,
                       #pred.rich = pred.rich,
             #          X = factor(seq(1,362,1)))

# 2012 dataset
# focus dataset on aggregated Family/Order arthropod groupings
w.arth.12.full <- wind.arth.df %>%
  filter(Year == "2012") %>%
  select(Block:plant_ID, Plot_code, GxE, total.abund,
         Gracilliaridae_miner:Spider) 

# subset of data where plants had at least one arthropod individual
w.arth.12.pos2 <- w.arth.12.full %>%
  filter(total.abund > 1) 

# for avoiding pseudoreplication while testing for wind effect. Take average abundance of each arthropod at the plot level.
#w.effect.12 <- w.arth.12.full %>%
 # select(-plant_ID, -Genotype) %>%
  #group_by(Block, Wind.Exposure, Plot_code) %>%
  #summarise_each(funs(mean))

# 2013 dataset
# same structure as 2012 dataset
w.arth.13.full <- wind.arth.df %>%
  filter(Year == "2013") %>%
  select(Block:plant_ID, Plot_code, GxE, total.abund,
         Gracilliaridae_miner:Spider) 

w.arth.13.pos2 <- w.arth.13.full %>%
  filter(total.abund > 1)

#w.effect.13 <- w.arth.13.full %>%
 # select(-plant_ID, -Genotype) %>%
  #group_by(Block, Wind.Exposure, Plot_code) %>%
  #summarise_each(funs(mean))

## Wind: arthropod abundance analysis ----

# GLMM
arth.abund.glmer <- glmer(total.abund ~ Wind.Exposure*Genotype*Year +
                             (1|Block) + 
                             (1|Block:Wind.Exposure) +
                             (1|X) +
                             (1|plant_ID),
                           data = wind.arth.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum",
                                            Year = "contr.sum"),
                           family = "poisson",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))
print(summary(arth.abund.glmer), correlation = TRUE)
overdisp_fun(arth.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(arth.abund.glmer) 

library(sjPlot)
sjp.glmer(update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Genotype), type = "fe")

# Likelihood ratio tests
(w.abund.3 <- drop1(arth.abund.glmer, test = "Chisq"))
(w.abund.2 <- drop1(update(arth.abund.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq"))
(w.abund.1 <- drop1(update(arth.abund.glmer, .~. -Wind.Exposure:Genotype:Year -Wind.Exposure:Genotype -Wind.Exposure:Year -Genotype:Year), test = "Chisq"))

#arth.abund.anova <- bind_rows(w.abund.1, w.abund.2, w.abund.3) %>%
 # mutate(Factor = c(rownames(w.abund.1), rownames(w.abund.2), rownames(w.abund.3))) %>%
#  filter(Df != "NA") %>%
 # select(Factor, Df, LRT, P = Pr(Chi))

#(arth.abund.anova <- anova.table(arth.abund.glmer, test = "Chisq", experiment = "wind"))
#plot(Effect("Height", arth.abund.glmer))

## Effects
Effect(c("Year"), arth.abund.glmer) # 1.4-fold more arthropods in 2013 vs. 2012. 
Effect(c("Wind.Exposure"), arth.abund.glmer) # 1.9 more arthropods on unexposed plants; 47% fewer arthropods on wind-exposed plants.
Effect(c("Genotype"), arth.abund.glmer) # arthropod abundance varied 4.7-fold among genotypes

w.abund.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.abund.glmer)) %>% mutate(response = "arthropod_abund")

## Calculate R2 for significant predictors
arth.abund.up <- update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

(arth.abund.R2 <- var.table(arth.abund.up, experiment = "wind"))

## Wind: arthropod richness analysis ----
plot(log(total.rich+1) ~ log(total.abund+1), wind.arth.df)

# GLMM
arth.rich.glmer <- glmer(total.rich ~ Wind.Exposure*Year*Genotype +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_ID),
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

arth.rich.glmer.P <- glmer(total.rich ~ Wind.Exposure + Genotype + Year + (1|Block) +  (1|Block:Wind.Exposure) + (1|plant_ID), data = wind.arth.df, family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
sjp.glmer(arth.rich.glmer.P, type = "fe")

# Likelihood ratio tests
(w.rich.3 <- drop1(arth.rich.glmer, test = "Chisq"))
(w.rich.2 <- drop1(update(arth.rich.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq"))
(w.rich.1 <- drop1(update(arth.rich.glmer, .~. -Wind.Exposure:Genotype:Year -Wind.Exposure:Genotype -Wind.Exposure:Year -Genotype:Year), test = "Chisq"))
#(arth.rich.anova <- anova.table(arth.rich.glmer, test = "Chisq", experiment = "wind"))

## Effects
Effect("Year", arth.rich.glmer) # 1.4-fold more arthropod species in 2013 vs. 2012
Effect("Wind.Exposure", arth.rich.glmer) # arthropod richness on exposed plants was 51% less.
Effect("Genotype", arth.rich.glmer) # arthropod richness varied 3.1-fold among willow genotypes

w.rich.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.rich.glmer)) %>% mutate(response = "arthropod_rich")

## Calculate R2
arth.rich.up <- update(arth.rich.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

(arth.rich.R2 <- var.table(arth.rich.up, experiment = "wind"))

## Wind: rarefied arthropod richness analysis ----
# need to filter data so that there are 2 or more individuals
hist(filter(wind.arth.df, total.abund > 1)$total.rarerich)

# GLMM
arth.rarerich.glmer <- lmer(total.rarerich ~ Wind.Exposure*Year+ Wind.Exposure*Genotype +
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           #(1|X) +
                           (1|plant_ID),
                         data = filter(wind.arth.df, total.abund > 1),
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"))
print(summary(arth.rarerich.glmer), correlation = TRUE)
#overdisp_fun(arth.rarerich.glmer) 
plot(arth.rarerich.glmer) 

## Likelihood ratio tests
car::Anova(arth.rarerich.glmer, test = "F", type = 2)
#(w.rrich.2 <- drop1(arth.rarerich.glmer, test = "Chisq"))
#(w.rrich.1 <- drop1(update(arth.rarerich.glmer, .~. -Wind.Exposure:Year -Wind.Exposure:Genotype), test = "Chisq")) # marginal genotype effect
#arth.rarerich.anova <- anova.table(arth.rarerich.glmer, test = "F", experiment = "wind")

## Calculate effects
Effect("Wind.Exposure", arth.rarerich.glmer) # rarefied richness of arthropods was 60% less on exposed plants
Effect("Genotype", arth.rarerich.glmer) # 5.9-fold variation in arthropod rarefied richness among genotypes

w.rrich.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.rarerich.glmer)) %>% mutate(response = "arthropod_rarerich")

# Calculate R2
arth.rarerich.up <- update(arth.rarerich.glmer, .~. -Wind.Exposure*Year -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))

(arth.rarerich.R2 <- var.table(arth.rarerich.up, experiment = "wind"))

## Wind: hellinger distance 2012 ----
# no G, wind, or GxE effect
w.12.hell.pos2 <- decostand(w.arth.12.pos2[ ,wind.arth.names], method = "hellinger")

# test G and GxE
w.12.hell.rda <- rda(w.12.hell.pos2 ~ Wind.Exposure*Genotype, data = w.arth.12.pos2) # no G or GxE
anova(w.12.hell.rda, by = "margin", permutations = how(block = w.arth.12.pos2$Block, nperm = 999)) # no GxE
anova(update(w.12.hell.rda, .~. -Wind.Exposure:Genotype), by = "margin", permutations = how(block = w.arth.12.pos2$Block, nperm = 999)) # no G

# test wind effect
w.12.plots.hell <- betadisper(vegdist(w.12.hell.pos2, method = "euclidean"),  w.arth.12.pos2$Plot_code, bias.adjust = TRUE)

w.12.plots.centr.hell <- data.frame(w.12.plots.hell$centroids, 
                               id = rownames(w.12.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.12.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.12.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.12.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Wind: Hellinger distance 2013 ----
# effect of wind only, but could be due to overdispersion
key.arths <- which(wind.arth.names %in% c("Tortricidiae_leaftier","Cecidomyiidae_gall"))
w.13.hell.pos2 <- decostand(w.arth.13.pos2[ ,wind.arth.names], method = "hellinger")

# test G and GxE
w.13.hell.rda <- rda(w.13.hell.pos2 ~ Wind.Exposure*Genotype, data = w.arth.13.pos2) 
plot(w.13.hell.rda, display = "sp") # key species are Cecidomyiids and Tortricids

anova(w.13.hell.rda, by = "margin",permutations = how(block = w.arth.13.pos2$Block, nperm = 999)) # no GxE
anova(update(w.13.hell.rda, .~. -Wind.Exposure:Genotype), by = "margin",permutations = how(block = w.arth.13.pos2$Block, nperm = 999)) # no G

# test wind effect
w.13.plots.hell <- betadisper(vegdist(w.13.hell.pos2, method = "euclidean"),  w.arth.13.pos2$Plot_code, bias.adjust = TRUE)

w.13.plots.centr.hell <- data.frame(w.13.plots.hell$centroids, 
                               id = rownames(w.13.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.13.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.13.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.13.plots.centr.hell$Block, nperm = 999)) # significant effect of wind exposure

w.hell.13.wind <- betadisper(vegdist(w.13.plots.hell$centroids, method = "euclidean"), w.13.plots.centr.hell$Wind.Exposure, bias.adjust = TRUE)
boxplot(w.hell.13.wind)
plot(w.hell.13.wind)
permutest(w.hell.13.wind, permutations = how(block = w.13.plots.centr.hell$Block, nperm = 999)) # significant effect, suggesting that the effect of wind on centroid location could be due to overdispersion.


w.hell.13.wind.rda <- rda(w.13.hell.pos2 ~ Block*Wind.Exposure, 
                          data = w.arth.13.pos2)
summary(w.hell.13.wind.rda)
plot(w.hell.13.wind.rda, display = c("bp","sp"))
ordiellipse(w.hell.13.wind.rda, groups = interaction(w.arth.13.pos2$Block, w.arth.13.pos2$Wind.Exposure))

w.hell.13.wind.nmds <- metaMDS(vegdist(w.13.plots.hell$centroids, method = "euclidean"))
plot(w.hell.13.wind.nmds)
ordiellipse(w.hell.13.wind.nmds, groups = w.13.plots.centr.hell$Wind.Exposure)
points(w.hell.13.wind.nmds, col = as.numeric(w.13.plots.centr.hell$Wind.Exposure))
text(w.hell.13.wind.nmds, labels = w.13.plots.centr.hell$Block)

## Wind: Relative abundance of Tortricidae ----

leaftier.rel.glmer <- glmer(Tortricidiae_leaftier/total.abund ~ 
                                Wind.Exposure + Genotype + 
                                (1|Block) + 
                                (1|Block:Wind.Exposure) +
                                (1|plant_ID),
                              weights = total.abund,
                              data = filter(wind.arth.df, Year == "2013"),
                              contrasts = list(Wind.Exposure = "contr.sum",
                                               Genotype = "contr.sum",
                                               Year = "contr.sum"),
                              family = "binomial",
                              control=glmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))
summary(leaftier.rel.glmer)
overdisp_fun(leaftier.rel.glmer) # no overdispersion
plot(leaftier.rel.glmer) 
drop1(leaftier.rel.glmer, test = "Chisq") # confirms large effect of wind exposure on relative abundance
plot(Effect("Wind.Exposure", leaftier.rel.glmer))

## Wind: Relative abundance of Cecidomyiids ----
gall.rel.glmer <- glmer(Cecidomyiidae_gall/total.abund ~ Wind.Exposure+Genotype +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|plant_ID),
                          weights = total.abund,
                          data = filter(wind.arth.df, Year == "2013"),
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum"),
                          family = "binomial",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(gall.rel.glmer)
overdisp_fun(gall.rel.glmer) # no overdispersion
plot(gall.rel.glmer) 
drop1(gall.rel.glmer, test = "Chisq") # confirms large effect of wind exposure on relative abundance, also points out that Genotypes vary in their relative abundance too.
plot(Effect("Wind.Exposure", gall.rel.glmer))

## Wind: Caloptilia analysis ----
hist(wind.arth.df$Gracilliaridae_miner)
with(wind.arth.df, table(Gracilliaridae_miner, Year))

# GLMM
LTF.abund.glmer <- glmer(Gracilliaridae_miner ~ Wind.Exposure*Year + 
                           Genotype + # only model without convergence
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           #(1|X) +
                           (1|plant_ID),
                         data = wind.arth.df,
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(LTF.abund.glmer)
overdisp_fun(LTF.abund.glmer) 
plot(LTF.abund.glmer) 

## Likelihood-ratio test
(LTF.2 <- drop1(LTF.abund.glmer, test = "Chisq")) # G effect
(LTF.2 <- drop1(update(LTF.abund.glmer, .~. -Wind.Exposure:Year), test = "Chisq")) # G and E effect, marginal Year effect
#(LTF.abund.anova <- anova.table(LTF.abund.glmer, test = "Chisq", experiment = "wind"))


## calculate variance components
LTF.abund.up <- update(LTF.abund.glmer, .~. -Wind.Exposure:Year
                       -Genotype + (1|Genotype)) 
(LTF.abund.R2 <- var.table(LTF.abund.up, experiment = "wind"))

## Functional trait model
# 2012
LTF.abund.2012 <- glmer(Gracilliaridae_miner ~ scale(Height) + scale(all.shoot.count) + scale(leaf_WC) + scale(leaf_trichome.density) +
                          (1|Genotype) +
                          (1|Block) + 
                          (1|Block:Wind.Exposure) +
                          #(1|X) +
                          (1|plant_ID),
                        data = filter(wind.arth.df[-nonsense, ], Year == "2012"),
                        family = "poisson",
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e5)))
summary(LTF.abund.2012)
overdisp_fun(LTF.abund.2012)
anova(LTF.abund.2012, update(LTF.abund.2012, .~. -(1|Genotype))) # genotype still has a significant effect.
plot(Effect("leaf_WC", LTF.abund.2012))
plot(Effect("leaf_C_N", LTF.abund.2012))

# 2013
LTF.abund.2013 <- glmer(Gracilliaridae_miner ~ scale(Height) + scale(leaf_WC) + scale(leaf_C_N) +
                          (1|Genotype) +
                          (1|Block) + 
                          (1|Block:Wind.Exposure),# +
                        #(1|X) +
                        #(1|plant_ID),
                        data = filter(wind.arth.df[-nonsense, ], Year == "2013"),
                        family = "poisson",
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e5)))
summary(LTF.abund.2013)
overdisp_fun(LTF.abund.2013)
anova(LTF.abund.2013, update(LTF.abund.2013, .~. -(1|Genotype))) # genotype still has a significant effect.
plot(Effect("leaf_WC", LTF.abund.2013))
plot(Effect("leaf_C_N", LTF.abund.2013))

## Wind: Tortricidae leaftier analysis ----
hist(wind.arth.df$Tortricidiae_leaftier)
with(wind.arth.df, table(Tortricidiae_leaftier, Year))

# GLMM
leaftier.abund.glmer <- glmer(Tortricidiae_leaftier/total.abund ~ 
                                #Wind.Exposure*Year + 
                                Wind.Exposure + Genotype + # only model with convergence
                                (1|Block) + 
                                (1|Block:Wind.Exposure) +
                                #(1|X) +
                                (1|plant_ID),
                              weights = total.abund,
                              data = filter(wind.arth.df, Year == "2013"),
                              contrasts = list(Wind.Exposure = "contr.sum",
                                               Genotype = "contr.sum",
                                               Year = "contr.sum"),
                              family = "binomial",#"poisson",
                              control=glmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))
summary(leaftier.abund.glmer)
overdisp_fun(leaftier.abund.glmer) # no overdispersion
plot(leaftier.abund.glmer) 

## Likelihood ratio test
(lt.2 <- drop1(leaftier.abund.glmer, test = "Chisq"))
(lt.1 <- drop1(update(leaftier.abund.glmer, .~. -Wind.Exposure:Year -Wind.Exposure:Genotype), test = "Chisq"))
#(leaftier.abund.anova <- anova.table(leaftier.abund.glmer, test = "Chisq", experiment = "wind"))

g.tort <- as.data.frame(Effect("Genotype", leaftier.abund.glmer)) %>%
  select(Genotype, tort.fit = fit)

## calculate variance components
leaftier.abund.up <- update(leaftier.abund.glmer, .~. -Wind.Exposure*Year
                            -Wind.Exposure*Genotype + Year + (1|Genotype))
(leaftier.abund.R2 <- var.table(leaftier.abund.up, experiment = "wind"))

## Wind: Cecidomyiidae galler analysis ----
with(wind.arth.df, table(Cecidomyiidae_gall, Year))

# GLMM
gall.abund.glmer <- glmer(Cecidomyiidae_gall/total.abund ~ Wind.Exposure+Genotype + #Year +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_ID),
                          weights = total.abund,
                          data = filter(wind.arth.df, Year == "2013"),
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "binomial", #"poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(gall.abund.glmer)
overdisp_fun(gall.abund.glmer) # no overdispersion
plot(gall.abund.glmer) 

# likelihood ratio test
(gall.1 <- drop1(gall.abund.glmer, test = "Chisq"))
#(gall.abund.anova <- anova.table(gall.abund.glmer, test = "Chisq", experiment = "wind"))

Effect("Wind.Exposure", gall.abund.glmer)
g.gall <- as.data.frame(Effect("Genotype", gall.abund.glmer)) %>%
  select(Genotype, Cecid_gall.fit = fit)

gall.abund.up <- update(gall.abund.glmer, .~. -Genotype + (1|Genotype))

(gall.abund.R2 <- var.table(gall.abund.up, experiment = "wind"))

## Functional trait models
# 2013
gall.abund.2013 <- glmer(Cecidomyiidae_gall ~ scale(Height) + scale(leaf_WC) + scale(leaf_C_N) +
                           (1|Genotype) +
                           (1|Block) + 
                           (1|Block:Wind.Exposure),# +
                         #(1|X) +
                         #(1|plant_ID),
                         data = filter(wind.arth.df[-nonsense, ], Year == "2013"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(gall.abund.2013)
overdisp_fun(gall.abund.2013)
anova(gall.abund.2013, update(gall.abund.2013, .~. -(1|Genotype))) # genotype no longer a significant predictor.
plot(Cecidomyiidae_gall ~ Height, data = filter(wind.arth.df[-nonsense, ], Year == "2013"))
plot(Cecidomyiidae_gall ~ leaf_WC, data = filter(wind.arth.df[-nonsense, ], Year == "2013")) # don't really trust this...
plot(Cecidomyiidae_gall ~ leaf_C_N, data = filter(wind.arth.df[-nonsense, ], Year == "2013")) # don't really trust this...

## Wind: Spider analysis ----
spider.abund.glmer <- glmer(Spider ~ Wind.Exposure*Year + 
                              Genotype +
                              (1|Block) + 
                              (1|Block:Wind.Exposure) +
                              #(1|X) +
                              (1|plant_ID),
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

## likelihood ratio tests
(sp.2 <- drop1(spider.abund.glmer, test = "Chisq"))
(sp.1 <- drop1(update(spider.abund.glmer, .~. -Wind.Exposure:Year), test = "Chisq"))
#spider.abund.anova <- anova.table(spider.abund.glmer, test = "Chisq", experiment = "wind")

Effect("Wind.Exposure", spider.abund.glmer)

g.spid <- as.data.frame(Effect("Genotype", spider.abund.glmer)) %>%
  select(Genotype, spid.fit = fit)

spider.abund.up <- update(spider.abund.glmer, .~. -Wind.Exposure:Year -Genotype)

(spider.abund.R2 <- var.table(spider.abund.up, experiment = "wind"))

## Wind: Aphididae analysis ----
with(wind.arth.df, table(Aphididae, Year)) # no aphids in 2013, therefore, I only analyze 2012

aphid.abund.glmer <- glmer(Aphididae ~ Wind.Exposure + Genotype +
                             (1|Block) + 
                             (1|Block:Wind.Exposure) +
                             #(1|X) +
                             (1|plant_ID),
                           data = filter(wind.arth.df, Year == "2012"),
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"),
                           family = "poisson",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))
summary(aphid.abund.glmer)
overdisp_fun(aphid.abund.glmer) # sig. overdisper, model with individual-level random effect
plot(aphid.abund.glmer) 

# likelihood ratio test
(aphid.1 <- drop1(aphid.abund.glmer, test = "Chisq")) # no wind exposure or genotype effect
#(aphid.abund.anova <- anova.table(aphid.abund.glmer, test = "Chisq", experiment = "wind")) # nothing significant

g.Aphids_nonAphis <- as.data.frame(Effect("Genotype", aphid.abund.glmer)) %>%
  select(Genotype, Aphids_nonAphis.fit = fit)

## Wind: plots ----
pd <- 0.2
ebar.w <- 0.05
l.size <- 1.5
alp <- 0.5
p.size <- 5
cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")

rich.Gord <- as.data.frame(Effect(c("Genotype"), arth.rich.glmer))  %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)]))

## richness
w.rE <- as.data.frame(Effect(c("Wind.Exposure"), arth.rich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Arthropod richness", limits = c(0,2), breaks = c(0,1,2)) + xlab(""); w.rE

w.rG <- as.data.frame(Effect(c("Genotype"), arth.rich.glmer))  %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0,2), breaks = c(0,1,2)) + xlab(""); w.rG

## abundance
w.aE <- as.data.frame(Effect(c("Wind.Exposure"), arth.abund.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Arthropod abundance", limits = c(0,3.1), breaks = c(0,1,2,3))+ xlab(""); w.aE

w.aG <- as.data.frame(Effect(c("Genotype"), arth.abund.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.Gord$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0, 3.1), breaks = c(0,1,2,3))+ xlab(""); w.aG

## rarefied richness
w.rrE <- as.data.frame(Effect(c("Wind.Exposure"), arth.rarerich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Rarefied richness", limits = c(0,0.75), breaks = c(0,0.25,0.5,0.75)) + xlab("Wind treatment"); w.rrE

w.rrG <- as.data.frame(Effect(c("Genotype"), arth.rarerich.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.Gord$fit)])) %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0,0.75), breaks = c(0,0.25,0.5,0.75)) + xlab("Willow genotype"); w.rrG # set lower confidence intrval for Genotype G to stop at 0.

w.arth.p <- plot_grid(w.rE, w.rG, w.aE, w.aG, w.rrE, w.rrG, ncol = 2, labels = "AUTO", align = 'hv')

save_plot(w.arth.p, filename = "fig_4_wind_arth_comm.png", base_height = 11, base_width = 8.5)

## composition
cp.keyE <- c( "#999999", "#E69F00","#FFFFFF", "#009E73")
w.key.E <- data.frame(Effect(c("Wind.Exposure"), LTF.abund.glmer), response = "leaf-mining moths", sig = "y") %>% 
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), leaftier.abund.glmer), response = "leaf-tiering moths", sig = "n")) %>%
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), gall.abund.glmer), response = "galling midges", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), spider.abund.glmer), response = "spiders", sig = "y")) %>%
  #bind_rows(., data.frame(Effect(c("Wind.Exposure"), aphid.abund.glmer), response = "aphids")) %>% 
  ggplot(aes(x = Wind.Exposure, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Wind treatment") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyE) + theme(legend.justification = c(1,0), legend.position = c(1,0)); w.key.E #+ geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w, position = position_dodge(width = pd)) # na.rm = TRUE, position = position_dodge(width = pd), 

cp.keyG <- c( "#999999", "#E69F00", "#56B4E9", "#FFFFFF") 
w.key.G <- data.frame(Effect(c("Genotype"), LTF.abund.glmer), response = "leaf-mining moths", sig = "y") %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), leaftier.abund.glmer), response = "leaf-tiering moths", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), gall.abund.glmer), response = "galling midges", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), spider.abund.glmer), response = "spiders", sig = "n")) %>%
  #bind_rows(., data.frame(Effect(c("Genotype"), aphid.abund.glmer), response = "aphids")) %>% 
  ggplot(aes(x = Genotype, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Willow genotype") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyG) + theme(legend.justification = c(1,0), legend.position = c(1,0)); w.key.G

#w.arth.p <- plot_grid(w.rE, w.rG, w.aE, w.aG, w.rrE, w.rrG, w.key.E, w.key.G, labels = "AUTO", ncol = 2, align = 'hv'); w.arth.p

#save_plot(w.arth.p, filename = "fig_wind_arth_test.png", base_height = 11, base_width = 8.5)

## manuscript plots
w.rGE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.rich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit, color = Genotype)) + geom_line(aes(group = Genotype), size = l.size) + geom_point(size = p.size)  + scale_color_manual(values = cbPal.10) + ylab("No. of species") + xlab("Wind exposure")+ theme(legend.position = "none") + ggtitle("Arthropods")

# Wind effect on community composition
levels(w.13.plots.centr.hell$Wind.Exposure) <- c("Exp.","Unexp.")
w.hell.wind.rda <- rda(w.13.plots.hell$centroids ~ Block + Wind.Exposure, data = w.13.plots.centr.hell) 
summary(w.hell.wind.rda)$cont$importance[ ,1:2] # first axis explains 46% of the variation, second axis explains 7%

ellip <- ordiellipse(w.hell.wind.rda, groups = w.13.plots.centr.hell$Wind.Exposure, draw = "polygon")

centroids.cap <- data.frame(scores(w.hell.wind.rda, display = "cn"), Wind.Exposure = levels(w.13.plots.centr.hell$Wind.Exposure))

sites.cap <- data.frame(scores(w.hell.wind.rda, choices = c(1,2), display = "sites"), droplevels(w.13.plots.centr.hell$Wind.Exposure), droplevels(w.13.plots.centr.hell$Block))
colnames(sites.cap)[3] <- "Wind.Exposure"
colnames(sites.cap)[4] <- "Block"

source('~/Documents/miscellaneous_R/veganCovEllipse.R')

# get data for ellipse. 
df_ell.cap <- data.frame() #sets up a data frame before running the function.
for(g in levels(sites.cap$Wind.Exposure)){
  df_ell.cap <- rbind(df_ell.cap, 
                      cbind(as.data.frame(
                        with(sites.cap[sites.cap$Wind.Exposure == g, ], 
                             veganCovEllipse(ellip[[g]]$cov, ellip[[g]]$center, ellip[[g]]$scale))), Wind.Exposure = g))
}

w.arth.ord.wind <- ggplot(data = df_ell.cap, aes(x = RDA1, y = RDA2, group = Wind.Exposure)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_text(data = sites.cap, aes(x = RDA1, y = RDA2, label = Block), color = "gray") +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap[11:12, ], 
            aes(x = RDA1, y = RDA2, label = Wind.Exposure), size = 5) +
  ylab("RDA 2 (7%)") + xlab("RDA 1 (46%)") + ggtitle("Arthropods")
w.arth.ord.wind

## ant-aphid: maximum Aphis abundance analysis ----

# GLMM
aa.Aphis.glmer <- glmer(aphid_Aphis ~ scale(Ant.mound.dist)*Genotype +
                          (1|plant_ID) + 
                          (1|Block/fact.Ant.mound.dist),
                        data = filter(aa.arth.df, Aphid.treatment == "aphid"),
                        contrasts = list(Genotype = "contr.sum"),
                        family = "poisson",
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e5)))
overdisp_fun(aa.Aphis.glmer)
summary(aa.Aphis.glmer)

## Likelihood ratio tests
(aa.Aphis.2 <- drop1(aa.Aphis.glmer, test = "Chisq"))
(aa.Aphis.2 <- drop1(update(aa.Aphis.glmer, .~. -scale(Ant.mound.dist):Genotype), test = "Chisq"))
#(aa.Aphis.anova <- anova.table(aa.Aphis.glmer, test = "Chisq", experiment = "ant-aphid")) 

## effects
plot(Effect("Ant.mound.dist", aa.Aphis.glmer))
plot(Effect("Genotype", aa.Aphis.glmer)) # genotypes varied 153-fold in mean aphid abundance.
Effect("Genotype", aa.Aphis.glmer)

aa.Aphis.G <- as.data.frame(Effect("Genotype", aa.Aphis.glmer)) %>% mutate(response = "Aphis_farinosa")

## Calculate R2
aa.Aphis.up <- update(aa.Aphis.glmer, .~. -scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.Aphis.R2 <- var.table(aa.Aphis.up, experiment = "ant-aphid"))

## ant-aphid: Formica obscuripes abundance analysis ----

# GLMM on abundance
aa.Fobs.abund.glmer <- glmer(ant_F_obscuripes ~ Aphid.treatment + scale(Ant.mound.dist) + (Aphid.treatment|Genotype) +
                               #(1|plant_ID) +
                               (1|Block/fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "poisson",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
summary(aa.Fobs.abund.glmer)
plot(aa.Fobs.abund.glmer) 
overdisp_fun(aa.Fobs.abund.glmer)

#plot(Effect(c("Aphid.treatment","Genotype"), aa.Fobs.abund.glmer))
#plot(Effect(c("aphid_Aphis","Ant.mound.dist"), aa.Fobs.abund.glmer))

(aa.Fobs.1 <- drop1(aa.Fobs.abund.glmer, test = "Chisq"))
anova(aa.Fobs.abund.glmer, update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype) + (1|Genotype))) # sig GxE
anova(update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype) + (1|Genotype)), update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype))) # test G effect

var.table(update(aa.Fobs.abund.glmer, .~. +(1|plant_ID)), "ant-aphid") # need to add individual-level random effect when calculating R2

## ant-aphid: arthropod abundance analysis ----

# GLMM
aa.arth.abund.glmer <- glmer(total.abund ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                            (1|plant_ID) +
                            (1|Block/fact.Ant.mound.dist),
                          data = aa.arth.df,
                          contrasts = list(Aphid.treatment = "contr.sum",
                                           #fact.Ant.mound.dist = "contr.sum",
                                           Genotype = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
print(summary(aa.arth.abund.glmer), correlation = TRUE)
overdisp_fun(aa.arth.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(aa.arth.abund.glmer) 

## Likelihood ratio tests
(aa.abund.3 <- drop1(aa.arth.abund.glmer, test = "Chisq"))
(aa.abund.2 <- drop1(update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aa.abund.1 <- drop1(update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq"))
#aa.arth.abund.anova <- anova.table(aa.arth.abund.glmer, test = "Chisq", experiment = "ant-aphid")

## calculate effects
#plot(Effect(c("Aphid.treatment"), aa.arth.abund.glmer))
plot(Effect("Genotype", aa.arth.abund.glmer)) # arthropod abundance varied 3.8-fold among the most disparate willow genotypes.
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.arth.abund.glmer)) # arthropod abundance increased 1.9-fold in the presence of aphids at increasing distances from ant mounds, but there was little effect of ant-mound distance when there were no aphids.

plot(Effect(c("Aphid.treatment","Genotype"), aa.arth.abund.glmer)) # GxAphid marginal interaction is predominantly influenced by genotypes T and J, which had much higher arthropod abundances in the presence of aphids.

aa.abund.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.arth.abund.glmer)) %>% mutate(response = "arthropod_abund")
aa.abund.ExE <- as.data.frame(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.arth.abund.glmer)) %>% mutate(response = "arthropod_abund")

# Calculate R2
aa.arth.abund.up <- update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (Aphid.treatment|Genotype))

summary(aa.arth.abund.up)

(aa.arth.abund.R2 <- var.table(aa.arth.abund.up, experiment = "ant-aphid"))

## ant-aphid: arthropod richness analysis ----
plot(total.rich ~ total.abund, aa.arth.df)

# GLMM
aa.arth.rich.glmer <- glmer(total.rich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                           #(1|plant_ID) +
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

## Likelihood ratio tests
(aa.rich.3 <- drop1(aa.arth.rich.glmer, test = "Chisq"))
(aa.rich.2 <- drop1(update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aa.rich.1 <- drop1(update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq")) # only Genotype effect

#aa.arth.rich.anova <- anova.table(aa.arth.rich.glmer, test = "Chisq", experiment = "ant-aphid")

## Effects
Effect("Genotype", aa.arth.rich.glmer) # arthropod richness varied 2.7-fold among willow genotypes.

aa.rich.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.arth.rich.glmer)) %>% mutate(response = "arthropod_rich") 

## Calculate R2
aa.arth.rich.up <- update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + (1|Genotype))

(aa.arth.rich.R2 <- var.table(aa.arth.rich.up, experiment = "ant-aphid"))

## ant-aphid: rarefied arthropod richness analysis ----

hist(filter(aa.arth.df, total.abund > 1)$total.rarerich)

# GLMM. Note that logit-transformation (common for proportion data) did not qualitatively affect the outcome or the residuals
aa.arth.rarerich.glmer <- lmer(total.rarerich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                              #(1|plant_ID) +
                              (1|Block/fact.Ant.mound.dist),
                            data = filter(aa.arth.df, total.abund > 1),
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"))
summary(aa.arth.rarerich.glmer)
plot(aa.arth.rarerich.glmer) 

## F tests
(aa.arth.rarerich.anova <- anova.table(aa.arth.rarerich.glmer, test = "F", type = 2, experiment = "ant-aphid")) # marginally significant aphid by genotype interactions

## Effects
plot(Effect(c("Aphid.treatment","Genotype"), aa.arth.rarerich.glmer))
plot(Effect(c("Aphid.treatment"), aa.arth.rarerich.glmer))
1-(0.5816897/0.6917055)

aa.rrich.GxE <- as.data.frame((Effect(c("Aphid.treatment","Genotype"), aa.arth.rarerich.glmer))) %>% mutate(response = "arthropod_rarerich")

## Calculate R2
aa.arth.rarerich.up <- update(aa.arth.rarerich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist) + Aphid.treatment + (Aphid.treatment|Genotype))

(aa.arth.rarerich.R2 <- var.table(aa.arth.rarerich.up, experiment = "ant-aphid"))

## ant-aphid: Hellinger distance ----
aa.hell <- decostand(aa.arth.12.pos[ ,aa.arth.names], method = "hellinger")

# test 3-way interaction
anova(rda(aa.hell ~ Ant.mound.dist*Aphid.treatment*Genotype, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# test 2-way interactions
anova(rda(aa.hell ~ (Ant.mound.dist + Aphid.treatment + Genotype)^2, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# which GxE are driving the GxEaphid effect?
aa.GxE.filt <- filter(aa.arth.12.pos, Genotype != "J")

aa.GxE.rda <- rda(decostand(aa.GxE.filt[ ,aa.arth.names], "hellinger") ~ Condition(Ant.mound.dist) + Aphid.treatment*Genotype, data = aa.GxE.filt)
plot(aa.GxE.rda)
anova(aa.GxE.rda, by = 'margin', permutations = how(block = aa.GxE.filt$Block, nperm = 999)) # removing Genotype J nullifies GxE effect 
anova(update(aa.GxE.rda, .~. -Aphid.treatment:Genotype), by = 'margin', permutations = how(block = aa.GxE.filt$Block, nperm = 999)) # still main effects of Genotype and aphid treatment though.


# test main effects
anova(rda(aa.hell ~ Ant.mound.dist + Aphid.treatment + Genotype, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# test Ant.mound.dist
aa.plots.hell <- betadisper(vegdist(aa.hell, method = "euclidean"), aa.arth.12.pos$Plot_code, bias.adjust = TRUE)

aa.plots.hell.centr <- data.frame(aa.plots.hell$centroids, 
                             id = rownames(aa.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Ant.mound.dist")) %>%
  mutate(Block = as.factor(Block),
         Ant.mound.dist = as.numeric(Ant.mound.dist))

anova(rda(aa.plots.hell$centroids ~ Block + Ant.mound.dist, data = aa.plots.hell.centr), by = "margin", permutations = how(block = aa.plots.hell.centr$Block, nperm = 999))

## assess assumptions for significant models. No evidence of overdispersion
aa.hell.geno <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$Genotype, bias.adjust = TRUE)
boxplot(aa.hell.geno)
plot(aa.hell.geno)
permutest(aa.hell.geno, permutations = how(block = aa.arth.12.pos$fact.Ant.mound.dist, nperm = 999))

aa.hell.GxE <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$GxE, bias.adjust = TRUE)
boxplot(aa.hell.GxE)
plot(aa.hell.GxE)
permutest(aa.hell.GxE, permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

aa.hell.aphid <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$Aphid.treatment, bias.adjust = TRUE)
boxplot(aa.hell.aphid)
plot(aa.hell.aphid)
permutest(aa.hell.aphid, permutations = how(block = aa.arth.12.pos$fact.Ant.mound.dist, nperm = 999))

## ant-aphid: Identify most abundant taxonomic groups ----
colSums(aa.arth.df[ ,aa.arth.names])
round(colSums(aa.arth.df[ ,aa.arth.names])/sum(colSums(aa.arth.df[ ,aa.arth.names]))*100,0)

## ant-aphid: Caloptilia analysis ----

# GLMM
aa.LTF.abund.glmer <- glmer(LTF_Caloptilia ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                              (1|plant_ID) +
                              (1|Block/fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.LTF.abund.glmer), correlation = TRUE)
overdisp_fun(aa.LTF.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(aa.LTF.abund.glmer) 

## likelihood-ratio tests
(LTF.2 <- drop1(aa.LTF.abund.glmer, test = "Chisq"))
(LTF.1 <- drop1(update(aa.LTF.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))
#(aa.LTF.abund.anova <- anova.table(aa.LTF.abund.glmer, test = "Chisq", experiment = "ant-aphid"))

aa.g.LTF <- as.data.frame(Effect(c("Genotype"), aa.LTF.abund.glmer)) %>%
  select(Genotype, aa.LTF.fit = fit)
plot(Effect(c("Ant.mound.dist"), aa.LTF.abund.glmer))
plot(Effect(c("Ant.mound.dist","Aphid.treatment"), aa.LTF.abund.glmer)) # expected response of LTF to distance and aphids

aa.LTF.abund.up <- update(aa.LTF.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (scale(Ant.mound.dist)|Genotype))

(aa.LTF.abund.R2 <- var.table(aa.LTF.abund.up, experiment = "ant-aphid"))

## ant-aphid: Aphids - non-Aphis analysis ----

# GLMM
aa.Aphids_nonAphis.glmer <- glmer(Aphididae ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                    (1|plant_ID) +
                                    (1|Block/fact.Ant.mound.dist),
                                  data = aa.arth.df,
                                  contrasts = list(Aphid.treatment = "contr.sum",
                                                   Genotype = "contr.sum"),
                                  family = "poisson",
                                  control=glmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e5)))
summary(aa.Aphids_nonAphis.glmer)
overdisp_fun(aa.Aphids_nonAphis.glmer) # accounted for overdispersion by modelling individual-level random effect.
plot(aa.Aphids_nonAphis.glmer) 

## likelihood ratio tests
(aphid.3 <- drop1(aa.Aphids_nonAphis.glmer, test = "Chisq"))
(aphid.2 <- drop1(update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aphid.1 <- drop1(update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))
#(aa.Aphids_nonAphis.anova <- anova.table(aa.Aphids_nonAphis.glmer, test = "Chisq", experiment = "ant-aphid"))

Effect("Genotype", aa.Aphids_nonAphis.glmer)

aa.g.Aphids_nonAphis <- as.data.frame(Effect("Genotype", aa.Aphids_nonAphis.glmer)) %>% # genotypes varied 35-fold in density of non-Aphis aphids.
  select(Genotype, aa.Aphids_nonAphis.fit = fit)
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.Aphids_nonAphis.glmer)) # marginal effect
plot(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer))

aa.Aphids_nonAphis.up <- update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + Aphid.treatment + (Aphid.treatment|Genotype))
summary(aa.Aphids_nonAphis.up)

(aa.Aphids_nonAphis.R2 <- var.table(aa.Aphids_nonAphis.up, experiment = "ant-aphid"))

## ant-aphid: leafhopper analysis ----
# GLMM
aa.leafhopper.glmer <- glmer(leafhopper ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                               #(1|plant_ID) +
                               (1|Block) +
                               (1|Block:fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "poisson",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
print(summary(aa.leafhopper.glmer), correlation = TRUE)
overdisp_fun(aa.leafhopper.glmer) # no overdispersion
plot(aa.leafhopper.glmer) 

# likelihood ratio tests
(lh.2 <- drop1(aa.leafhopper.glmer, test = "Chisq")) 
(lh.1 <- drop1(aa.leafhopper.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype, test = "Chisq")) 

#(aa.leafhopper.anova <- anova.table(aa.leafhopper.glmer, test = "Chisq", experiment = "ant-aphid"))

plot(Effect("Genotype",aa.leafhopper.glmer))

aa.leafhopper.up <- update(aa.leafhopper.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + (1|Genotype))

(aa.leafhopper.R2 <- var.table(aa.leafhopper.up, experiment = "ant-aphid"))

## ant-aphid: spiders analysis ----
# GLMM
aa.spiders.glmer <- glmer(spiders ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                            #(1|plant_ID) +
                            (1|Block) +
                            (1|Block:fact.Ant.mound.dist),
                          data = aa.arth.df,
                          contrasts = list(Aphid.treatment = "contr.sum",
                                           Genotype = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
print(summary(aa.spiders.glmer), correlation = TRUE)
overdisp_fun(aa.spiders.glmer) # no overdispersion
plot(aa.spiders.glmer) 

# likelihood ratio tests
(spid.2 <- drop1(aa.spiders.glmer, test = "Chisq"))
(spid.1 <- drop1(update(aa.spiders.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))
#(aa.spiders.anova <- anova.table(aa.spiders.glmer, test = "Chisq", experiment = "ant-aphid"))

aa.g.spid <- as.data.frame(Effect("Genotype", aa.spiders.glmer)) %>%
  select(Genotype, aa.spid.fit = fit)

aa.spiders.up <- update(aa.spiders.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + Aphid.treatment + (Aphid.treatment|Genotype))

(aa.spiders.R2 <- var.table(aa.spiders.up, experiment = "ant-aphid"))

## ant-aphid: black ant analysis ----

# GLMM
aa.ant_black.glmer <- glmer(ant_black ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                              #(1|plant_ID) +
                              (1|Block) +
                              (1|Block:fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.ant_black.glmer), correlation = TRUE)
overdisp_fun(aa.ant_black.glmer) # no evidence of overdispersion
plot(aa.ant_black.glmer) 

# likelihood ratio tests
(ant.2 <- drop1(aa.ant_black.glmer, test = "Chisq")) 
(ant.1 <- drop1(update(aa.ant_black.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq")) 
#(aa.ant_black.anova <- anova.table(aa.ant_black.glmer, test = "Chisq", experiment = "ant-aphid"))

plot(Effect("Aphid.treatment",aa.ant_black.glmer))
plot(Effect("Genotype",aa.ant_black.glmer))

aa.ant_black.up <- update(aa.ant_black.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + Aphid.treatment + (1|Genotype) + (1|plant_ID)) # need to add to account for overdispersion

(aa.ant_black.R2 <- var.table(aa.ant_black.up, experiment = "ant-aphid"))

## ant-aphid: leaftier Tortricidae analysis ----

# GLMM
aa.leaftier_Tortricid.glmer <- glmer(leaftier_Tortricid ~  scale(Ant.mound.dist)*Aphid.treatment + Genotype +
                                       #(1|plant_ID) +
                                       (1|Block) +
                                       (1|Block:fact.Ant.mound.dist),
                                     data = aa.arth.df,
                                     contrasts = list(Aphid.treatment = "contr.sum",
                                                      Genotype = "contr.sum"),
                                     family = "poisson",
                                     control=glmerControl(optimizer="bobyqa",
                                                          optCtrl=list(maxfun=2e5)))
summary(aa.leaftier_Tortricid.glmer)
overdisp_fun(aa.leaftier_Tortricid.glmer) # no strong evidence of overdispersion
plot(aa.leaftier_Tortricid.glmer) 

# likelihood ratio tests
(lt.2 <- drop1(aa.leaftier_Tortricid.glmer, test = "Chisq"))
(lt.1 <- drop1(update(aa.leaftier_Tortricid.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment), test = "Chisq"))
#(aa.leaftier_Tortricid.anova <- anova.table(aa.leaftier_Tortricid.glmer, test = "Chisq", experiment = "ant-aphid"))

plot(Effect(c("Ant.mound.dist"),aa.leaftier_Tortricid.glmer))
aa.g.tort <- as.data.frame(Effect(c("Genotype"),aa.leaftier_Tortricid.glmer)) %>%
  select(Genotype, aa.tort.fit = fit)
plot(Effect(c("Ant.mound.dist","Aphid.treatment"),aa.leaftier_Tortricid.glmer ))

aa.leaftier_Tortricid.up <- update(aa.leaftier_Tortricid.glmer, .~. -Genotype + (1|Genotype) + (1|plant_ID))

(aa.leaftier_Tortricid.R2 <- var.table(aa.leaftier_Tortricid.up, experiment = "ant-aphid"))

##### ant-aphid: manuscript plots ----
pd <- 0.2
#cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")
p.size <- 5
l.size <- 1.5

## richness
rich.sum <- as.data.frame(Effect(c("Genotype"), aa.arth.rich.glmer)) %>%
  mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)])) # re-order Genotype according to increasing richness

aa.r <- rich.sum %>% ggplot(aes(x = Genotype, y = fit)) +  geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + ylab("Arthropod richness") + xlab("Willow genotype")

## rarefied richness
aa.rr <- as.data.frame(Effect(c("Aphid.treatment"), aa.arth.rarerich.glmer)) %>% mutate(Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = fit, shape = Aphid.treatment, fill = Aphid.treatment))  + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd/2)+ geom_point(size = p.size) + ylab("Rarefied richness") + xlab("Aphid treatment") + scale_shape_manual(values = c(23,21), guide = "none") + scale_fill_manual(values = c("gray50", "white"), guide = "none")

## arthropod abundance
aa.aG <- as.data.frame(Effect(c("Genotype"), aa.arth.abund.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.sum$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) +  geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + xlab("Willow genotype") + scale_y_continuous(name = "Arthropod abundance", breaks = c(0,2,4,6,8,10), limits = c(0,10.5))

aa.arth.abund.fact <- glmer(total.abund ~ ord.Ant.mound.dist*Aphid.treatment + Genotype + (1|plant_ID) + (1|Block/fact.Ant.mound.dist), data = aa.arth.df, contrasts = list(Aphid.treatment = "contr.sum", Genotype = "contr.sum"), family = "poisson", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # model ant mound distance as an ordered factor to make plotting easier. Same qualitative result as when I model it as a continuous predictor. Also, I removed interaction with willow genotype since it was not significant and it allowed the model to converge.

aa.aExE <- as.data.frame(Effect(c("Aphid.treatment","ord.Ant.mound.dist"), aa.arth.abund.fact)) %>% mutate(ord.Ant.mound.dist = ordered(as.numeric(as.character(ord.Ant.mound.dist))), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = ord.Ant.mound.dist, y = fit, fill = Aphid.treatment, shape = Aphid.treatment)) + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd, position = position_dodge(width = 0.5)) + geom_point(size = p.size, position = position_dodge(width = 0.5)) + ylab("No. of individuals") + xlab("Distance from ant mound (m)") + scale_shape_manual(name = "Aphid treatment", values = c(23,21)) + scale_fill_manual(name = "Aphid treatment", values = c("gray50", "white")) + scale_y_continuous(name = "Arthropod abundance", breaks = c(0,2,4,6,8,10), limits = c(0,10.5)) + theme(legend.justification = c(0,1), legend.position = c(0,1)); aa.aExE

## ordination of community compositino
aa.hell.GxE.rda <- rda(aa.hell ~ Condition(Block) + Condition(Ant.mound.dist) +
                         GxE, 
                       data = aa.arth.12.pos)
summary(aa.hell.GxE.rda)$cont$importance[ ,1:2] 

aa.hell.GxE.df <- data.frame(scores(aa.hell.GxE.rda, display = "cn")) %>%
  mutate(factor = levels(aa.arth.12.pos$GxE)) %>%
  separate(factor, into = c("Genotype","Aphid.treatment")) %>%
  mutate(Genotype = as.factor(Genotype),
         Aphid.treatment = as.factor(Aphid.treatment),
         GxE.sig.lab = ifelse(Genotype == "J","yes","no"))

aa.hell.GxE.df.arrows <- data.frame(x = aa.hell.GxE.df$RDA1[11:20], y = aa.hell.GxE.df$RDA2[11:20], xend = aa.hell.GxE.df$RDA1[1:10], yend = aa.hell.GxE.df$RDA2[1:10], Genotype = aa.hell.GxE.df$Genotype[11:20]) %>% mutate(GxE.sig.lab = ifelse(Genotype == "J","yes","no"))

aa.hell.GxE.sites <- data.frame(scores(aa.hell.GxE.rda, display = "sites"))

(aa.arth.ord.GxE <- ggplot(data = aa.hell.GxE.sites, aes(x = RDA1, y = RDA2)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_point(shape = 1, color = "gray") +
  geom_segment(data = aa.hell.GxE.df.arrows, aes(x = x, xend = xend, y = y, yend = yend, linetype = GxE.sig.lab), color = "black") + #, color = Genotype
  geom_point(data = aa.hell.GxE.df, aes(x = RDA1, y = RDA2, fill = Aphid.treatment, shape = Aphid.treatment), color = "black", size = 5) + 
  geom_text(data = aa.hell.GxE.df, 
            aes(x = RDA1, y = RDA2, label = Genotype), size = 3) +
  scale_fill_manual(values = c("gray50","white"))+ #cbPal.10) +
  scale_color_manual(values = c("gray50","white")) +# cbPal.10) +
  scale_shape_manual(values = c(23,21)) +
  scale_linetype_manual(values = c("dotted","solid"), guide = "none") +
  ylab("RDA 2 (3%)") + xlab("RDA 1 (9%)") + theme(legend.position = "none"))

## visualize non-Aphis aphids driving GxE and which genotype
aa.Aphid.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer)) %>% mutate(sig.GxE = ifelse(Genotype == "J", "yes","no"), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = fit, group = Genotype, fill = Aphid.treatment, shape = Aphid.treatment)) + geom_line(aes(linetype = sig.GxE), size = 1) + geom_point(size = p.size) + xlab("Aphid treatment") + scale_fill_manual(values = c("gray50","white"), guide = "none") + scale_shape_manual(values = c(23,21), guide = "none") + geom_text(aes(label = Genotype), size = 3) + scale_linetype_manual(values = c("dotted","solid"), guid = "none") + scale_y_continuous(name = "Abundance of other aphids", breaks = c(0,1,2), limits = c(0,2)); aa.Aphid.GxE

## establish multi-panel plot
aa.arth.p <- plot_grid(aa.r, aa.aG, aa.aExE, aa.rr, aa.arth.ord.GxE, aa.Aphid.GxE, labels = "AUTO", ncol = 3, align = 'hv'); aa.arth.p

save_plot(aa.arth.p, filename = "fig_1_aa_arth_comm.png", base_height = 8.5, base_width = 11)

##### ant-aphid: Aphis farinosa and Formica obscuripes responses ----

## Aphis farinosa
aa.Aphis.p <- as.data.frame(Effect("Genotype", aa.Aphis.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.sum$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + ylab(expression(paste(italic("Aphis farinosa")," abundance"))) + xlab("Willow genotype") + scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1,10,20, 30), limits = c(0.001,30), labels = c(0.001, 0.01, 0.1, 1,10,20, 30)); aa.Aphis.p

aa.Fobs.p <- aa.arth.df %>% group_by(Genotype, Aphid.treatment) %>% summarise(mean.Fobs = mean(ant_F_obscuripes))  %>% ungroup() %>% mutate(Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = mean.Fobs, group = Genotype, shape = Aphid.treatment, fill = Aphid.treatment)) + geom_line() + geom_point(size = p.size) + geom_text(aes(label = Genotype), size = 3) + scale_shape_manual(values = c(23,21), guide = "none") + scale_fill_manual(values = c("gray50", "white"), guide = "none") + xlab("Aphid treatment") + ylab(expression(paste(italic("Formica obscuripes")," abundance"))) + scale_y_continuous(breaks = c(0,0.2, 0.4, 0.6), limits = c(0,0.6)); aa.Fobs.p #+ scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1), limits = c(0,1))

## ant-aphid: plots for supp mat -----
# consider turning one significant genotype to black solid line and circle, while everything else is grey with white circle and different line types.
cp.keyEaa <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#009E73", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")
aa.key.E <- data.frame(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer), response = "aphids") %>% 
  ggplot(aes(x = Aphid.treatment, y = fit, group = Genotype, color = Genotype, fill = Genotype)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Genotype", values = cbPal.10) + ylab("No. of aphid individuals") + xlab("Aphid treatment")  + scale_fill_manual(name = "Genotype", values = cp.keyEaa) + scale_y_log10() + theme(legend.justification = c(1,0), legend.position = c(1,0)) #+ geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w, position = position_dodge(width = pd)) # na.rm = TRUE, position = position_dodge(width = pd), 

cp.keyGaa <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#FFFFFF") # consider altering colors to match other wind graphs
aa.key.G <- data.frame(Effect(c("Genotype"), aa.LTF.abund.glmer), response = "leaf-mining moths") %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), aa.leaftier_Tortricid.glmer), response = "leaf-tiering moths")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), aa.ant_black.glmer), response = "ants")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), aa.spiders.glmer), response = "spiders")) %>%
  #bind_rows(., data.frame(Effect(c("Genotype"), aa.Aphids_nonAphis.glmer), response = "aphids")) %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), aa.leafhopper.glmer), response = "leafhoppers")) %>% 
  ggplot(aes(x = Genotype, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Willow genotype") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyGaa) + theme(legend.justification = c(1,0), legend.position = c(1,0))

## ant-aphid: functional trait models ----

# 2012
aa.arth.abund.2012 <- glmer(total.abund ~ Height + all.shoot.count + leaf_WC + leaf_trichome.density + Aphid.treatment*Ant.mound.dist + (1|Genotype) +
                           (1|plant_ID) + 
                           (1|Block) + (1|Block:fact.Ant.mound.dist),
                         data = filter(aa.arth.df, Height > 0, all.shoot.count > 0, leaf_WC > 0, leaf_trichome.density != "NA"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
overdisp_fun(aa.arth.abund.2012) # need to add individual-level random effect
lm.beta.lmer(aa.arth.abund.2012) # biggest effect from plant height but not by much.
summary(aa.arth.abund.2012) 
drop1(aa.arth.abund.2012, test = "Chisq") # Height is significant and trichome density is marginal
anova(aa.arth.abund.2012, update(aa.arth.abund.2012, .~. -(1|Genotype))) # still a highly significant effect of willow genotype, suggesting we didn't capture the key traits mediating variation in arthropod abundance.

## Wind: Plots of G and E effects on euclidean distance ----
# GxE
w.13.euc.GxE.rda <- rda(w.13.euc.log ~ Condition(Block) + 
                           GxE, data = w.arth.13.full)
summary(w.13.euc.GxE.rda)$cont$importance[ ,1:2] 

w.13.euc.GxE.df <- data.frame(scores(w.13.euc.GxE.rda, display = "cn")) %>%
  mutate(factor = levels(droplevels(w.arth.13.full$GxE))) %>%
  separate(factor, into = c("Genotype","Wind.Exposure")) 

plot(w.13.euc.GxE.rda, scaling = 1, type = "n")
#points(w.13.euc.GxE.rda, display = "sites", col = "gray50")
arrows(x1 = w.13.euc.GxE.df$RDA1[1:10], y1 = w.13.euc.GxE.df$RDA2[1:10], x0 = w.13.euc.GxE.df$RDA1[11:20], y0 = w.13.euc.GxE.df$RDA2[11:20], length = 0.1, lwd = 1.5)
text(x = w.13.euc.GxE.df$RDA1, y = w.13.euc.GxE.df$RDA2, labels = w.13.euc.GxE.df$Genotype, cex = 1.5)
plot(w.13.euc.GxE.rda, scaling = 1, display = "sp")

# Genotype effect
w.13.euc.pres <- w.arth.13.full[ ,wind.arth.names]
w.13.euc.pres[w.13.euc.pres > 0] <- 1

w.13.euc.geno.rda <- rda(w.13.euc.pres ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = w.arth.13.full)
summary(w.13.euc.geno.rda)$cont$importance[ ,1:2] 

plot(w.13.euc.geno.rda, display = "sp", scaling = 1)
points(w.13.euc.geno.rda, display = "sites", col = "gray50")
ordiellipse(w.13.euc.geno.rda, groups = w.arth.13.full$Genotype, col = "gray50",draw = "polygon")
text(w.13.euc.geno.rda, display = "cn", labels = levels(w.arth.13.full$Genotype))

## Wind: Plots of G and E effects on hellinger distances ----

# GxE
table(w.arth.13.pos2$Genotype, w.arth.13.pos2$Wind.Exposure)

w.13.hell.GxE.rda <- rda(w.13.hell.pos2 ~ Condition(Block) + 
                           GxE, data = w.arth.13.pos2)
summary(w.13.hell.GxE.rda)$cont$importance[ ,1:2] 

w.13.hell.GxE.df <- data.frame(scores(w.13.hell.GxE.rda, display = "cn", scaling = 3)) %>%
  mutate(factor = levels(droplevels(w.arth.13.pos2$GxE))) %>%
  separate(factor, into = c("Genotype","Wind.Exposure")) 

plot(w.13.hell.GxE.rda, display = "sp", scaling = 2)#type = "n")
#points(w.13.hell.GxE.rda, display = "sites", col = "gray50")
arrows(x1 = w.13.hell.GxE.df$RDA1[1:7], y1 = w.13.hell.GxE.df$RDA2[1:6], x0 = w.13.hell.GxE.df$RDA1[c(8,10:12,14:16)], y0 = w.13.hell.GxE.df$RDA2[c(8,10:12,14:16)], length = 0.1, lwd = 1.5)
text(x = w.13.hell.GxE.df$RDA1[-c(9,13,17)], y = w.13.hell.GxE.df$RDA2[-c(9,13,17)], labels = w.13.hell.GxE.df$Genotype[-c(9,13,17)], cex = 1.5)
text(x = w.13.hell.GxE.df$RDA1[c(9,13,17)], y = w.13.hell.GxE.df$RDA2[c(9,13,17)], labels = w.13.hell.GxE.df$Genotype[c(9,13,17)], cex = 1.5, col = "gray50")


# Genotype effect
w.13.hell.geno.rda <- rda(w.13.hell.pos2 ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = w.arth.13.pos2)
summary(w.13.hell.geno.rda)$cont$importance[ ,1:2] 

plot(w.13.hell.geno.rda, type = "n")
points(w.13.hell.geno.rda, display = "sites", col = "gray50")
ordiellipse(w.13.hell.geno.rda, groups = w.arth.13.pos2$Genotype, col = "gray50",draw = "polygon")
text(w.13.hell.geno.rda, display = "cn", labels = levels(w.arth.13.pos2$Genotype))

## Wind effect
w.13.hell.wind.rda <- rda(w.13.plots.hell$centroids ~ Block + Wind.Exposure, data = w.13.plots.centr.hell)
summary(w.13.hell.wind.rda)$cont$importance[ ,1:2] 

plot(w.13.hell.wind.rda, type = "n")
points(w.13.hell.wind.rda, display = "sites", col = "gray50")
ordiellipse(w.13.hell.wind.rda, groups = w.13.plots.centr.hell$Wind.Exposure, draw = "polygon", col = "gray50")
text(x = scores(w.13.hell.wind.rda, display = "cn")[11:12,1], y = scores(w.13.hell.wind.rda, display = "cn")[11:12,2], labels = c("Exposed","Unexposed"))

## Ant-aphid: Plots of G and E effects on hellinger distances ----
aa.hell.GxE.rda <- rda(aa.hell ~ Condition(Block) + Condition(Ant.mound.dist) +
                           GxE, data = aa.arth.12.pos)
summary(aa.hell.GxE.rda)$cont$importance[ ,1:2] 

aa.hell.GxE.df <- data.frame(scores(aa.hell.GxE.rda, display = "cn")) %>%
  mutate(factor = levels(aa.arth.12.pos$GxE)) %>%
  separate(factor, into = c("Genotype","Aphid.treatment")) %>%
  mutate(Genotype = as.factor(Genotype),
         Aphid.treatment = as.factor(Aphid.treatment))


plot(aa.hell.GxE.rda, type = "n")
points(aa.hell.GxE.rda, display = "sites", col = "gray50")
arrows(x1 = aa.hell.GxE.df$RDA1[1:10], y1 = aa.hell.GxE.df$RDA2[1:10], x0 = aa.hell.GxE.df$RDA1[11:20], y0 = aa.hell.GxE.df$RDA2[11:20], length = 0.1, lwd = 1.5)
text(x = aa.hell.GxE.df$RDA1, y = aa.hell.GxE.df$RDA2, labels = aa.hell.GxE.df$Genotype, cex = 1.5)

## Maximum aphid abundance as the predictor of F. obscuripes abundance ----
aa.Fobs.abund.aphis <- glmer(ant_F_obscuripes > 0 ~ scale(aphid_Aphis) + Genotype +
                               #(1|plant_ID) +
                               (1|Block/fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "binomial",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
summary(aa.Fobs.abund.aphis)
plot(aa.Fobs.abund.aphis) 

drop1(aa.Fobs.abund.aphis, test = "Chisq")

plot(Effect("aphid_Aphis", aa.Fobs.abund.aphis))

## Print Summary ----
w.arth.means <- bind_rows(w.abund.GxE, w.rich.GxE, w.rrich.GxE) %>% mutate(experiment = "wind")

aa.arth.means <- bind_rows(aa.abund.GxE, aa.rich.GxE, aa.rrich.GxE, aa.Aphis.G, aa.ant.GxE, aa.ant.ExE) %>% mutate(experiment = "ant_aphid")

arth.means <- bind_rows(w.arth.means, aa.arth.means)

arth.anovas <- bind_rows(arth.abund.anova,
                          arth.rich.anova,
                          arth.rarerich.anova,
                          LTF.abund.anova,
                          gall.abund.anova,
                          leaftier.abund.anova,
                          spider.abund.anova,
                          aphid.abund.anova,
                          aa.arth.abund.anova,
                          aa.arth.rich.anova,
                          aa.arth.rarerich.anova,
                          aa.Aphis.anova,
                          aa.Fobs.abund.anova,
                          aa.LTF.abund.anova,
                          aa.Aphids_nonAphis.anova,
                          aa.leafhopper.anova,
                          aa.spiders.anova,
                          aa.ant_black.anova,
                          aa.leaftier_Tortricid.anova) %>%
  select(experiment, response, term, test, df, Df.res, statistic, p.value, Pr..Chisq.)

arth.R2s <- bind_rows(arth.abund.R2,
                      arth.rich.R2,
                      arth.rarerich.R2,
                      LTF.abund.R2,
                      gall.abund.R2,
                      leaftier.abund.R2,
                      spider.abund.R2,
                      #aphid.abund.R2, # nothing significant
                      aa.arth.abund.R2,
                      aa.arth.rich.R2,
                      aa.arth.rarerich.R2,
                      aa.Aphis.R2,
                      aa.Fobs.abund.R2,
                      aa.LTF.abund.R2,
                      aa.Aphids_nonAphis.R2,
                      aa.leafhopper.R2,
                      aa.spiders.R2,
                      aa.ant_black.R2,
                      aa.leaftier_Tortricid.R2)

write.csv(arth.means, "arth.means.csv")
write.csv(arth.anovas, "arth.anovas.csv")
write.csv(arth.R2s, "arth.R2s.csv")

##### probably not needed ----

## Arthropod correlations across experiments
tort.df <- left_join(g.tort, aa.g.tort)
plot(tort.fit ~ aa.tort.fit, tort.df)
with(tort.df, cor.test(tort.fit, aa.tort.fit))

LTF.df <- left_join(g.LTF, aa.g.LTF) %>%
  left_join(., g.LTF.2012)
plot(LTF.fit ~ aa.LTF.fit, LTF.df)
plot(LTF.fit.2012 ~ aa.LTF.fit, LTF.df)
plot(LTF.fit ~ LTF.fit.2012, LTF.df)

with(LTF.df, cor.test(LTF.fit, aa.LTF.fit))
with(LTF.df, cor.test(LTF.fit.2012, aa.LTF.fit))
with(LTF.df, cor.test(LTF.fit.2012, LTF.fit))

aphid.df <- left_join(g.Aphids_nonAphis, aa.g.Aphids_nonAphis)
plot(Aphids_nonAphis.fit ~ aa.Aphids_nonAphis.fit, aphid.df)
text(x = aphid.df$aa.Aphids_nonAphis.fit, y = aphid.df$Aphids_nonAphis.fit, labels = aphid.df$Genotype)
with(aphid.df, cor.test(Aphids_nonAphis.fit, aa.Aphids_nonAphis.fit))

spid.df <- left_join(g.spid, aa.g.spid)
plot(spid.fit ~ aa.spid.fit, spid.df)

## Old analyses of key taxa ----


######## END  ###############################################
#### Wind 2012: community analysis ----
library(mvabund)
w.13.mv <- mvabund(w.arth.13.full[ ,wind.arth.names])

w.13.glm <- manyglm(w.13.mv ~ Wind.Exposure*Genotype, data = w.arth.13.full, family = "negative.binomial")
plot(w.13.glm, which = 1:3)

anova.manyglm(w.13.glm, block = w.arth.13.full$Block)


# separate NMDS for each year suggest we have insufficient data, so we pooled data across years
metaMDS(decostand(w.arth.12.pos[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")
metaMDS(decostand(w.arth.13.pos[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")

w.arth.all <- wind.arth.df %>%
  select(Block, Wind.Exposure, Genotype, plant_ID, Gracilliaridae_miner:Spider) %>%
  group_by(Block, Wind.Exposure, Genotype, plant_ID) %>%
  summarise_each(funs(sum))
w.arth.all$total <- rowSums(w.arth.all[ ,wind.arth.names])

w.arth.all.sub <- filter(w.arth.all, total > 0)

metaMDS(decostand(w.arth.all.sub[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")

## Hellinger distance
with(w.arth.12.pos, table(Genotype, Wind.Exposure)) # most made up from unexposed plots...

# for avoiding pseudoreplication while testing for wind effect
w.geno.13 <- w.arth.13.full %>%
  select(-plant_ID, -Wind.Exposure, -Plot_code) %>%
  group_by(Block, Genotype) %>%
  summarise_each(funs(sum)) %>%
  filter(total.abund > 0)

w.arth.pos <- wind.arth.df %>%
  filter(wind.arth.df, total.abund > 0)
w.13.mds <- metaMDS(decostand(w.effect.13[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")
plot(w.13.mds, display = c("sites"))
points(w.13.mds, display = "sites", col = as.numeric(w.effect.13$Wind.Exposure))
points(w.13.mds, display = "sites", col = as.numeric(w.effect.13$Genotype))

w.13.mds.geno <- metaMDS(decostand(w.geno.13[ ,wind.arth.names], method = "hellinger"), distance = "euclidean", k = 1)
plot(w.13.mds.geno, display = c("sites"))
points(w.13.mds.geno, display = "sites", col = as.numeric(w.geno.13$Genotype))

# GxE effect
w.12.rda.comp <- rda(decostand(w.arth.12.pos[ ,wind.arth.names], method = "hellinger") ~ Genotype*Wind.Exposure, data = w.arth.12.pos)
anova(w.12.rda.comp, by = "margin", permutations = how(block = w.arth.12.pos$Block)) # no GxE effect

anova(update(w.12.rda.comp, .~. -Wind.Exposure:Genotype), by = "margin", permutations = how(block = w.arth.12.pos$Block)) # no GxE effect

# wind effect, summarized at plot level
w.12.rda.wind.comp <- rda(decostand(w.effect.12[ ,wind.arth.names], method = "hellinger") ~ Wind.Exposure, data = w.effect.12)
anova(w.12.rda.wind.comp, by = "margin", permutations = how(block = w.effect.12$Block)) # marginal wind effect
plot(w.12.rda.wind.comp, display = c("bp","sp","sites"))
plot(w.12.rda.wind.comp, display = c("sites"), type = "n")
text(w.12.rda.wind.comp, display = c("sites"), labels = w.effect.12$Wind.Exposure)

plot(w.12.rda.wind.comp, display = "sites")
ordiellipse(w.12.rda.wind.comp, groups = w.effect.12$Wind.Exposure, 
            kind = "sd", draw = "polygon", 
            col= "gray50", #"gainsboro", 
            border = NA, label = T)


w.12.wind.nmds <- metaMDS(decostand(w.effect.12[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")
# plots of species never make sense to me in NMDS, they usually have the least abundant species on the far ends of the plots.
plot(w.12.wind.nmds, display = "sites", type = "n")
text(w.12.wind.nmds, display = "sites", labels = w.effect.12$Wind.Exposure)

# test Genotype
table(w.arth.12.pos$Genotype) # still likely have enough reps per genotype

w.geno.12 <- filter(w.geno.12, total.abund > 0)

w.12.rda.geno.comp <- rda(decostand(w.geno.12[ ,wind.arth.names], method = "hellinger") ~ Genotype, data = w.geno.12)
anova(w.12.rda.geno.comp, by = "margin", permutations = how(block = w.geno.12$Block)) # Genotype explains 15% of the variance

plot(w.12.rda.geno.comp, display = c("bp","sp","sites"))
plot(w.12.rda.geno.comp, display = c("sites"), type = "n")
text(w.12.rda.geno.comp, display = c("sites"), labels = w.geno.12$Genotype)

plot(w.12.rda.geno.comp, display = "sites", type = "n")
ordiellipse(w.12.rda.geno.comp, groups = w.geno.12$Genotype, 
            kind = "sd", draw = "polygon", 
            col= "gray50", #"gainsboro", 
            border = NA, label = T)

## Wind 2013: community analysis ----

## Euclidean distance - joint absenses from sites influence analysis
w.comm.13.full <- w.arth.13.full[ ,wind.arth.names]

# test GxE
w.13.rda <- rda(log(w.comm.13.full+1) ~ Wind.Exposure*Genotype, data = w.arth.13.full) # log transform to remove influence of outliers
anova(w.13.rda, by = "margin", permutations = how(blocks = w.arth.13.full$Block)) # no GxE

# test G
w.13.rda.geno <- update(w.13.rda, .~. -Wind.Exposure*Genotype + Genotype)
anova(w.13.rda.geno, by = "margin", permutations = how(blocks = w.arth.13.full$Plot_code)) # G effect (sometimes marginal). explains 7% of the variance.
plot(w.13.rda.geno, display = c("bp","sp","sites"))

# test wind exposure. Data summarized at plot level
w.13.rda.wind <- rda(log(w.effect.13[ ,wind.arth.names]+1) ~ Wind.Exposure, data = w.effect.13)
anova(w.13.rda.wind, by = "margin", permutations = how(block = w.effect.13$Block)) # wind effect, explains 20% of the variance
plot(w.13.rda.wind, display = c("bp","sp","sites"))

# test Genotype, summarized at block level
w.geno.13 <- filter(w.geno.13, total.abund > 0)
w.13.rda.geno2 <- rda(log(w.geno.13[ ,wind.arth.names]+1) ~ Genotype, data = w.geno.13)
anova(w.13.rda.geno2, by = "margin", permutations = how(block = w.geno.13$Block)) # genotype effect, explains 14% of the variance
plot(w.13.rda.geno2, display = c("bp","sp","sites"))

## Hellinger distance

# GxE effect
w.13.rda.comp <- rda(decostand(w.arth.13.pos[ ,wind.arth.names], method = "hellinger") ~ Genotype*Wind.Exposure, data = w.arth.13.pos)
anova(w.13.rda.comp, by = "margin", permutations = how(block = w.arth.13.pos$Block)) # no GxE effect

# wind effect, summarized at plot level
w.13.rda.wind.comp <- rda(decostand(w.effect.13[ ,wind.arth.names], method = "hellinger") ~ Wind.Exposure, data = w.effect.13)
anova(w.13.rda.wind.comp, by = "margin", permutations = how(block = w.effect.13$Block)) # wind explains 26% of the variance in composition
plot(w.13.rda.wind.comp, display = c("bp","sp","sites"))
plot(w.13.rda.wind.comp, display = c("sites"), type = "n")
text(w.13.rda.wind.comp, display = c("sites"), labels = w.effect.13$Wind.Exposure)

plot(w.13.rda.wind.comp, display = c("sites"), type = "n")
text(w.13.rda.wind.comp, display = "sites", labels = w.effect.13$Wind.Exposure)
ordiellipse(w.13.rda.wind.comp, groups = w.effect.13$Wind.Exposure)

w.13.wind.nmds <- metaMDS(decostand(w.effect.13[ ,wind.arth.names], method = "hellinger"), distance = "euclidean")
# plots of species never make sense to me in NMDS, they usually have the least abundant species on the far ends of the plots.
plot(w.13.wind.nmds, display = "sites", type = "n")
text(w.13.wind.nmds, display = "sites", labels = w.effect.13$Wind.Exposure)

# test Genotype
table(w.arth.13.pos$Genotype) # still likely have enough reps per genotype

w.13.rda.geno.comp <- rda(decostand(w.arth.13.pos[ ,wind.arth.names], method = "hellinger") ~ Genotype, data = w.arth.13.pos)
anova(w.13.rda.geno.comp, by = "margin", permutations = how(block = w.arth.13.pos$Plot_code)) # no genotype effect

plot(w.13.rda.geno.comp, display = c("bp","sp","sites"))
plot(w.13.rda.geno.comp, display = c("sites"), type = "n")
text(w.13.rda.geno.comp, display = c("sites"), labels = w.arth.13.pos$Genotype)




# dissimilarity matrix, full dataset

w.dis.13 <- vegdist(w.arth.13.pos[ ,wind.arth.names],
                    method = "horn")
w.dis.13.sub <- vegdist(w.effect.13[ ,wind.arth.names],
                        method = "horn")

## Note that adonis tests terms sequentially in the model...
# Testing interaction: block level as strata
# non-significant effect
adonis(w.dis.13.pres ~ Wind.Exposure*Genotype,
       data = w.arth.13.full,#w.arth.13.pos,
       strata = w.arth.13.full$Block, #w.arth.13.pos$Block,
       method = "euclidean")
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
## ant-aphid: herbivore abundance analysis ----
hist(aa.arth.df$herb.abund.nonAphis)

# GLMM
aa.herb.nonAphis.glmer <- glmer(herb.abund.nonAphis ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                  (1|plant_ID) +
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

(aa.herb.nonAphis.anova <- anova.table(aa.herb.nonAphis.glmer, test = "Chisq", experiment = "ant-aphid"))

plot(Effect(c("Ant.mound.dist","Aphid.treatment"), aa.herb.nonAphis.glmer))
as.data.frame(Effect(c("Ant.mound.dist","Aphid.treatment"), aa.herb.nonAphis.glmer))

aa.herb.nonAphis.up <- update(aa.herb.nonAphis.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (Aphid.treatment|Genotype) -(1|plant_ID))

plotFEsim(FEsim(update(aa.herb.nonAphis.up, .~. -1)))
plotREsim(REsim(aa.herb.nonAphis.up))
with(aa.arth.df, interaction.plot(Aphid.treatment, Genotype, herb.abund.nonAphis))
with(aa.arth.df, interaction.plot(Ant.mound.dist, Aphid.treatment, herb.abund.nonAphis))

aa.herb.nonAphis.R2 <- var.table(aa.herb.nonAphis.up, experiment = "ant-aphid")

## ant-aphid: non-aphid herbivores
aa.arth.df$herb_nonAphid <- rowSums(
  select(aa.arth.df, psyllid,leafhopper,froghopper,sawfly_larva,gall_R_salicisbattatus,grasshopper,leaftier_Tortricid,LTF_Caloptilia)
)

aa.herb_nonAPHID.glmer <- glmer(herb_nonAphid ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                  (1|plant_ID) +
                                  (1|Block/fact.Ant.mound.dist),
                                data = aa.arth.df,
                                contrasts = list(Aphid.treatment = "contr.sum",
                                                 Genotype = "contr.sum"),
                                family = "poisson",
                                control=glmerControl(optimizer="bobyqa",
                                                     optCtrl=list(maxfun=2e5)))
print(summary(aa.herb_nonAPHID.glmer), correlation = TRUE)
overdisp_fun(aa.herb_nonAPHID.glmer) # overdispersed
plot(aa.herb_nonAPHID.glmer) 

aa.herb_nonAPHID.anova <- anova.table(aa.herb_nonAPHID.glmer, test = "Chisq", experiment = "ant-aphid")

plot(Effect(c("Ant.mound.dist","Aphid.treatment"), aa.herb_nonAPHID.glmer))

## ant-aphid: predator abundance analysis ----

# GLMM
aa.pred.nonFobs.glmer <- glmer(pred.abund.nonFobs ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                 #(1|plant_ID) +
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

plot(Effect("Aphid.treatment",aa.pred.nonFobs.glmer))
plot(Effect("Genotype",aa.pred.nonFobs.glmer))

aa.pred.nonFobs.up <- update(aa.pred.nonFobs.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + Aphid.treatment + (1|Genotype))

plotFEsim(FEsim(aa.pred.nonFobs.up))
plotREsim(REsim(aa.pred.nonFobs.up))
with(aa.arth.df, plot(pred.abund.nonFobs ~ Aphid.treatment + Genotype))
with(aa.arth.df, interaction.plot(Ant.mound.dist, Aphid.treatment, pred.abund.nonFobs))

aa.pred.nonFobs.R2 <- var.table(aa.pred.nonFobs.up, experiment = "ant-aphid")

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
  select(Year:plant_ID, Gracilliaridae_miner:Spider, wind.arth.nonzero, wind.arth.rich) %>%
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
  select(Year:plant_ID, Gracilliaridae_miner:Spider) %>%
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

## more old ----

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

## old again----
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

## Wind: community analysis ----

w.rda <- rda(w.comm ~ Wind.Exposure + Genotype + Condition(Block), data = filter(wind.arth.df, Year == "2013"))
anova(w.rda, by = "margin")
plot(w.rda, display = c("sp","bp"))

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

## Wind: genetic correlation analysis ----
tort.g <- as.data.frame(Effect("Genotype", leaftier.abund.glmer)) %>% mutate(arthropod = "Tortricidae_leaftier", experiment = "wind")
cecid.g <- as.data.frame(Effect("Genotype", gall.abund.glmer)) %>% mutate(arthropod = "Cecidomyiidae_gall", experiment = "wind")
calop.g <- as.data.frame(Effect("Genotype", LTF.abund.glmer)) %>% mutate(arthropod = "Caloptilia", experiment = "wind")

herb.g <- bind_rows(tort.g, cecid.g, calop.g)

herb.g.corr <- tidyr::spread(select(herb.g, Genotype, arthropod, fit), arthropod, fit)

scatterplotMatrix(herb.g.corr[ ,-1])
corr.test(herb.g.corr[ ,-1], method = "spearman") # little correlation in responses among genotypes.
ggplot(herb.g, aes(x = Genotype, y = fit, color = arthropod)) +
  geom_point(position = position_dodge()) 

## Wind: herbivore abundance analysis ----
hist(wind.arth.df$herb.abund)
table(wind.arth.df$herb.abund)

# GLMM
herb.abund.glmer <- glmer(herb.abund ~ (Wind.Exposure + Year + Genotype)^2 + # solves the problem with the model being identified
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|X) +
                            (1|plant_ID),
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
(herb.abund.anova <- anova.table(herb.abund.glmer, test = "Chisq", experiment = "wind"))


## Calculate R2 for significant predictors
herb.abund.up <- update(herb.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + (1|Genotype) - (1|X)) # removing individual-level effect, because the R2 calculation already accounts for overdispersion.

(herb.abund.R2 <- var.table(herb.abund.up, "wind"))

## Wind: predator abundance analysis ----

# GLMM
pred.abund.glmer <- glmer(pred.abund ~ Wind.Exposure*Year + Genotype +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_ID),
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

## Examine correlations between plant traits and willow community ----
nonsense <- with(wind.arth.df, which(Height > 80 | leaf_WC < 0))
colSums(filter(wind.arth.df, Year == "2012")[ ,wind.arth.names])
# 2012
corr.test(x = select(filter(wind.arth.df[-nonsense, ], Year == "2012"), Height:all.shoot.count, leaf_trichome.density, leaf_WC),
          y = select(filter(wind.arth.df[-nonsense, ], Year == "2012"), total.abund,total.rich, Tortricidiae_leaftier, Aphididae, Gracilliaridae_miner), 
          adjust = "none",
          method = "spearman")

colSums(filter(wind.arth.df, Year == "2013")[ ,wind.arth.names])
plot(Tortricidiae_leaftier ~ Genotype, filter(wind.arth.df, Year == "2013"))
# 2013
corr.test(x = select(filter(wind.arth.df[-nonsense, ], Year == "2013"), Height:all.shoot.count, leaf_WC, SLA, leaf_C_N),
          y = select(filter(wind.arth.df[-nonsense, ], Year == "2013"), total.abund, total.rich, Cecidomyiidae_gall, Tortricidiae_leaftier, Gracilliaridae_miner, Spider),
          adjust = "none",
          method = "spearman")

# variance inflation factor analysis. Note that choice of response variable does not matter.
vif(lm(total.abund ~ Height + all.shoot.count + leaf_WC + leaf_trichome.density, data = filter(wind.arth.df[-nonsense, ], Year == "2012"))) # hmm, relatively low vif among plant growth traits. Still, it made sense to remove at least average shoot length which was highly correlated with plant height.
vif(lm(total.abund ~ Height + leaf_WC + leaf_C_N, data = filter(wind.arth.df[-nonsense, ], Year == "2013"))) # still low vif values overall with all of the traits. But, I felt justified in removing SLA since it was highly correlated with leaf water content and we have data on it from both years. I also removed the other correlated growth traits since height gives an overall better picture of plant size and we have data on it from the common garden experiment.

## Wind community analyses

## ant-aphid: Aphis growth rates analysis ----
# issues with model convergence.
# calculate average aphid growth rate for each plant.
aa.aphid.GR.sum <- aa.aphid.GR %>%
  group_by(Block, Genotype, Ant.mound.dist, Aphid.treatment, plant_ID) %>%
  dplyr::summarise(mean.Aphid.GR = mean(Aphis.growth.rate, na.rm = TRUE)) %>%
  mutate(fact.Ant.mound.dist = as.factor(Ant.mound.dist))
hist(aa.aphid.GR.sum$mean.Aphid.GR) # note though that growth rates are virtually all negative.

aphid.GR.lmer <- lmer(mean.Aphid.GR ~ scale(Ant.mound.dist) +
                        #(1|Genotype) + 
                        (1|Block/fact.Ant.mound.dist),
                      aa.aphid.GR.sum[1:141, ],
                      contrasts = list(Genotype = "contr.sum"))
summary(aphid.GR.lmer)
plot(aphid.GR.lmer)

aphid.GR.anova <- anova.table(aphid.GR.lmer, test = "F", type = 3, experiment = "ant-aphid") # effect of genotype on aphid growth rate.

aphid.GR.up <- update(aphid.GR.lmer, .~. -scale(Ant.mound.dist)*Genotype + (1|Genotype))

aphid.GR.R2 <- var.table(aphid.GR.up, experiment = "ant-aphid")

## Wind: euclidean distance 2012 ----
# only a marginal effect of wind exposure

w.12.euc.log <- log(w.arth.12.full[ ,wind.arth.names]+1)

# test G and GxE. 
adonis(w.12.euc.log ~ Wind.Exposure*Genotype, 
       data = w.arth.12.full, 
       method = "euclidean", 
       permutations = how(block = w.arth.12.full$Block, nperm = 999)) # no G or GxE

# test wind effect
w.12.plots <- betadisper(vegdist(w.12.euc.log, method = "euclidean"),  w.arth.12.full$Plot_code, bias.adjust = TRUE)

w.12.plots.centr <- data.frame(w.12.plots$centroids, 
                               id = rownames(w.12.plots$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.12.plots$centroids ~ Block + Wind.Exposure, 
       data = w.12.plots.centr, 
       method = "euclidean", 
       permutations = how(block = w.12.plots.centr$Block, nperm = 999)) # marginal effect of wind exposure

w.euc.12.wind <- betadisper(vegdist(w.12.plots$centroids, method = "euclidean"), w.12.plots.centr$Wind.Exposure, bias.adjust = TRUE)
boxplot(w.euc.12.wind)
plot(w.euc.12.wind) # likely a dispersion effect
anova(w.euc.12.wind) # marginal effect, suggesting that marginal effect of wind on centroid location could be due to overdispersion.

##  Wind: euclidean distance 2013 ----
# no G, wind, or GxE effects
w.13.euc.log <- log(w.arth.13.full[ ,wind.arth.names]+1)
w.13.euc.hell <- decostand(w.arth.13.full[ ,wind.arth.names], method = "hellinger")

# test G and GxE
adonis(w.13.euc.pres ~ Wind.Exposure*Genotype, 
       data = w.arth.13.full, 
       method = "euclidean", 
       permutations = how(block = w.arth.13.full$Block, nperm = 9999)) # no GxE

w.13.euc.geno <- betadisper(vegdist(w.13.euc.pres, method = "euclidean"), w.arth.13.full$Genotype, bias.adjust = TRUE)
boxplot(w.13.euc.geno) # doesn't seem like there is too much variation in dispersion
plot(w.13.euc.geno)
permutest(w.13.euc.geno, permutations = how(block = w.arth.13.full$Plot_code, nperm = 999)) # used plot code as blocking factor to account for effects of wind exposure. Significant dispersion between genotypes, suggesting that genotype effect on centroid location may be due to genotype effect.

# test wind effect
w.13.plots <- betadisper(vegdist(w.13.euc.pres, method = "euclidean"),  w.arth.13.full$Plot_code, bias.adjust = TRUE)

w.13.plots.centr <- data.frame(w.13.plots$centroids, 
                               id = rownames(w.13.plots$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.13.plots$centroids ~ Block + Wind.Exposure, 
       data = w.13.plots.centr,
       method = "euclidean", 
       permutations = how(block = w.13.plots.centr$Block, nperm = 999)) # no effect of wind exposure

## ant-aphid: Euclidean distance analysis ----
aa.euc <- log(aa.arth.df[ ,aa.arth.names]+1)

# test: Genotype, Aphid.treatment:Genotype, and 3-way interaction 
adonis(aa.euc ~ Ant.mound.dist*Aphid.treatment*Genotype, 
       data = aa.arth.df,
       method = "euclidean",
       permutations = how(block = aa.arth.df$Block, nperm = 999))

# shuffle order to test: Aphid.treatment and Aphid.treatment:Ant.mound.dist
adonis(aa.euc ~ Genotype*Ant.mound.dist*Aphid.treatment, 
       data = aa.arth.df,
       method = "euclidean",
       permutations = how(block = aa.arth.df$Block, nperm = 999))

# shuffle order to test: Genotype:Ant.mound.dist. Need to test Ant.mound.dist separately to account for appropriate degrees of freedom
adonis(aa.euc ~ Aphid.treatment*Genotype*Ant.mound.dist, 
       data = aa.arth.df,
       method = "euclidean",
       permutations = how(block = aa.arth.df$Block, nperm = 999))

# test Ant.mound.dist
aa.plots <- betadisper(vegdist(aa.euc, method = "euclidean"), aa.arth.df$Plot_code, bias.adjust = TRUE)

aa.plots.centr <- data.frame(aa.plots$centroids, 
                             id = rownames(aa.plots$centroids)) %>%
  separate(col = id, into = c("Block","Ant.mound.dist")) %>%
  mutate(Block = as.factor(Block),
         Ant.mound.dist = as.numeric(Ant.mound.dist))

adonis(aa.plots$centroids ~ Block + Ant.mound.dist, 
       data = aa.plots.centr, 
       method = "euclidean", 
       permutations = how(block = aa.plots.centr$Block, nperm = 999))

## assess assumptions for significant models. No evidence of overdispersion
aa.euc.geno <- betadisper(vegdist(aa.euc, "euclidean"), aa.arth.df$Genotype, bias.adjust = TRUE)
boxplot(aa.euc.geno)
plot(aa.euc.geno)
permutest(aa.euc.geno, permutations = how(block = aa.arth.df$fact.Ant.mound.dist, nperm = 999)) # significant overdispersion

aa.euc.GxE <- betadisper(vegdist(aa.euc, "euclidean"), aa.arth.df$GxE, bias.adjust = TRUE)
boxplot(aa.euc.GxE)
plot(aa.euc.GxE)
permutest(aa.euc.GxE, permutations = how(block = aa.arth.df$Block, nperm = 999)) # significant overdispersion

aa.euc.ExE <- betadisper(vegdist(aa.euc, "euclidean"), with(aa.arth.df, interaction(Aphid.treatment, Ant.mound.dist)), bias.adjust = TRUE)
boxplot(aa.euc.ExE)
plot(aa.euc.ExE)
permutest(aa.euc.ExE, permutations = how(block = aa.arth.df$Block, nperm = 999)) # significant overdispersion
#RsquareAdj(aa.euc.12)

#anova(aa.euc.12, by = "margin", permutations = how(block = aa.arth.df$Block, nperm = 999)) # 3-way interaction: not sig.

#anova(update(aa.euc.12, .~. -Genotype:Aphid.treatment:Ant.mound.dist), by = "margin", permutations = how(block = aa.arth.df$Block, nperm = 999)) # test 2-way interactions. GxAphid and Aphid-by-Ant distance are signficant

#anova(update(aa.euc.12, .~. -Genotype*Aphid.treatment*Ant.mound.dist + Genotype + Aphid.treatment + Condition(Ant.mound.dist)), by = "margin", permutations = how(block = aa.arth.df$Block, nperm = 999)) # Genotype significant and aphid treatment is marginal

#aa.euc.12.Anteff <- rda(log(aa.effect.12[ ,aa.arth.names]+1) ~ Ant.mound.dist + Condition(Block), data = aa.effect.12)
#RsquareAdj(aa.euc.12.Anteff)

#anova(aa.euc.12.Anteff, by = "margin", permutations = how(block = aa.effect.12$Block, nperm = 999)) # test Wind: not sig.