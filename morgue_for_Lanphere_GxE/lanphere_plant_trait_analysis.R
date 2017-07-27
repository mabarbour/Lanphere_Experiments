############################################
## Description: This script analyzes data for plant traits of the Lanphere experiments.
## Code author: Matt Barbour
## Email: barbour@zoology.ubc.ca
############################################

#### Load required libraries and fucnctions ----
source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must be loaded before dplyr
library(dplyr)
library(tidyr)
#library(reshape2)
#library(ggplot2)
library(lme4)
library(effects) # calculating mean and confidence intervals of treatment and genotype effects.
library(psych) # for correlation tests
#library(car) # for Anova function
#library(broom) # for tidying up model outputs
library(RCurl) # get code directly from github

## source and parse github code
script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/model_diagnostic_functions.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#### Load required data sets ----

## wind plant traits
w.trait.df <- read.csv('~/Lanphere_Experiments/final_data/wind_trait_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block), # necessary for properly modelling as a random effect.
         Year = as.factor(Year),
         #plant.position = as.character(plant.position),
         X = as.factor(X)) # necessary for accounting for overdispersion
         #treat.tmp = ifelse(treatment == "Exposed","E","U"),
         #Plot_code = paste(Block, treat.tmp, sep = "."),
         #Sample = paste(Block, treat.tmp, genotype, plant.position, sep = ""))
glimpse(w.trait.df)
summary(w.trait.df$Height)

## wind soil properties
#w.soil <- read.csv('final_data/wind_soil_df.csv') %>%
 # tbl_df() 
#summary(w.soil$avg.moisture.vwc) # moisture ranges from 0.04 to 0.20

## Root C:N
#rootCN <- read.csv('Output/RootCNdataLanphereDunes.csv') %>%
#  tbl_df() 

# discuss with angelica on limiting to biologically reasonable values
#table(rootCN$Nitrogen_percent) # 3 zeros
#table(rootCN$Carbon_percent) # 3 zeros, 1 value at 115.7%
#table(rootCN$Nitrogen_mg) # 3 zeros and 3 values close to zero
#table(rootCN$Carbon_mg)
#table(rootCN.df$root_CN)

# need to split up plant ID into Block, treatment, genotype and position
#rootCN$TreatGeno <- gsub("[[:digit:]]", "", rootCN$Sample)
#rootCN$Block__Position <- gsub("[^[:digit:]]", "_", rootCN$Sample)

#rootCN.df <- rootCN %>%
  #separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  #filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) %>% 
  #separate(Block__Position, into = c("Block","Plant.Position")) %>%
  #mutate(Genotype = as.factor(Genotype),
     #    Wind.Exposure = as.factor(Wind.Exposure),
    #     Block = as.factor(Block),
         #sample.id = paste(Block, Wind.Exposure, Genotype, Plant.Position, sep = ""),
   #      Year = "2013",
         #Sample = interaction(Block, Wind.Exposure, Genotype, Plant.Position),
  #       root_CN = Carbon_mg/Nitrogen_mg) %>%
 # filter(Nitrogen_percent > 0 & Carbon_percent > 0 & Carbon_percent < 100 & root_CN < 700) # restricting to biologically reasonable values
#rootCN.df$Year <- as.factor(rootCN.df$Year)

#table(rootCN.df$Nitrogen_percent) # 3 zeros
#table(rootCN.df$Carbon_percent) # 3 zeros, 1 value at 115.7%
#table(rootCN.df$Nitrogen_mg) # 3 zeros and 3 values close to zero
#table(rootCN.df$Carbon_mg)
#table(rootCN.df$root_CN)

## combine wind datasets
#w.trait.df <- left_join(w.trait.df, select(w.soil, Plot_code, Total.N:percent.TOM, avg.moisture.vwc, nut.PC1), by = "Plot_code") %>%
 # left_join(., select(rootCN.df, Sample, Year, root_CN), by = c("Sample","Year"))
#glimpse(w.trait.df)
#w.trait.df$Year <- as.factor(w.trait.df$Year)

## ant-aphid above-ground plant traits
aa.trait.df <- read.csv('~/Lanphere_Experiments/final_data/ant_aphid_trait_df.csv') %>%
  tbl_df() %>%
  filter(Year == "2012") %>%
  mutate(fact.Ant.mound.dist = as.factor(Ant.mound.dist), # necessary for modelling as a random effect
         Block = as.factor(Block),
         X = as.factor(X),# necessary for modelling individual-level random effects to account for overdispersion.
         Plot_code = paste(Block, Ant.mound.dist, sep = "_")) %>%
  select(-Year)
glimpse(aa.trait.df)
summary(aa.trait.df$Height)

#aa.soil <- read.csv('final_data/ant_aphid_soil_2013.csv') %>%
 # tbl_df() 
#glimpse(aa.soil)
#summary(aa.soil$moisture.vwc.July9_13) # soil moisture only ranges from 0.03 to 0.05

## combine ant-aphid datasets
#aa.df <- left_join(aa.trait.df, select(aa.soil, plot.id, moisture.vwc.July9_13), by = "Plot_code")
#glimpse(aa.df)

#### Examine correlation among plant traits ----
## Wind: phenotypic correlations
# 2012
scatterplotMatrix( ~ Height + all.shoot.avg.length + all.shoot.count + leaf_trichome.density + leaf_WC, data = filter(w.trait.df, Year == "2012"))
corr.test(x = select(filter(w.trait.df, Year == "2012", Wind.Exposure == "Exposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))
corr.test(x = select(filter(w.trait.df, Year == "2012", Wind.Exposure == "Unexposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))

# 2013
scatterplotMatrix( ~ Height + all.shoot.avg.length + all.shoot.count + leaf_WC + leaf_C_N + SLA, data = filter(w.trait.df, Year == "2013"))
corr.test(x = select(filter(w.trait.df, Year == "2013", Wind.Exposure == "Exposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_WC, SLA, leaf_C_N))
corr.test(x = select(filter(w.trait.df, Year == "2013", Wind.Exposure == "Unexposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_WC, SLA, leaf_C_N))

## Ant-aphid: phenotypic correlations
scatterplotMatrix( ~ Height + mature.shoot.avg.length + all.shoot.count + leaf_trichome.density + leaf_WC, data = filter(aa.trait.df, Year == "2012"))
corr.test(x = select(filter(aa.trait.df, Year == "2012"), Height, mature.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))

## Wind: plant height analysis ----

# LMMM. 
height.lmer <- lmer(Height ~ Wind.Exposure*Year*Genotype + #scale(avg.moisture.vwc) + 
                        (1|Block) + 
                        (1|Block:Wind.Exposure) +
                        (1|plant_ID),
                    data = w.trait.df,
                    contrasts = list(Wind.Exposure = "contr.sum",
                                     Year = "contr.sum",
                                     Genotype = "contr.sum"))
summary(height.lmer)
plot(height.lmer) # diamond shaped...

# Anova table
(height.anova <- anova.table(height.lmer, type = 2, experiment = "wind"))

# Effects
Effect(c("Wind.Exposure","Genotype","Year"), height.lmer)
Effect(c("Wind.Exposure","Year"), height.lmer) # effect of unexposed plots increased from 1.2-fold in 2012 to 2-fold in 2013.
Effect(c("Year"), height.lmer) # plants were 39% shorter in 2013 vs. 2012. 
Effect(c("Wind.Exposure"), height.lmer) # plants were 29% shorter in unexposed plots
Effect(c("Genotype"), height.lmer) # plants varied 2-fold in height among the most disparate Genotypes.

w.height.resp <- as.data.frame(Effect(c("Genotype"), height.lmer)) %>% bind_rows(., as.data.frame(Effect(c("Wind.Exposure","Year"), height.lmer))) %>% mutate(response = "height")

# calculate R2
height.up <- update(height.lmer, .~. -Wind.Exposure*Year*Genotype + (1|Genotype) + Wind.Exposure*Year) # simplify model
lm.beta.lmer(height.up)

(height.R2 <- var.table(height.up, experiment = "wind"))

## Wind: shoot count analysis ----

# GLMM
shoot.count.glmer <- glmer(all.shoot.count ~ Wind.Exposure*Genotype*Year + #scale(avg.moisture.vwc) + # moisture was the best predictor of the soil traits
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           (1|plant_ID),
                         data = w.trait.df,
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(shoot.count.glmer)
overdisp_fun(shoot.count.glmer) # no overdispersion
plot(shoot.count.glmer)

## Likelihood-ratio tests
(shoot.count.3 <- drop1(shoot.count.glmer, test = "Chisq"))
(shoot.count.2 <- drop1(update(shoot.count.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq"))
(shoot.count.1 <- drop1(update(shoot.count.glmer, .~. -Wind.Exposure:Genotype:Year -Wind.Exposure:Genotype -Wind.Exposure:Year -Genotype:Year), test = "Chisq"))
#(shoot.count.anova <- anova.table(shoot.count.glmer, test = "Chisq", type = 2, experiment = "wind"))

## Calculate effects
Effect(c("Wind.Exposure","Year"), shoot.count.glmer) # effect of unexposed plots increased from 1.1-fold in 2012 to 1.6-fold in 2013.
Effect(c("Year"), shoot.count.glmer)
plot(Effect(c("Year","Genotype"), shoot.count.glmer) )
Effect(c("Wind.Exposure"), shoot.count.glmer) # plants produced 22% fewer shoots.
Effect(c("Genotype"), shoot.count.glmer) # plants varied nearly 2.3-fold in shoot count among the most disparate Genotypes.
shoot.count.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), shoot.count.glmer)) %>% mutate(response = "shoot.count")

shoot.count.GxEy <- as.data.frame(Effect(c("Year","Genotype"), shoot.count.glmer)) %>% mutate(response = "shoot.count")
shoot.count.GxEy %>% ggplot(aes(x = Year, y = fit, color = Genotype, group = Genotype)) + geom_line()

# calculate variance components for simplified model
shoot.count.up <- update(shoot.count.glmer, .~. -Wind.Exposure*Year*Genotype + (Year|Genotype) + Wind.Exposure*Year)

(shoot.count.R2 <- var.table(shoot.count.up, experiment = "wind"))

## Wind: average shoot length analysis ----

# LMM. 
shoot.length.lmer <- lmer(log(all.shoot.avg.length) ~ Wind.Exposure*Genotype*Year +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|plant_ID),
                          data = w.trait.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"))
summary(shoot.length.lmer)
plot(shoot.length.lmer) # looks okay

# Kenward-Roger tests 
(shoot.length.anova <- anova.table(shoot.length.lmer, type = 2, experiment = "wind"))

# effects
Effect(c("Wind.Exposure"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants shoots grew 28% shorter.
Effect(c("Year"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants shoots grew 1.7-fold longer in 2012.
Effect(c("Genotype"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants varied 2.2-fold in shoot length among the most disparate Genotypes.

plot(Effect(c("Wind.Exposure","Genotype"), shoot.length.lmer, transformation = list(link = log, inverse = exp)))

## Calculate R2 for significant predictors
# no interaction effects were significant so I didn't calculate R2
shoot.length.up <- update(shoot.length.lmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

(shoot.length.R2 <- var.table(shoot.length.up, experiment = "wind"))

## Wind: leaf water content ----

leaf_WC.lmer <- lmer(log(leaf_WC) ~ Wind.Exposure*Genotype*Year + #scale(avg.moisture.vwc) +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|plant_ID),
                         data = w.trait.df,
                         contrast = list(Wind.Exposure = "contr.sum",
                                         Genotype = "contr.sum",
                                         Year = "contr.sum"))
summary(leaf_WC.lmer)  # note that plant_ID explains the least amount of variance
plot(leaf_WC.lmer) # looks good

# Kenward-Roger test
(leaf_WC.anova <- anova.table(leaf_WC.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype","Year"), leaf_WC.lmer, transformation = list(link = log, inverse = exp)))
plot(Effect(c("Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp)))

WC.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_WC")

## Calculate R2
leaf_WC.up <- update(leaf_WC.lmer, .~. -Wind.Exposure*Genotype*Year + Year + (0+Year|Genotype))

(leaf_WC.R2 <- var.table(leaf_WC.up, experiment = "wind"))

## Wind: trichome density analysis ----
leaf_trichome.density.glmer <- glmer(
  leaf_trichome.density ~ Wind.Exposure*Genotype +
    (1|Block) + 
    (1|Block:Wind.Exposure) +
    (1|plant_ID), # need individual-level random effect to account for overdispersion
  data = w.trait.df, 
  family = "poisson",
  contrasts = list(Wind.Exposure = "contr.sum",
                   Genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(leaf_trichome.density.glmer)
overdisp_fun(leaf_trichome.density.glmer)
plot(leaf_trichome.density.glmer) 

# Likelihood-ratio tests
(td.2 <- drop1(leaf_trichome.density.glmer, test = "Chisq"))
(td.1 <- drop1(update(leaf_trichome.density.glmer, .~. -Wind.Exposure:Genotype), test = "Chisq"))
#trichome.density.anova <- anova.table(leaf_trichome.density.glmer, test = "Chisq", experiment = "wind")

# effects
plot(Effect("Genotype", leaf_trichome.density.glmer)) # willows varied 46-fold in leaf trichome density.

td.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_trichome.density.glmer)) %>% mutate(response = "leaf_trichome.density")

## Calculate R2
leaf_trichome.density.up <- update(leaf_trichome.density.glmer, .~. -Wind.Exposure*Genotype + (1|Genotype)) 

(leaf_trichome.density.R2 <- var.table(leaf_trichome.density.up, experiment = "wind"))

## Wind: leaf C:N analysis ----

# LMM
leaf_CN.lmer <- lmer(log(leaf_C_N) ~ Wind.Exposure*Genotype +
                            (1|Block) + (1|Block:Wind.Exposure),
                          data = w.trait.df,
                     contrasts = list(Wind.Exposure = "contr.sum",
                                      Genotype = "contr.sum"))
print(summary(leaf_CN.lmer), correlation = TRUE)
plot(leaf_CN.lmer) # looks pretty good

## Kenwar-Roger test
(leaf_CN.anova <- anova.table(leaf_CN.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype"), leaf_CN.lmer, transformation = list(link = log, inverse = exp))) # willows varied 1.6-fold in leaf C:N

CN.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_CN.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_CN")

## Calculate R2
leaf_CN.up <- update(leaf_CN.lmer, .~. -Wind.Exposure*Genotype + (1|Genotype)) 

(leaf_CN.R2 <- var.table(leaf_CN.up, experiment = "wind"))

## Wind: SLA analysis ----
hist(w.trait.df$SLA)

# LMM
SLA.lmer <- lmer(SLA ~ Wind.Exposure*Genotype +
                            (1|Block) + (1|Block:Wind.Exposure),
                          data = w.trait.df,
                      contrasts = list(Wind.Exposure = "contr.sum",
                                       Genotype = "contr.sum"))
summary(SLA.lmer)
plot(SLA.lmer) # residuals look okay

## Kenward-Roger test
(SLA.anova <- anova.table(SLA.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype"), SLA.lmer)) # Genotypes varied 1.5-fold in SLA.

SLA.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), SLA.lmer)) %>% mutate(response = "SLA")

# Calculate R2
SLA.up <- update(SLA.lmer, .~. -Wind.Exposure*Genotype + (1|Genotype))

(SLA.R2 <- var.table(SLA.up, experiment = "wind"))

## Wind: Root C:N analysis ----
rootCN.lmer <- lmer(log(root_CN) ~ Wind.Exposure*Genotype + (1|Block) + (1|Block:Wind.Exposure),
                    data = w.trait.df, # results are qualitatively the same if I filter C:N less than 100
                    contrasts = list(Wind.Exposure = "contr.sum",
                                     Genotype = "contr.sum")) # restricting to biologically reasonable values
summary(rootCN.lmer)
plot(rootCN.lmer)

(rootCN.anova <- anova.table(rootCN.lmer, type = 2, experiment = "wind"))

rootCN.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), rootCN.lmer)) %>% mutate(response = "root_CN")
write.csv(rootCN.GxE, "rootCN.means.csv")

rootCN.up <- update(rootCN.lmer, .~. -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))

(rootCN.R2 <- var.table(rootCN.up, experiment = "wind"))

## Wind: plots ----
library(cowplot)
cbPal.10 <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66")
ebar.w <- 0.05
l.size <- 1.5
alp <- 0.5
p.size <- 5

height.Gord <- as.data.frame(Effect("Genotype", height.lmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)]))

# height
w.h <- as.data.frame(Effect("Genotype", height.lmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + ylab("Plant height (cm)") + xlab("Willow Genotype") + scale_y_continuous(limits = c(9,35)); w.h

w.h.ExE <- as.data.frame(Effect(c("Year","Wind.Exposure"), height.lmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit, fill = Year, shape = Year)) + geom_line(aes(group = Year), linetype = "dotted")  + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w)+ geom_point(size = p.size)  + xlab("Wind Exposure") + theme(legend.justification = c(1,0), legend.position = c(1,0)) + scale_y_continuous(name = "Plant height (cm)", limits = c(9,35))  + scale_fill_manual(values = c("black","gray50")) + scale_shape_manual(values = c(21,24)); w.h.ExE

# shoot count
#as.data.frame(Effect(c("Genotype","Year"), shoot.length.lmer, transformation = list(link = log, inverse = exp))) %>% ggplot(aes(x = Year, y = fit)) + geom_point() + geom_line(aes(group = Genotype, color = Genotype)) #+ geom_errorbar(aes(ymax = upper, ymin = lower))

#as.data.frame(Effect(c("Year","Wind.Exposure"), shoot.count.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point() + geom_line(aes(group = Year)) + geom_errorbar(aes(ymax = upper, ymin = lower))

# trichomes
#w.trich <- as.data.frame(Effect("Genotype", leaf_trichome.density.glmer)) %>% ggplot(aes(x = Genotype, y = fit/1.1)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + scale_y_log10() + ylab(expression(paste("Trichome density (no. ",cm^-2,")"))) + xlab("Willow Genotype")

# leaf C:N
w.CN <- as.data.frame(Effect("Genotype", leaf_CN.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(height.Gord$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + scale_y_continuous(name = "Leaf C:N", limits = c(25,95)) + xlab("Willow Genotype"); w.CN #+ scale_y_log10() 

# root C:N
w.rootCN <- as.data.frame(Effect("Genotype", rootCN.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(height.Gord$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + scale_y_continuous(name = "Root C:N", limits = c(25,95)) + xlab("Willow Genotype"); w.rootCN #+ scale_y_log10() 
# leaf water content
#as.data.frame(Effect(c("Year","Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% ggplot(aes(x = Year, y = fit, group = Genotype, color = Genotype))  + geom_line() + geom_point(size = p.size) + ylab("Leaf water content") + xlab("Willow Genotype") + theme(legend.justification = c(1,0), legend.position = c(1,0)) + scale_color_manual(name = "Genotype", values = cbPal.10)#+ scale_linetype_discrete(name = "Genotype")

#w.WC <- as.data.frame(Effect(c("Year","Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% ggplot(aes(x = Genotype, y = fit, group = Year, shape = Year, fill = Year, color = Year)) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 0.5)) + geom_point(size = p.size, position = position_dodge(width = 0.5)) + scale_fill_manual(values = c("black","gray50"))+ scale_color_manual(values = c("black","gray50")) + scale_shape_manual(values = c(21,24)) + ylab("Leaf water content") + xlab("Willow Genotype") + theme(legend.justification = c(0,1), legend.position = c(0,1)) 

# run lanphere_soil_analysis script before plot_grid
w.soil.traits.p <- plot_grid(w.moist.p, w.N.p, w.h.ExE, w.h, w.CN, w.rootCN, labels = "AUTO", ncol = 2, align = 'hv'); w.soil.traits.p

save_plot("fig_6_wind_soil_traits.png", w.soil.traits.p, base_height = 11, base_width = 8.5)

## Ant-aphid: plant height analysis ----

# LMM
aa.height.lmer <- lmer(Height ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + 
                           (1|Block) + (1|Block:fact.Ant.mound.dist),
                         data = aa.trait.df, #filter(aa.df, Year == "2012"),
                         contrasts = list(Aphid.treatment = "contr.sum",
                                          Genotype = "contr.sum"))
summary(aa.height.lmer)
plot(aa.height.lmer)

## Kenward-Roger test
(aa.height.anova <- anova.table(aa.height.lmer, type = 2, experiment = "ant-aphid"))

# effects
plot(Effect("Genotype", aa.height.lmer)) # willows varied 2-fold in height among the most disparate genotypes.

aa.height.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.height.lmer)) %>% mutate(response = "height")


## Calculate R2
aa.height.up <- update(aa.height.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.height.R2 <- var.table(aa.height.up, experiment = "ant-aphid"))

## Ant-aphid: Shoot count analysis ----

# GLMM
aa.shoot.count.glmer <- glmer(
  all.shoot.count ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + 
    (1|Block) + (1|Block:fact.Ant.mound.dist),
                             data = aa.trait.df, #filter(aa.df, Year == "2012"),
                             family = "poisson",
                             contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
plot(aa.shoot.count.glmer) 
summary(aa.shoot.count.glmer)
overdisp_fun(aa.shoot.count.glmer) # no overdispersion

## Likelihood rato tests
(aa.shoot.count.3 <- drop1(aa.shoot.count.glmer, test = "Chisq"))
(aa.shoot.count.2 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.treatment:scale(Ant.mound.dist):Genotype), test = "Chisq"))
(aa.shoot.count.1 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.treatment:scale(Ant.mound.dist):Genotype -Aphid.treatment:scale(Ant.mound.dist) -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq"))
#aa.shoot.count.anova <- anova.table(aa.shoot.count.glmer, test = "Chisq", experiment = "ant-aphid")

# effects
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.shoot.count.glmer)) # don't quite understand these effects
plot(Effect(c("Aphid.treatment"), aa.shoot.count.glmer))
Effect(c("Genotype"), aa.shoot.count.glmer) # 1.8-fold variation

aa.shoot.count.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.shoot.count.glmer)) %>% mutate(response = "shoot.count")
aa.shoot.count.ExE <- as.data.frame(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.shoot.count.glmer)) %>% mutate(response = "shoot.count")

## Calculate R2 for significant predictors
aa.shoot.count.up <- update(aa.shoot.count.glmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + Aphid.treatment*scale(Ant.mound.dist) + (1|Genotype))

(aa.shoot.count.R2 <- var.table(aa.shoot.count.up, experiment = "ant-aphid"))

## Ant-aphid: Average shoot length analysis ----
hist(aa.df$mature.shoot.avg.length)

aa.mature.shoot.avg.length.lmer <- lmer(mature.shoot.avg.length ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Block) + (1|Block:Ant.mound.dist),
                                 data = aa.trait.df, #filter(aa.df, Year == "2012"),
                                 contrasts = list(Aphid.treatment = "contr.sum",
                                                  Genotype = "contr.sum"))
plot(aa.mature.shoot.avg.length.lmer) # not great...
summary(aa.mature.shoot.avg.length.lmer)

## Kenward-Roger
(aa.mature.shoot.avg.length.anova <- anova.table(aa.mature.shoot.avg.length.lmer, type = 2, experiment = "ant-aphid"))

## effects
Effect("Genotype", aa.mature.shoot.avg.length.lmer) # 2.5-fold

aa.shoot.length.GxE <- as.data.frame((Effect(c("Aphid.treatment","Genotype"), aa.mature.shoot.avg.length.lmer))) %>% mutate(response = "shoot.length")

## Calculate R2
aa.mature.shoot.avg.length.up <- update(aa.mature.shoot.avg.length.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.mature.shoot.avg.length.R2 <- var.table(aa.mature.shoot.avg.length.up, experiment = "ant-aphid"))

## Ant-Aphid: leaf water content analysis ----
hist(aa.df$leaf_WC)

## LMM. Had to fix a less complex model because we didn't have enough replication of all treatment combinations
aa.leaf_WC.lmer <- lmer(log(leaf_WC) ~ Aphid.treatment*Genotype*scale(Ant.mound.dist) + #Aphid.treatment*scale(Ant.mound.dist) + 
                          #scale(Ant.mound.dist)*Genotype + 
                          (1|Block) + (1|Block:fact.Ant.mound.dist), 
                        data = aa.trait.df, #filter(aa.df, Year == "2012"),
                        contrasts = list(Aphid.treatment = "contr.sum",
                                         Genotype = "contr.sum"))
summary(aa.leaf_WC.lmer)
plot(aa.leaf_WC.lmer) # looks okay

## Kenward-Roger test
(aa.leaf_WC.anova <- anova.table(aa.leaf_WC.lmer, type = 2, experiment = "ant-aphid")) # nothing significant, but genotype is the closest.

## effects
plot(Effect("Genotype",aa.leaf_WC.lmer, transformation = list(link = log, inverse = exp))) # no significant effect

aa.WC.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"),aa.leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_WC")

## Calculate R2
# Genotype was not significant, but I'm still going to calculate its R2.
aa.leaf_WC.up <- update(aa.leaf_WC.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist) - scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.leaf_WC.R2 <- var.table(aa.leaf_WC.lmer, experiment = "ant-aphid"))

## Ant-aphid: leaf trichome density ----

# GLMM. Again, needed to simplify due to lack of replication
aa.leaf_trichome.density.glmer <- glmer(leaf_trichome.density ~ (Aphid.treatment + Genotype + scale(Ant.mound.dist))^2 +
                                          (1|Block) + 
                                          (1|X) +
                                          (1|Block:Ant.mound.dist),
                                        data = aa.trait.df, 
                                        family = "poisson",
                                        contrasts = list(Aphid.treatment = "contr.sum",
                                                         Genotype = "contr.sum"),
                                        control=glmerControl(optimizer="bobyqa",
                                                             optCtrl=list(maxfun=2e5)))
summary(aa.leaf_trichome.density.glmer)
plot(aa.leaf_trichome.density.glmer) 
overdisp_fun(aa.leaf_trichome.density.glmer) # overdispersed so I modelled an individual-level random effect

## Likelihood ratio tests
(aa.td.2 <- drop1(aa.leaf_trichome.density.glmer, test = "Chisq"))
(aa.td.1 <- drop1(update(aa.leaf_trichome.density.glmer, .~. -(Aphid.treatment + Genotype + scale(Ant.mound.dist))^2 + Aphid.treatment + Genotype + scale(Ant.mound.dist)), test = "Chisq"))
#(aa.trichome.density.anova <- anova.table(aa.leaf_trichome.density.glmer, test = "Chisq", experiment = "ant-aphid"))

with(filter(aa.trait.df, leaf_trichome.density != "NA"), interaction.plot(x.factor = Aphid.treatment, response = leaf_trichome.density, trace.factor = Genotype))
table(filter(aa.trait.df, leaf_trichome.density != "NA")$Genotype, filter(aa.trait.df, leaf_trichome.density != "NA")$Aphid.treatment)
## effects
plot(Effect(c("Genotype","Aphid.treatment"), aa.leaf_trichome.density.glmer))
plot(Effect("Genotype",aa.leaf_trichome.density.glmer)) # 30-fold variation in leaf trichome density

aa.td.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"),aa.leaf_trichome.density.glmer)) %>% mutate(response = "leaf_trichome.density")

## Calculate R2
aa.leaf_trichome.density.up <- update(aa.leaf_trichome.density.glmer, .~. -Aphid.treatment*scale(Ant.mound.dist) - scale(Ant.mound.dist)*Genotype + (1|Genotype)) 

(aa.leaf_trichome.density.R2 <- var.table(aa.leaf_trichome.density.glmer, experiment = "ant-aphid"))

## Ant-aphid: plots ----
pd <- 0.2
#cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")
p.size <- 5
l.size <- 1.5

# height
aa.h <- as.data.frame(Effect("Genotype", aa.height.lmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.sum$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size, shape = 22, fill = "black") + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + scale_y_continuous(name = "Plant height (cm)", breaks = c(20,30,40,50,60), limits = c(15,60)) + xlab("Willow genotype"); aa.h

# shoot count
#aa.sc <- as.data.frame(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.shoot.count.glmer)) %>% ggplot(aes(x = Ant.mound.dist, group = Aphid.treatment, color = Aphid.treatment)) + geom_line(aes(y = fit)) + geom_ribbon(aes(ymax = upper, ymin = lower), color = "gray", fill = "gray", alpha = 0.25) + scale_x_continuous(name = "Distance from ant mound (m)", breaks = c(1,6,12)) + ylab("No. of shoots")

#as.data.frame(Effect(c("Year","treatment"), shoot.count.glmer)) %>% ggplot(aes(x = treatment, y = fit)) + geom_point() + geom_line(aes(group = Year)) + geom_errorbar(aes(ymax = upper, ymin = lower))

# trichomes
aa.leaf_trichome.density.glmer.noG <- glmer(leaf_trichome.density ~ (Aphid.treatment + Genotype + scale(Ant.mound.dist))^2 + (1|Block) + (1|X) + (1|Block:Ant.mound.dist), data = filter(aa.trait.df, Genotype != "G"), # removed Genotype G because we didn't have data on trichomes in the aphid treatment for this plant
                                            family = "poisson",contrasts = list(Aphid.treatment = "contr.Sum", Genotype = "contr.Sum"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #
summary(aa.leaf_trichome.density.glmer.noG) # Genotypes T, S, and L are driving the interactive effect


aa.trich <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.leaf_trichome.density.glmer.noG)) %>% mutate(GxE.sig.lab = ifelse(Genotype %in% c("T","S","L"), "yes","no"), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = fit/1.1, shape = Aphid.treatment, fill = Aphid.treatment, group = Genotype)) + geom_line(aes(linetype = GxE.sig.lab)) + geom_point(size = p.size) + geom_text(aes(label = Genotype)) + ylab(expression(paste("Trichome density (no. ",cm^-2,")"))) + xlab("Aphid treatment")+ scale_shape_manual(values = c(23,21), guide = "none") + scale_fill_manual(values = c("gray50", "white"), guide = "none") + scale_linetype_manual(values =c("dotted","solid"), guide = "none"); aa.trich

aa.mechs <- plot_grid(aa.Aphis.p, aa.Fobs.p, aa.h, aa.trich, labels = "AUTO", ncol = 2, align = 'hv'); aa.mechs

save_plot("fig_2_antaphid_mechs.png", aa.mechs, base_height = 8.5, base_width = 8.5)

## Calculating direct and indirect effects of wind exposure and Genotype ----
# Always estimating this for the simplified models

# height
height.effect <- lmer(scale(Height) ~ scale(as.numeric(treatment))*scale(as.numeric(Year)) + scale(avg.moisture.vwc) + 
                        (1|genotype) + 
                        (1|Block) + 
                        (1|Block:treatment) +
                        (1|plant_ID),
                      data = w.trait.df[-height.outlier, ])#,
                      #contrasts = list(treatment = "contr.sum",
                                       #Year = "contr.sum"))
fixef(height.effect)
summary(height.effect) # SD genotype = 0.36
# direct effects: treatment = -0.28; moisture = 0.19
# wind effect on moisture = -0.41 (ref: lanphere_soil_analysis.R)
# wind indirect effect = 0.19*-0.41 = -0.08. Total treatment effect estimated as -0.35, so this is pretty darn close to a clean partition (-0.35 ~ -0.28 -0.08).

# shoot count. reran as lmer to calculate beta-coefficients 
shoot.count.effect <- lmer(scale(all.shoot.count) ~ scale(as.numeric(treatment))*scale(as.numeric(Year)) + scale(avg.moisture.vwc) + 
                             (Year|genotype) + 
                             (1|Block) + 
                             (1|Block:treatment) +
                             (1|plant_ID),
                           data = w.trait.df,
                           contrasts = list(treatment = "contr.sum",
                                            Year = "contr.sum"))
fixef(shoot.count.effect)
summary(shoot.count.effect) # genotype int SD = 0.37; genotype slope SD = 0.22
# direct effects: treatment = -0.21; moisture = 0.12
# indirect effects of wind: 0.12*-0.41 = -0.05

# shoot length
shoot.length.effect <- lmer(scale(all.shoot.avg.length) ~ scale(as.numeric(treatment)) + scale(as.numeric(Year)) + scale(avg.moisture.vwc) + 
                             (1|genotype) + 
                             (1|Block) + 
                             (1|Block:treatment) +
                             (1|plant_ID),
                           data = w.trait.df)
fixef(shoot.length.effect)
Anova(shoot.length.effect, type = 2, test = "F") # no effect of soil moisture
# direct effects: treatment = -0.21; moisture = 0.0.05
# indirect effects of wind: 0.05*-0.41 = -0.02

# leaf water content
leaf_WC.effect <- lmer(scale(leaf_WC) ~ Year + (Year|genotype) + scale(avg.moisture.vwc) + # scale(as.numeric(Year)) 
                       (1|Block) + 
                       (1|Block:treatment) +
                       (1|plant_ID),
                     data = w.trait.df[-WC_nonsense, ],
                     contrasts = list(Year = "contr.sum"))
summary(leaf_WC.effect) # genotype int SD = 0.30; genotype slope SD = 0.36
fixef(leaf_WC.effect) # -0.13 for the direct effect of soil moisture.
var.table(leaf_WC.effect, experiment = "wind") # Soil moisture explains 1.6% of the variance, absorbing the variation explained by Block and some of the residual variance.

# leaf C:N
leaf_CN.effect <- lmer(scale(leaf_C_N) ~ (1|genotype) + 
                         (1|Block) + 
                         (1|Block:treatment),
                       data = w.trait.df)
summary(leaf_CN.effect) # genotype int SD = 0.49

# root C:N. No apparent nutrient effect when the variation is restricted to biologically reasonable root_CN values.
rootCN.effect <- lmer(scale(root_CN) ~ scale(as.numeric(treatment)) + scale(nut.PC1) + (1|genotype) + (1|Block) + (1|Block:treatment),
                    data = filter(w.trait.df, root_CN < 100)) # restricting to biologically reasonable values
summary(rootCN.effect)
fixef(rootCN.effect)
Anova(rootCN.effect, type = 2, test = "F")

(rootCN.anova <- anova.table(rootCN.lmer, type = 2, experiment = "wind"))

rootCN.GxE <- as.data.frame(Effect(c("treatment","genotype"), rootCN.lmer)) %>% mutate(response = "root_CN")
write.csv(rootCN.GxE, "rootCN.means.csv")

rootCN.up <- update(rootCN.lmer, .~. -treatment*genotype + treatment + (1|genotype))

## Print Summary ----
w.trait.means <- bind_rows(height.GxE, shoot.count.GxE, shoot.length.GxE, td.GxE, WC.GxE, SLA.GxE, CN.GxE) %>% mutate(experiment = "wind")

aa.trait.means <- bind_rows(aa.height.GxE, aa.shoot.count.GxE, aa.shoot.length.GxE, aa.td.GxE, aa.WC.GxE) %>% mutate(experiment = "ant_aphid")

trait.means <- bind_rows(w.trait.means, aa.trait.means)

trait.anovas <- bind_rows(height.anova, 
                          shoot.count.anova, 
                          shoot.length.anova,
                          leaf_WC.anova,
                          SLA.anova,
                          trichome.density.anova,
                          leaf_CN.anova,
                          aa.height.anova,
                          aa.shoot.count.anova,
                          aa.mature.shoot.avg.length.anova,
                          aa.leaf_WC.anova,
                          aa.trichome.density.anova) %>%
  select(experiment, response, term, test, df, Df.res, statistic, p.value, Pr..Chisq.)

trait.R2s <- bind_rows(height.R2,
                       shoot.count.R2,
                       shoot.length.R2,
                       leaf_WC.R2,
                       SLA.R2,
                       leaf_trichome.density.R2,
                       leaf_CN.R2,
                       aa.height.R2,
                       aa.shoot.count.R2,
                       aa.mature.shoot.avg.length.R2,
                       aa.leaf_WC.R2,
                       aa.leaf_trichome.density.R2)

write.csv(trait.means, "trait.means.csv")
write.csv(trait.anovas, "trait.anovas.csv")
write.csv(trait.R2s, "trait.R2s.csv")

## testing plot for R2s ----
ggplot(trait.R2s, aes(x = response, y = var_percent, fill = Factors)) + geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() + 
  facet_grid(~experiment) 

## Old ----
## Wind: Plant architecture 2012

# explore correlations among architecture traits
arch.traits <- c("Height","all.shoot.avg.length","all.shoot.count")

# 2012 
scatterplotMatrix(filter(w.trait.df, Year == "2012")[ ,arch.traits])
height.outlier <- which(w.trait.df$Height > 80)  # not a typo according to data entry, but definitely seems to be an outlier. Maybe it was written down incorrectly in the field, because it is almost twice the size of the other plants. I've decided to remove it from this analysis.
scatterplotMatrix(filter(w.trait.df, Year == "2012")[-height.outlier,
                                                  arch.traits])
scatterplotMatrix(log(filter(w.trait.df, Year == "2012")[-height.outlier,
                                                      arch.traits])) # logging doesn't appear to linearize these relationships anymore.
corr.test(filter(w.trait.df, Year == "2012")[-height.outlier,arch.traits]) # strong correlation between total shoot length and height in 2012

# PCA: 2012
arch.pca.2012.df <- filter(w.trait.df, Year == "2012")[-height.outlier, ] %>%
  select(Block, genotype, treatment, 
         Height, all.shoot.avg.length, all.shoot.count) %>%
  na.omit()

arch.PCA.2012 <- princomp(arch.pca.2012.df[ ,arch.traits], cor = TRUE)
summary(arch.PCA.2012) # only first component has an eigenvalue > 1
plot(arch.PCA.2012)
arch.PCA.2012$loadings # positive values of PC1 indicate higher levels of plant biomass.
arch.PCA.2012$scores[ ,"Comp.1"]

arch.pca.2012.df <- mutate(arch.pca.2012.df, arch.PC1 = arch.PCA.2012$scores[ ,"Comp.1"], arch.PC2=arch.PCA.2012$scores[ ,"Comp.2"]) 

hist(arch.pca.2012.df$arch.PC1) # normally distributed
hist(arch.pca.2012.df$arch.PC2)

## Analysis

## arch PC1 2012 
with(arch.pca.2012.df, interaction.plot(treatment, genotype, arch.PC1)) # looks like main effect of wind and genotype

pc1.2012.lmer <- lmer(arch.PC1 ~ treatment*genotype +
                        (1|Block) + (1|Block:treatment),
                      data = arch.pca.2012.df,
                      contrasts = list(treatment = "contr.sum",
                                       genotype = "contr.sum"))
summary(pc1.2012.lmer)
plot(pc1.2012.lmer) # looks pretty good

#pc1.2012.mains <- update(pc1.2012.lmer, .~. -treatment:genotype)
#pc1.2012.geno <- 
# test for GxE
anova(pc1.2012.lmer,
      update(pc1.2012.lmer, .~. -treatment*genotype) + genotype),
update(pc1.2012.lmer, .~. -(treatment|genotype)))

anova(pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(pc1.2012.lmer, xvar = "treatment")

# random model
pc1.2012.rand <- lmer(arch.PC1 ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = arch.pca.2012.df)
pc1.2012.vars <- as.data.frame(VarCorr(pc1.2012.rand))

## Wind: Leaf quality 2012 

# explore correlations among architecture traits
LQ.traits <- c("leaf_WC","leaf_trichome.density")

# 2012 
scatterplotMatrix(filter(w.trait.df, Year == "2012")[ ,LQ.traits])
WC.outlier <- which(w.trait.df$leaf_WC > 5)  # I've decided to remove it from this analysis, because it is unreasonably large and was likely a measurement error
scatterplotMatrix(filter(w.trait.df, Year == "2012")[-WC.outlier,
                                                  LQ.traits])
scatterplotMatrix(log(filter(w.trait.df, Year == "2012")[-WC.outlier,
                                                      LQ.traits]+1)) # logging doesn't appear to linearize these relationships anymore.
corr.test(filter(w.trait.df, Year == "2012")[-WC.outlier,LQ.traits]) # positive correlation between leaf_WC and trichome density

# PCA: 2012
LQ.pca.2012.df <- filter(w.trait.df, Year == "2012")[-WC.outlier, ] %>%
  select(Block, genotype, treatment, 
         leaf_WC, leaf_trichome.density) %>%
  mutate(log_trichome.density = log(leaf_trichome.density+1),
         log_WC = log(leaf_WC)) %>%
  na.omit()

log.LQ.traits <- c("log_WC","log_trichome.density")

LQ.PCA.2012 <- princomp(LQ.pca.2012.df[ ,log.LQ.traits], cor = TRUE)
summary(LQ.PCA.2012) # only first component has an eigenvalue > 1
plot(LQ.PCA.2012)
LQ.PCA.2012$loadings # positive values of PC1 indicate higher levels of leaf water content and trichome densities.
LQ.PCA.2012$scores[ ,"Comp.1"]

LQ.pca.2012.df <- mutate(LQ.pca.2012.df, LQ.PC1 = LQ.PCA.2012$scores[ ,"Comp.1"]) 

hist(LQ.pca.2012.df$LQ.PC1) # much better distribution with log transformations

## Analysis

# LQ PC1 
LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ treatment + (treatment|genotype) +
                           (1|Block) + (1|Block:treatment),
                         data = LQ.pca.2012.df)
summary(LQ.pc1.2012.lmer)
plot(LQ.pc1.2012.lmer) # looks okay

# test for GxE
anova(LQ.pc1.2012.lmer,
      update(LQ.pc1.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(LQ.pc1.2012.lmer, .~. -(treatment|genotype)))

anova(LQ.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc1.2012.lmer, xvar = "treatment")

# random model
LQ.pc1.2012.rand <- lmer(LQ.PC1 ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = LQ.pca.2012.df)
summary(LQ.pc1.2012.rand )
LQ.pc1.2012.vars <- as.data.frame(VarCorr(LQ.pc1.2012.rand))

# WC

## Wind plant architecture 2013

# 2013 
scatterplotMatrix(filter(w.trait.df, Year == "2013")[ ,arch.traits]) # no outliers
scatterplotMatrix(log(filter(w.trait.df, Year == "2013")[ ,arch.traits])) # logging appears to help a bit
corr.test(filter(w.trait.df, Year == "2013")[ ,arch.traits]) # all traits are now highly correlated
corr.test(log(filter(w.trait.df, Year == "2013")[ ,arch.traits]))

# PCA: 2013
arch.pca.2013.df <- w.trait.df %>%
  filter(Year == "2013") %>%
  select(Block, genotype, treatment, 
         Height, all.shoot.avg.length, all.shoot.count) %>%
  na.omit()

arch.PCA.2013 <- princomp(log(arch.pca.2013.df[ ,arch.traits]), cor = TRUE)
summary(arch.PCA.2013) # only first component has an eigenvalue > 1
plot(arch.PCA.2013)
arch.PCA.2013$loadings # positive values of PC1 indicate higher levels of productivity.
arch.PCA.2013$scores[ ,"Comp.1"]

arch.pca.2013.df <- mutate(arch.pca.2013.df, arch.PC1 = -1*arch.PCA.2013$scores[ ,"Comp.1"]) # multiplied by -1 so that positive values of PC1 indicate higher levels of productivity

hist(arch.pca.2013.df$arch.PC1) # normally distributed with log-transformation

## Analysis

# arch PC1 2013 
with(arch.pca.2013.df, interaction.plot(treatment, genotype, arch.PC1))

pc1.2013.lmer <- lmer(arch.PC1 ~ treatment + (treatment|genotype) +
                        (1|Block) + (1|Block:treatment),
                      data = arch.pca.2013.df)
summary(pc1.2013.lmer)
plot(pc1.2013.lmer) # looks pretty good

# test for GxE
anova(pc1.2013.lmer,
      update(pc1.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(pc1.2013.lmer, .~. -(treatment|genotype)))

anova(pc1.2013.lmer, ddf = "Kenward-Roger")
visreg(pc1.2013.lmer, xvar = "treatment")

# random model
pc1.2013.rand <- lmer(arch.PC1 ~  (1|treatment) + (1|genotype) + (1|Block) + (1|Block:treatment), data = arch.pca.2013.df)
summary(pc1.2013.rand)
0.5305/(0.3961+0.5186+0.2041+0.5305+1.0331) # 20% of the variance for wind exposure
0.2041/(0.3961+0.5186+0.2041+0.5305+1.0331) # 8% of the variance for genotype.
pc1.2013.vars <- as.data.frame(VarCorr(pc1.2013.rand))

# plant height
with(arch.pca.2013.df, interaction.plot(treatment, genotype, Height))
height.2013.lmer <- lmer(Height ~ treatment + (treatment|genotype) +
                           (1|Block) + (1|Block:treatment),
                         data = filter(w.trait.df, Year == "2013"))
summary(height.2013.lmer)
summary(update(height.2013.lmer, .~. -1))
1-(11.648/22.62)
coef(update(height.2013.lmer, .~. -1))

plot(height.2013.lmer) # looks okay, note that logging qualitatively alters the conclusion about GxE.

# test for GxE
anova(height.2013.lmer,
      update(height.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(height.2013.lmer, .~. -(treatment|genotype)))

anova(height.2013.lmer, ddf = "Kenward-Roger")
visreg(height.2013.lmer, xvar = "treatment")

# random model
height.2013.rand <- lmer(Height ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = arch.pca.2013.df)
summary(height.2013.rand)
56.295/(19.083+16.579+1.176+23.105+56.295+64.671) # wind exposure explained 31% of the variance in plant height
1.176/(19.083+16.579+1.176+23.105+56.295+64.671) # genotype explained <1% of the variance in plant height
23.105/(19.083+16.579+1.176+23.105+56.295+64.671) # GxE explained 13% of the variance in plant height.
height.2013.vars <- as.data.frame(VarCorr(height.2013.rand))

# shoot count
hist(arch.pca.2013.df$all.shoot.count)
hist(sqrt(arch.pca.2013.df$all.shoot.count))
hist(log(arch.pca.2013.df$all.shoot.count))

shoot.count.2013.lmer <- lmer(all.shoot.count ~ treatment + (treatment|genotype) +
                                (1|Block) + (1|Block:treatment),
                              data = arch.pca.2013.df)
summary(shoot.count.2013.lmer)
summary(update(shoot.count.2013.lmer, .~. -1))
1-(3.4655/5.7586)
plot(shoot.count.2013.lmer)  

# test for GxE
anova(shoot.count.2013.lmer,
      update(shoot.count.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(shoot.count.2013.lmer, .~. -(treatment|genotype))) # marginally close to GxE (should divide p-value by 2)

anova(shoot.count.2013.lmer, ddf = "Kenward-Roger")
visreg(shoot.count.2013.lmer, xvar = "treatment")

# random model
shoot.count.2013.rand <- lmer(all.shoot.count ~  (1|treatment) + (1|genotype) + (1|Block) + (1|Block:treatment), data = arch.pca.2013.df)
summary(shoot.count.2013.rand)
2.5598/(1.7568+1.3476+0.8464+2.5598+7.9179) # wind exposure explained 18% of the variance
0.8464/(1.7568+1.3476+0.8464+2.5598+7.9179) # genotype explained 6% of the variance
shoot.count.2013.vars <- as.data.frame(VarCorr(shoot.count.2013.rand))

# average shoot length
hist(filter(w.trait.df, Year == "2013")$all.shoot.avg.length)
avg.shoot.length.2013.lmer <- lmer(log(all.shoot.avg.length) ~
                                     treatment + (treatment|genotype) +
                                     (1|Block) + (1|Block:treatment),
                                   data = filter(w.trait.df, Year == "2013"))
summary(avg.shoot.length.2013.lmer)
summary(lmer(all.shoot.avg.length ~
               treatment + (treatment|genotype) - 1 +
               (1|Block) + (1|Block:treatment),
             data = filter(w.trait.df, Year == "2013")))
1-(0.6018/0.8429) # 29% shorter in wind exposed plots
coef(lmer(all.shoot.avg.length ~
            treatment + (treatment|genotype) +
            (1|Block) + (1|Block:treatment),
          data = filter(w.trait.df, Year == "2013")))
0.8221796/0.3622077
plot(avg.shoot.length.2013.lmer) # much better with log-transformation

# test for GxE
anova(avg.shoot.length.2013.lmer,
      update(avg.shoot.length.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(avg.shoot.length.2013.lmer, .~. -(treatment|genotype)))

anova(avg.shoot.length.2013.lmer, ddf = "Kenward-Roger")
visreg(avg.shoot.length.2013.lmer, xvar = "treatment")

# random model
avg.shoot.length.2013.rand <- lmer(log(all.shoot.avg.length) ~  (1|treatment) + (1|genotype) + (1|Block) + (1|Block:treatment), data = filter(w.trait.df, Year == "2013"))
summary(avg.shoot.length.2013.rand)
0.06784/(0.09846+0.29431+0.06784+0.04201+0.29892) # genotype explained 8% of the total variance
0.04201/(0.09846+0.29431+0.06784+0.04201+0.29892) # wind exposure explained 5% of the total variance
avg.shoot.length.2013.vars <- as.data.frame(VarCorr(avg.shoot.length.2013.rand))

# shoot length
shoot.length.2013.lmer <- lmer(log(all.shoot.total.length) ~ treatment + (treatment|genotype) +
                                 (1|Block) + (1|Block:treatment),
                               data = arch.pca.2013.df)
summary(shoot.length.2013.lmer)
plot(shoot.length.2013.lmer) # much better with log-transformation

# test for GxE
anova(shoot.length.2013.lmer,
      update(shoot.length.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(shoot.length.2013.lmer, .~. -(treatment|genotype)))

anova(shoot.length.2013.lmer, ddf = "Kenward-Roger")
visreg(shoot.length.2013.lmer, xvar = "treatment")

# random model
shoot.length.2013.rand <- lmer(all.shoot.total.length ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = arch.pca.2013.df)
summary(shoot.length.2013.rand)
shoot.length.2013.vars <- as.data.frame(VarCorr(shoot.length.2013.rand))

mean(c(31,18,5)) # 18% of variance for wind exposure
mean(c(8,6,13)) # 9% of variance for G and/or GxE


## Wind: Leaf quality traits 2013 
LQ.traits.2013 <- c("leaf_WC","SLA","leaf_C_N","leaf_C","leaf_N",
                    "larva.wet.wt.exp1", "larva.wet.wt.exp2")

# explore correlations
scatterplotMatrix(filter(w.trait.df, Year == "2013")[ ,LQ.traits.2013]) # non-linearity between C:N and N.
WC_nonsense <- which(filter(w.trait.df, Year == "2013")[ ,"leaf_WC"] < 0) # nonsensical values for leaf WC. Essentially, they suggest that the wet leaf mass is less than the dry leaf mass, which is not possible. Therefore, I'm removing them from further analysis
scatterplotMatrix(filter(w.trait.df, Year == "2013")[-WC_nonsense,LQ.traits.2013])
scatterplotMatrix(log(filter(w.trait.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]+1))
corr.test(filter(w.trait.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]) # larva.wet.wt.exp1 has the lowest sample size, but has a strong positive correlation with larva.wet.wt.exp2; therefore, I'm only going to use the latter for the PCA analysis.
corr.test(log(filter(w.trait.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]+1))

# PCA
LQ.pca.2013.df <- filter(w.trait.df, Year == "2013")[-WC_nonsense, ] %>%
  select(Block, genotype, treatment, 
         leaf_WC, SLA, leaf_C_N, larva.wet.wt.exp2) %>%
  mutate(log_CN = log(leaf_C_N),
         log_SLA = log(SLA),
         log_WC = log(leaf_WC),
         log_larva.wet = log(larva.wet.wt.exp2 + 1)) %>%
  na.omit()

pca.LQ.traits.2013 <- c("leaf_WC", "SLA", "leaf_C_N", "larva.wet.wt.exp2")
pca.log.LQ.traits.2013 <- c("log_WC", "SLA", "log_CN", "log_larva.wet")

scatterplotMatrix(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013])
corr.test(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013])

LQ.PCA.2013 <- princomp(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013], cor = TRUE)
summary(LQ.PCA.2013) # first two components have an eigenvalue > 1
plot(LQ.PCA.2013)
LQ.PCA.2013$loadings # Positive values of PC1 have higher leaf water content, whereas lower values have higher SLA and C:N. Positive value of PC2 have higher Spodoptera growth rates, whereas lower values have high C:N.
LQ.PCA.2013$scores[ ,c("Comp.1","Comp.2")]

LQ.pca.2013.df <- mutate(LQ.pca.2013.df, 
                         LQ.PC1 = LQ.PCA.2013$scores[ ,"Comp.1"],
                         LQ.PC2 = LQ.PCA.2013$scores[ ,"Comp.2"]) 

hist(LQ.pca.2013.df$LQ.PC1) # normally distributed
hist(LQ.pca.2013.df$LQ.PC2) # roughly normal

## Analysis
# note that I use the filtered w.trait.df for the individual traits to preserve the original sample sizes.

# LQ PC1
LQ.pc1.2013.lmer <- lmer(LQ.PC1 ~ treatment + (treatment|genotype) +
                           (1|Block) + (1|Block:treatment),
                         data = LQ.pca.2013.df)
summary(LQ.pc1.2013.lmer)
plot(LQ.pc1.2013.lmer) # looks okay

# test for GxE
anova(LQ.pc1.2013.lmer,
      update(LQ.pc1.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(LQ.pc1.2013.lmer, .~. -(treatment|genotype)))

anova(LQ.pc1.2013.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc1.2013.lmer, xvar = "treatment")

# random model
LQ.pc1.2013.rand <- lmer(LQ.PC1 ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = LQ.pca.2013.df)
summary(LQ.pc1.2013.rand )
LQ.pc1.2013.vars <- as.data.frame(VarCorr(LQ.pc1.2013.rand))

# WC
leaf_WC.2013.lmer <- lmer(log(leaf_WC) ~ treatment + (treatment|genotype) +
                            (1|Block) + (1|Block:treatment),
                          data = filter(w.trait.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_WC.2013.lmer)
plot(leaf_WC.2013.lmer)  # better with log-transformation

# test for GxE
anova(leaf_WC.2013.lmer,
      update(leaf_WC.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(leaf_WC.2013.lmer, .~. -(treatment|genotype)))

anova(leaf_WC.2013.lmer, ddf = "Kenward-Roger")
visreg(leaf_WC.2013.lmer, xvar = "treatment")

# random model
leaf_WC.2013.rand <- lmer(log(leaf_WC) ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = filter(w.trait.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_WC.2013.rand )
leaf_WC.2013.vars <- as.data.frame(VarCorr(leaf_WC.2013.rand))

# LQ PC2
LQ.pc2.2013.lmer <- lmer(LQ.PC2 ~ treatment + (treatment|genotype) +
                           (1|Block) + (1|Block:treatment),
                         data = LQ.pca.2013.df)
summary(LQ.pc2.2013.lmer)
plot(LQ.pc2.2013.lmer) # looks okay

# test for GxE
anova(LQ.pc2.2013.lmer,
      update(LQ.pc2.2013.lmer, .~. -(treatment|genotype) + (1|genotype)), # Genotype effect would be significant when dividing p-value by 2, since this test is on the boundary.
      update(LQ.pc2.2013.lmer, .~. -(treatment|genotype)))

anova(LQ.pc2.2013.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc2.2013.lmer, xvar = "treatment")

# random model
LQ.pc2.2013.rand <- lmer(LQ.PC2 ~  (1|treatment) + (treatment|genotype) + (1|Block) + (1|Block:treatment), data = LQ.pca.2013.df)
summary(LQ.pc2.2013.rand )
LQ.pc2.2013.vars <- as.data.frame(VarCorr(LQ.pc2.2013.rand))

## Ant-aphid 2012: architecture 
aa.arch.traits <- c("Height","all.shoot.count","mature.shoot.total.length")

# explore correlations
scatterplotMatrix(filter(aa.df, Year == "2012")[ ,aa.arch.traits])

# PCA: 2012
aa.arch.pca.2012.df <- filter(aa.df, Year == "2012") %>%
  select(Block, Genotype, Aphid.Treatment, Ant.Mound.Dist, fact.Ant.Mound.Dist,
         Height, mature.shoot.total.length, all.shoot.count) %>%
  na.omit()

aa.arch.PCA.2012 <- princomp(aa.arch.pca.2012.df[ ,aa.arch.traits], cor = TRUE)
summary(aa.arch.PCA.2012) # only first component has an eigenvalue > 1
plot(aa.arch.PCA.2012)
aa.arch.PCA.2012$loadings 
aa.arch.PCA.2012$scores[ ,"Comp.1"]

aa.arch.pca.2012.df <- mutate(aa.arch.pca.2012.df, arch.PC1 = -1*aa.arch.PCA.2012$scores[ ,"Comp.1"]) # multiplied by -1 so that positive values of PC1 indicate higher levels of productivity.

hist(aa.arch.pca.2012.df$arch.PC1) # roughly normal distribution

## Analysis

# arch PC1 2012 
aa.pc1.2012.lmer <- lmer(arch.PC1 ~ Aphid.Treatment*Ant.Mound.Dist + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data = aa.arch.pca.2012.df)
plot(aa.pc1.2012.lmer) # looks pretty good
summary(aa.pc1.2012.lmer)

# test for GxE
anova(aa.pc1.2012.lmer,
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)),
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)),
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

anova(aa.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.pc1.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.pc1.2012.rand <- lmer(arch.PC1 ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data = aa.arch.pca.2012.df)
summary(aa.pc1.2012.rand)
aa.pc1.2012.vars <- as.data.frame(VarCorr(aa.pc1.2012.rand))


## Ant-aphid 2012: leaf quality 

# explore correlations
scatterplotMatrix(filter(aa.df, Year == "2012")[ ,LQ.traits])
corr.test(filter(aa.df, Year == "2012")[ ,LQ.traits])

# PCA
aa.LQ.pca.2012.df <- filter(aa.df, Year == "2012") %>%
  select(Block, Genotype, Aphid.Treatment, 
         Ant.Mound.Dist, fact.Ant.Mound.Dist, 
         leaf_WC, leaf_trichome.density) %>%
  mutate(log_WC = log(leaf_WC),
         log_trichomes = log(leaf_trichome.density + 1)) %>%
  na.omit()

# decided from correlation plots that the non-logged transformed relationships appear to be more linear, so I've decided to go with those ones.
aa.pca.log.LQ.traits <- c("log_WC", "log_trichomes")
scatterplotMatrix(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])
corr.test(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])

aa.LQ.PCA <- princomp(aa.LQ.pca.2012.df[ ,c("leaf_WC","leaf_trichome.density")], cor = TRUE) 
summary(aa.LQ.PCA) # only first component has an eigenvalue > 1
plot(aa.LQ.PCA)
aa.LQ.PCA$loadings # PC1 represents the correlation between leaf water content and trichome density.
aa.LQ.PCA$scores[ ,"Comp.1"]

aa.LQ.PCA.df <- mutate(aa.LQ.pca.2012.df, 
                       LQ.PC1 = -1*aa.LQ.PCA$scores[ ,"Comp.1"]) # multiply by -1 so positive values indicate higher leaf water content and higher trichome density

hist(aa.LQ.PCA.df$LQ.PC1) # a bit skewed, but transformation should help
hist(log(aa.LQ.PCA.df$LQ.PC1+2))

## Analysis

with(aa.LQ.PCA.df, table(Genotype, Ant.Mound.Dist, Aphid.Treatment)) # testing 3-way interaction is unwise
# 2-way interactions mostly unbalanced
with(aa.LQ.PCA.df, table(Genotype, Ant.Mound.Dist)) 
with(aa.LQ.PCA.df, table(Genotype, Aphid.Treatment)) 
with(aa.LQ.PCA.df, table(Aphid.Treatment, Ant.Mound.Dist)) # this one isn't too bad.
# 1-way effects
with(aa.LQ.PCA.df, table(Genotype)) 

contr.sum(aa.LQ.PCA.df$Genotype)

## Analysis

# LQ PC1. Ran without 3-way interaction because there were convergence problems
aa.LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ Aphid.Treatment*Ant.Mound.Dist + 
                              (Aphid.Treatment + Ant.Mound.Dist|Genotype) +
                              (1|Block) + (1|Block:Ant.Mound.Dist),
                            data =  aa.LQ.PCA.df)
plot(aa.LQ.pc1.2012.lmer) # looks pretty good
summary(aa.LQ.pc1.2012.lmer)

# test for GxE
anova(aa.LQ.pc1.2012.lmer,
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype) + (Ant.Mound.Dist|Genotype)), # replaced with Ant.Mound.Dist because it had the highest variance of the possible interactions
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype) + (1|Genotype)), # likely significance of Genotype after dividing p-value by 2 to account for boundary effect.
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype)))

anova(aa.LQ.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.LQ.pc1.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.LQ.pc1.2012.rand <- lmer(LQ.PC1 ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                              (Aphid.Treatment + Ant.Mound.Dist|Genotype) +
                              (1|Block) + (1|Block:Ant.Mound.Dist),
                            data =  aa.LQ.PCA.df)
summary(aa.LQ.pc1.2012.rand)

aa.LQ.pc1.2012.vars <- as.data.frame(VarCorr(aa.LQ.pc1.2012.rand))

## Survival analysis ----

# explore 2-way interactions
with(w.trait.df, table(Dead, genotype, treatment))
with(w.trait.df, table(Dead, genotype, Year))
with(w.trait.df, table(Dead, treatment, Year))

# explore main effects
with(w.trait.df, table(Dead, treatment))
with(w.trait.df, table(Dead, genotype))
with(w.trait.df, table(Dead, Year))

# GLMM. Modelled genotype as random effect because the model was having difficulty converging.
surv.glmer <- glmer(Dead ~ treatment + Year + genotype +
                      (1|Block) + 
                      (1|Block:treatment) +
                      (1|plant_ID),
                    data = w.trait.df,
                    contrasts = list(treatment = "contr.sum",
                                     genotype = "contr.sum",
                                     Year = "contr.sum"),
                    family = "binomial",
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=2e5)))
summary(surv.glmer)
# overdisp_fun() # overdispersion test unnecessary, since the response is binary.
plot(surv.glmer) 

# Likelihood ratio tests for main effects
drop1(surv.glmer, test = "Chisq") # no ExE
drop1(update(surv.glmer, .~. -treatment:Year), test = "Chisq") # only Year effect

# no G or GxE effect
anova(surv.glmer, 
      update(surv.glmer, .~. -(treatment|genotype) + (1|genotype)),
      update(surv.glmer, .~. -(treatment|genotype)))

## Wind: larva wet weight exp. 1 analysis ----
hist(w.trait.df$larva.wet.wt.exp1)

with(filter(w.trait.df, larva.wet.wt.exp1 != "NA"),
     interaction.plot(treatment, genotype, larva.wet.wt.exp1))

# convergence problem if I included GxE, so I only modelled the main effects
larva.wet.wt.exp1.glmer <- glmer(
  larva.wet.wt.exp1 ~ treatment + genotype +
    (1|Block) + (1|Block:treatment),
  data = w.trait.df, family = "poisson",
  contrasts = list(treatment = "contr.sum",
                   genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(larva.wet.wt.exp1.glmer)
overdisp_fun(larva.wet.wt.exp1.glmer) # no overdispersion
plot(larva.wet.wt.exp1.glmer) # weird residuals...

## Likelihood ratio tests
drop1(larva.wet.wt.exp1.glmer, test = "Chisq") # no G or E

## Wind: larva wet weight exp. 2 analysis ----
hist(w.trait.df$larva.wet.wt.exp2)

with(filter(w.trait.df, larva.wet.wt.exp2 != "NA"),
     interaction.plot(treatment, genotype, larva.wet.wt.exp2))

larva.wet.wt.exp2.glmer <- glmer(
  larva.wet.wt.exp2 ~ treatment*genotype +
    (1|Block) + (1|Block:treatment),
  data = w.trait.df, family = "poisson",
  contrasts = list(treatment = "contr.sum",
                   genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(larva.wet.wt.exp2.glmer)
overdisp_fun(larva.wet.wt.exp2.glmer) # no overdispersion
plot(larva.wet.wt.exp2.glmer) # weird residuals...

## Likelihood ratio tests
drop1(larva.wet.wt.exp2.glmer, test = "Chisq") # no GxE
drop1(update(larva.wet.wt.exp2.glmer, .~. -treatment:genotype), 
      test = "Chisq") # no genotype or treatment effect

#### Examine correlation between soil properties and plant traits ---- 
## Wind experiment

# get data at plot level to avoid pseudoreplication
trait.soil.plots <- w.trait.df %>%
  select(Block, treatment, Plot_code, Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC, SLA, leaf_C_N, root_CN, Total.N, nut.PC1, percent.TOM, avg.moisture.vwc) %>%
  group_by(Block, treatment, Plot_code) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>%
  ungroup() 

# examine correlations. Soil moisture appears to affect plant height, number of shoots, and leaf water content.
corr.test(x = select(trait.soil.plots, Total.N:avg.moisture.vwc),
          y = select(trait.soil.plots, Height:leaf_C_N, root_CN), adjust = "none")
with(trait.soil.plots, cor.test(avg.moisture.vwc, Height))
with(trait.soil.plots, cor.test(avg.moisture.vwc, all.shoot.count))
with(trait.soil.plots, cor.test(avg.moisture.vwc, leaf_WC))
with(trait.soil.plots, cor.test(nut.PC1, root_CN))

plot(Height ~ avg.moisture.vwc, trait.soil.plots)
plot(all.shoot.count ~ avg.moisture.vwc, trait.soil.plots)
plot(leaf_WC ~ avg.moisture.vwc, trait.soil.plots)
plot(root_CN ~ nut.PC1, trait.soil.plots) # quite variable still...
summary(lm(scale(root_CN) ~ scale(nut.PC1), trait.soil.plots)) # no longer significant when I restrict root_CN to values < 100.
summary(lm(scale(leaf_WC) ~ scale(avg.moisture.vwc), trait.soil.plots))

# ant-aphid experiment
aa.trait.soil.plots <- aa.df %>%
  select(Block, Ant.Mound.Dist, plot.id, Height, mature.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC, moisture.vwc.July9_13) %>%
  group_by(Block, Ant.Mound.Dist, plot.id) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>%
  ungroup()

# no signficant correlations. Likely due to the small range in soil moisture.
corr.test(x = select(aa.trait.soil.plots, moisture.vwc.July9_13),
          y = select(aa.trait.soil.plots, Height:leaf_WC), adjust = "none")