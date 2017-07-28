############################################
## Description: This script analyzes soil data of Lanphere experiments (mostly wind).
## Code author: Matt Barbour
## Email: barbour@zoology.ubc.ca
############################################

#### load libraries ----
#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must be before dplyr
library(dplyr)
library(lme4) # for mixed-effect models
library(car) # test significance of fixed effects
library(RLRsim) # test significance of random effects
library(psych) # correlation analysis
library(effects) # calculating effects

library(RCurl) # for loading github files directly.

script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/model_diagnostic_functions.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#### upload data ----
w.soil <- read.csv('final_data/wind_soil_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) # necessary for modelling as a random effect model

aa.soil <- read.csv('final_data/ant_aphid_soil_2013.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) # necessary for modelling as a random effect

#### Ant-aphid: soil moisture analysis ----
aa.soil.lmer <- lmer(moisture.vwc.July9_13 ~ Ant.mound.dist + (1|Block), data = aa.soil)
summary(aa.soil.lmer)
Anova(aa.soil.lmer, test = "F") # no effect of ant mound distance on soil
var.table(aa.soil.lmer, experiment = "ant-aphid")
exactRLRT(aa.soil.lmer) # no effect of block either

#### Wind: correlation among soil properties ----

# Soil moisture of plots is highly correlated between years, so I feel justified simply analyzing the average between these years (avg.moisture.vwc).
plot(avg.moisture.vwc.2013 ~ avg.moisture.vwc.2012, w.soil)
with(w.soil, cor.test(avg.moisture.vwc.2013, avg.moisture.vwc.2012)) 
props <- colnames(select(w.soil, Total.N, nut.PC1, percent.TOM, avg.moisture.vwc))
scatterplotMatrix(w.soil[ ,props])
corr.test(w.soil[ ,props]) # the only soil properties that are not highly correlated are Total.N and avg.moisture.vwc
plot(log(avg.moisture.vwc) ~ log(percent.TOM), w.soil)

#### Wind: Soil properties PCA ----
props.all <- colnames(select(w.soil, NO3.N:percent.TOM, avg.moisture.vwc))
library(missMDA)

soil.mat <- as.matrix(w.soil[ ,props.all])

soil.nb <- estim_ncpPCA(soil.mat, ncp.max = 5) # 2 PCs
w.soil.imp <- imputePCA(soil.mat, ncp = 2, scale = TRUE)
 
props.PCA <- princomp(w.soil.imp$completeObs, cor = TRUE)
biplot(props.PCA)
summary(props.PCA) # 39% of variance explained by PC1 (14% for PC2)
props.PCA$loadings # positive values of PC1 indicate higher supply rates of micronutrients, but lower supply of nitrogen compounds
hist(props.PCA$scores[ ,"Comp.1"]) # look at distribution

props.lmer <- lmer(props.PCA$scores[ ,"Comp.2"] ~ Wind.Exposure + (1|Block), w.soil)
summary(props.lmer)
plot(props.lmer)

Anova(props.lmer, test = "F") 

with(filter(w.soil, percent.TOM > 0), interaction.plot(x.factor = Wind.Exposure, trace.factor = Block, response = props.PCA$scores[ ,"Comp.1"]))

#### Wind: Total Nitrogen analysis ----

N.lmer <- lmer(Total.N ~ Wind.Exposure + (1|Block), w.soil)
summary(N.lmer)
plot(N.lmer)

Anova(N.lmer, test = "F") # marginally significant effect of wind
var.table(N.lmer, experiment = "ant-aphid") # wind explains 17% of the variance
exactRLRT(N.lmer) # no significant effect of block

(N.ef <- Effect("Wind.Exposure", N.lmer)) # 1.4-fold higher N content in exposed plots.
plot(N.ef)

# calculate the standardized beta-coefficient of wind exposure on soil N. This is useful for calculating the direct and indirect effects of N on plant traits.
N.std <- lmer(scale(Total.N) ~ scale(as.numeric(Wind.Exposure)) + (1|Block), w.soil)
summary(N.std) # beta = -0.42

#### Wind: Nutrient PC1 analysis ----
nutPC1.lmer <- lmer(nut.PC1 ~ Wind.Exposure + (1|Block), w.soil)
summary(nutPC1.lmer)
plot(nutPC1.lmer)

Anova(nutPC1.lmer, test = "F") # no effect of Wind.Exposure
exactRLRT(nutPC1.lmer) # no significant effect of block

##### Wind: Percent organic matter analysis ----
TOM.lmer <- lmer(percent.TOM ~ Wind.Exposure + (1|Block), w.soil)
summary(TOM.lmer)
plot(TOM.lmer)

Anova(TOM.lmer, test = "F") # no effect of wind
exactRLRT(TOM.lmer) # no significant effect of block

#### Wind: Soil moisture content analysis ----
moist.lmer <- lmer(log(avg.moisture.vwc) ~ Wind.Exposure + (1|Block), w.soil)
summary(moist.lmer)
plot(moist.lmer) # a bit better with log-transformation

Anova(moist.lmer, test = "F") # marginal effect of wind
exactRLRT(moist.lmer) # no significant effect of block

(moist.ef <- Effect("Wind.Exposure", moist.lmer, transformation = list(link = log, inverse = exp))) # 1.4-fold more soil moisture in unexposed plots
plot(moist.ef)

# calculate the standardized beta-coefficient of wind exposure on soil moisture. This is useful for calculating the direct and indirect effects of wind on plant traits.
moist.std <- lmer(scale(avg.moisture.vwc) ~ scale(as.numeric(Wind.Exposure)) + (1|Block), w.soil)
summary(moist.std) # positive effect of unexposed plots: beta = 0.41

#### Wind: Soil Plots
w.moist.p <- as.data.frame(Effect("Wind.Exposure", moist.lmer, transformation = list(link = log, inverse = exp))) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + ylab(expression(paste("Soil moisture (",m^3, m^-3,")"))) + xlab("Wind Exposure"); w.moist.p

w.N.p <- as.data.frame(Effect("Wind.Exposure", N.lmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + ylab(expression(paste("Soil N (",mu,"g","/10", cm^2, "/11 days)"  ))) + xlab("Wind Exposure"); w.N.p
