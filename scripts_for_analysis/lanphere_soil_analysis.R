############################################
## Description: This script analyzes soil data of Lanphere experiments (mostly wind).
## Code author: Matt Barbour
## Email: matthew.a.barbour@gmail.com
############################################

#### load libraries ----
source('scripts_for_analysis/required_libraries.R')

#### upload data ----
w.soil <- read.csv('final_data/wind_soil_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) # necessary for modelling as a random effect model

aa.soil <- read.csv('final_data/ant_aphid_soil_2013.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) # necessary for modelling as a random effect


#### Wind: correlation among soil properties ----
props <- colnames(select(w.soil, Total.N, nut.PC1, percent.TOM, avg.moisture.vwc))
scatterplotMatrix(w.soil[ ,props])
corr.test(w.soil[ ,props]) # the only soil properties that are not highly correlated are Total.N and avg.moisture.vwc
plot(log(avg.moisture.vwc) ~ log(percent.TOM), w.soil)

#### Wind: Total Nitrogen analysis ----

N.lmer <- lmer(Total.N ~ Wind.Exposure + (1|Block), w.soil)
summary(N.lmer)

Anova(N.lmer, test = "F") # marginally significant effect of wind
var.table(N.lmer, experiment = "ant-aphid") # wind explains 17% of the variance

(N.ef <- Effect("Wind.Exposure", N.lmer)) # 1.4-fold higher N content in exposed plots.
plot(N.ef)

# calculate the standardized beta-coefficient of wind exposure on soil N. This is useful for calculating the direct and indirect effects of N on plant traits.
N.std <- lmer(scale(Total.N) ~ scale(as.numeric(Wind.Exposure)) + (1|Block), w.soil)
summary(N.std) # beta = -0.42

#### Wind: Nutrient PC1 analysis ----
nutPC1.lmer <- lmer(nut.PC1 ~ Wind.Exposure + (1|Block), w.soil)
summary(nutPC1.lmer)

Anova(nutPC1.lmer, test = "F") # no effect of Wind.Exposure

##### Wind: Percent organic matter analysis ----
TOM.lmer <- lmer(percent.TOM ~ Wind.Exposure + (1|Block), w.soil)
summary(TOM.lmer)

Anova(TOM.lmer, test = "F") # no effect of wind

#### Wind: Soil moisture content analysis ----
moist.lmer <- lmer(log(avg.moisture.vwc) ~ Wind.Exposure + (1|Block), w.soil)
summary(moist.lmer)

Anova(moist.lmer, test = "F") # marginal effect of wind

(moist.ef <- Effect("Wind.Exposure", moist.lmer, transformation = list(link = log, inverse = exp))) # 1.4-fold more soil moisture in unexposed plots
plot(moist.ef)

# calculate the standardized beta-coefficient of wind exposure on soil moisture. This is useful for calculating the direct and indirect effects of wind on plant traits.
moist.std <- lmer(scale(avg.moisture.vwc) ~ scale(as.numeric(Wind.Exposure)) + (1|Block), w.soil)
summary(moist.std) # positive effect of unexposed plots: beta = 0.41

#### Wind: Soil Plots ----
w.moist.p <- as.data.frame(Effect("Wind.Exposure", moist.lmer, transformation = list(link = log, inverse = exp))) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + ylab(expression(paste("Soil moisture (",m^3, m^-3,")"))) + xlab("Wind Exposure"); w.moist.p

w.N.p <- as.data.frame(Effect("Wind.Exposure", N.lmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w*2) + ylab(expression(paste("Soil N (",mu,"g","/10", cm^2, "/11 days)"  ))) + xlab("Wind Exposure"); w.N.p
