## load packages
#library("reshape")
library(dplyr)
library(tidyr)
library(lme4)
library(car)

## upload data

# total organic matter
windTOM <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Wind - Soil Total Organic Matter.csv") %>%
  tbl_df() %>%
  filter(BAD.Sample == "n") %>%
  mutate(OvenWt = Cruc....Oven.Dry.Wt. - Crucible.Wt.,
         IgniteWt = Cruc....Ignited.Wt. - Crucible.Wt.,
         PercTOM = (OvenWt - IgniteWt)/OvenWt,
         Block = as.factor(Wind.Site.No.)) %>%
  select(Block, Wind.Treatment, PercTOM)

# Percent Total Organic Matter doesn't significantly vary among sites, although there is a trend for greater amount in unexposed treatments, which makes sense.
PercTOM.lmer <- lmer(PercTOM ~ Wind.Treatment + (1|Block), data = windTOM)
summary(PercTOM.lmer)
Anova(PercTOM.lmer, test.statistic = "F")
plot(PercTOM.lmer)

# still no effect in a linear model.
PercTOM.lm <- lm(PercTOM ~ Wind.Treatment, data = windTOM)
summary(PercTOM.lm)
anova(PercTOM.lm)


## soil moisture, temperature, and electrical conductivity for Sep 18
wind.sep.18.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Lanphere Dunes - Soil Moisture-Temp-EC - Sep 18 2012.csv") %>%
  tbl_df() %>%
  filter(Treatment %in% c("unexposed","exposed")) %>%
  select(Measurement.Time:EC.2) %>%
  group_by(Site, Treatment) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>%
  mutate(Moisture = (moisture + moisture.1 + moisture.2)/3,
         Temp = (temp + temp.1 + temp.2)/3,
         Elect.Conduct = (EC + EC.1 + EC.2)/3) %>%
  select(Block = Site, Treatment, Moisture, Temp, Elect.Conduct)

# moisture analysis. Trend for more moisture in unexposed plots, but this is not significant
with(wind.sep.18.2012, interaction.plot(Treatment, Block, Moisture, fun=mean, col=1:10))
sep.18.moist.lmer <- lmer(Moisture ~ Treatment + (1|Block), data = wind.sep.18.2012)
summary(sep.18.moist.lmer) # no variation by block
Anova(sep.18.moist.lmer, test.statistic = "F")

sep.18.moist.lm <- lm(Moisture ~ Treatment, data = wind.sep.18.2012)
summary(sep.18.moist.lm) # no quite sig, but a trend for greater moisture in unexposed plots.
plot(sep.18.moist.lm)

# temperature analysis. 
with(wind.sep.18.2012, interaction.plot(Treatment, Block, Temp, fun=mean, col=1:10))
plot(Temp ~ Block, wind.sep.18.2012) # Strong block effect which probably has to do with the day getting warmer.
sep.18.temp.lmer <- lmer(Temp ~ Treatment + (1|Block), data = wind.sep.18.2012)
summary(sep.18.temp.lmer) # huge block variation
Anova(sep.18.temp.lmer, test.statistic = "F")

# Electrical conductivity.
with(wind.sep.18.2012, interaction.plot(Treatment, Block, Elect.Conduct, fun=mean, col=1:10))
sep.18.EC.lmer <- lmer(Elect.Conduct ~ Treatment + (1|Block), data = wind.sep.18.2012)
summary(sep.18.EC.lmer) # huge block variation
Anova(sep.18.EC.lmer, test.statistic = "F")

## soil moisture, temperature, and electrical conductivity for Sep 22
wind.sep.22.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Lanphere Dunes - Soil Moisture-Temp-EC - Sep 22 2012.csv") %>%
  tbl_df() %>%
  filter(Treatment %in% c("unexposed","exposed")) %>%
  select(Wind.Site:EC.2) %>%
  group_by(Wind.Site, Treatment) %>%
  summarise_each(funs(mean)) %>%
  mutate(Moisture = (moisture + moisture.1 + moisture.2)/3,
         Temp = (temp + temp.1 + temp.2)/3,
         Elect.Conduct = (EC + EC.1 + EC.2)/3) %>%
  select(Block = Wind.Site, Treatment, Moisture, Temp, Elect.Conduct)

# moisture analysis. Trend for more moisture in unexposed plots, but this is not significant
with(wind.sep.22.2012, interaction.plot(Treatment, Block, Moisture, fun=mean, col=1:10))
sep.22.moist.lmer <- lmer(Moisture ~ Treatment + (1|Block), data = wind.sep.22.2012)
summary(sep.22.moist.lmer) # no variation by block
Anova(sep.22.moist.lmer, test.statistic = "F")

sep.22.moist.lm <- lm(Moisture ~ Treatment, data = wind.sep.22.2012)
summary(sep.22.moist.lm) # no quite sig, but a trend for greater moisture in unexposed plots.
plot(sep.22.moist.lm)

# temperature analysis. 
with(wind.sep.22.2012, interaction.plot(Treatment, Block, Temp, fun=mean, col=1:10))
plot(Temp ~ Block, wind.sep.22.2012) # Strong block effect which probably has to do with the day getting warmer.
sep.22.temp.lmer <- lmer(Temp ~ Treatment + (1|Block), data = wind.sep.22.2012)
summary(sep.22.temp.lmer) # huge block variation
Anova(sep.22.temp.lmer, test.statistic = "F")

# Electrical conductivity.
with(wind.sep.22.2012, interaction.plot(Treatment, Block, Elect.Conduct, fun=mean, col=1:10))
sep.22.EC.lmer <- lmer(Elect.Conduct ~ Treatment + (1|Block), data = wind.sep.22.2012)
summary(sep.22.EC.lmer) # huge block variation
Anova(sep.22.EC.lmer, test.statistic = "F")


