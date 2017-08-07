############################################
## Description: This script manages data for soil characteristics of Lanphere experiments (mostly wind).
## Code author: Matt Barbour
## Email: barbour@zoology.ubc.ca
############################################

#### Load required libraries ----
library(dplyr)
library(tidyr)
library(psych) # for correlations
library(car) # for visualizing correlations
library(missMDA) # for imputing missing values in PCA

#### Manage soil nutrient data ----
prs <- read.csv('manage_raw_data/raw_data/PRS_probe_data_2012.csv') %>%
  tbl_df() %>%
  select(Sample.Label:Cd) %>%
  mutate(wind.treat = ifelse(Treatment == "e", "exposed", "unexposed"),
         plot.id = paste(Wind.Site, wind.treat, sep = "_")) %>%
  select(Block = Wind.Site, Treatment = wind.treat, plot.id, Total.N:Cd)
glimpse(prs)

#### Manage total organic matter ----
windTOM <- read.csv('manage_raw_data/raw_data/Wind_Soil_Total_Organic_Matter.csv') %>%
  tbl_df() %>%
  filter(BAD.Sample == "n") %>%
  mutate(OvenWt = Cruc....Oven.Dry.Wt. - Crucible.Wt.,
         IgniteWt = Cruc....Ignited.Wt. - Crucible.Wt.,
         PercTOM = (OvenWt - IgniteWt)/OvenWt,
         plot.id = paste(Wind.Site.No., Wind.Treatment, sep = "_")) %>%
  select(plot.id, percent.TOM = PercTOM)
glimpse(windTOM)

#### Soil moisture, temperature, and electrical conductivity ----

## September 18th, 2012
wind.sep.18.2012 <- read.csv('manage_raw_data/raw_data/Lanphere_Dunes_Soil_Moisture-Temp-EC_Sep_18_2012.csv') %>%
  tbl_df() %>%
  filter(Treatment %in% c("unexposed","exposed")) %>%
  select(Measurement.Time:EC.2) %>%
  group_by(Site, Treatment) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  mutate(moisture.vwc.Sep18_2012 = (moisture + moisture.1 + moisture.2)/3,
         temp.C.Sep18_2012 = (temp + temp.1 + temp.2)/3,
         EC.Sep18_2012 = (EC + EC.1 + EC.2)/3,
         plot.id = paste(Site, Treatment, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.Sep18_2012, temp.C.Sep18_2012, EC.Sep18_2012)
glimpse(wind.sep.18.2012)

## Sep 22, 2012
wind.sep.22.2012 <- read.csv('manage_raw_data/raw_data/Lanphere_Dunes_Soil_Moisture-Temp-EC_Sep_22_2012.csv') %>%
  tbl_df() %>%
  filter(Treatment %in% c("unexposed","exposed")) %>%
  select(Wind.Site:EC.2) %>%
  group_by(Wind.Site, Treatment) %>%
  summarise_each(funs(mean)) %>%
  mutate(moisture.vwc.Sep22_2012 = (moisture + moisture.1 + moisture.2)/3,
         temp.C.Sep22_2012 = (temp + temp.1 + temp.2)/3,
         EC.Sep22_2012 = (EC + EC.1 + EC.2)/3,
         plot.id = paste(Wind.Site, Treatment, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.Sep22_2012, temp.C.Sep22_2012, EC.Sep22_2012)
glimpse(wind.sep.22.2012)

## 2013 (both wind and ant-aphid experiments)
soil.2013 <- read.csv('manage_raw_data/raw_data/wind_soil_2013.csv', skip = 1) %>%
  tbl_df() %>%
  select(-survey.time) %>%
  group_by(experiment, date, block, treatment) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) 

aa.soil.2013 <- soil.2013 %>%
  filter(experiment == "ant_aphid") %>%
  mutate(plot.id = paste(block, treatment, sep = "_")) %>%
  ungroup() %>%
  select(Block = block, Ant.mound.dist = treatment, Plot_code = plot.id, moisture.vwc.July9_13 = moisture.vwc) # note that temperature and electrical conductivity data are not available. 

wind.soil.July7_13 <- soil.2013 %>%
  filter(experiment == "wind", date == "07-Jul-13") %>%
  mutate(Treat = ifelse(treatment == "E", "exposed", "unexposed"),
         plot.id = paste(block, Treat, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.July7_13 = moisture.vwc, temp.C.July7_13 = temp.C,
         EC.July7_13 = electrical.conductivity)

wind.soil.July9_13 <- soil.2013 %>%
  filter(experiment == "wind", date == "09-Jul-13") %>%
  mutate(Treat = ifelse(treatment == "E", "exposed", "unexposed"),
         plot.id = paste(block, Treat, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.July9_13 = moisture.vwc, temp.C.July9_13 = temp.C,
         EC.July9_13 = electrical.conductivity)

wind.soil.earlyJuly_13 <- soil.2013 %>%
  filter(experiment == "wind", date == "before July 7, 2013") %>%
  mutate(Treat = ifelse(treatment == "E", "exposed", "unexposed"),
         plot.id = paste(block, Treat, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.earlyJuly_13 = moisture.vwc, temp.C.earlyJuly_13 = temp.C,
         EC.earlyJuly_13 = electrical.conductivity)

#### Join all of the wind data together ----
wind.soil.df <- left_join(prs, windTOM) %>%
  left_join(., wind.sep.18.2012) %>%
  left_join(., wind.sep.22.2012) %>%
  left_join(., wind.soil.earlyJuly_13) %>%
  left_join(., wind.soil.July7_13) %>%
  left_join(., wind.soil.July9_13) %>%
  mutate(avg.moisture.vwc.2012 = (moisture.vwc.Sep18_2012 + moisture.vwc.Sep22_2012)/2,
         avg.moisture.vwc.2013 = (moisture.vwc.earlyJuly_13 + moisture.vwc.July7_13 + moisture.vwc.July9_13)/3,
         avg.moisture.vwc = (avg.moisture.vwc.2012 + avg.moisture.vwc.2013)/2,
         treat.tmp = ifelse(Treatment == "exposed", "E", "U"),
         Plot_code = paste(Block, treat.tmp, sep = ".")) %>%
  select(Block, Treatment, Plot_code, Total.N:percent.TOM, avg.moisture.vwc.2012, avg.moisture.vwc.2013, avg.moisture.vwc)

with(wind.soil.df, cor.test(avg.moisture.vwc.2012, avg.moisture.vwc.2013)) # r = 0.93, t = 10.908, df = 18, P < 0.001. Since plot-level soil moisture was highly correlated between years, I decided to just take the average for future analyses.

#### Create Principal Components Scores of Soil Nutrients ----
nuts <- colnames(select(wind.soil.df, NO3.N:Cd)) # get nutrient names

scatterplotMatrix(wind.soil.df[ ,nuts]) # visualize nutrient correlations. A lot to examine...
corr.test(wind.soil.df[ ,nuts]) # the most apparent trend is that all of the micronutrients (e.g. Ca to Cd) are negatively correlated with both Nitrogen compounds.

nut.PCA <- princomp(wind.soil.df[ ,nuts], cor = TRUE) # set cor = TRUE, because although they are all in the same units it is not necessarily true that the magnitudes are equivalent in their importance to the plant.
biplot(nut.PCA)
summary(nut.PCA) # 34% of variance explained by PC1 (16% for PC2)
nut.PCA$loadings # positive values of PC1 indicate higher supply rates of micronutrients, but lower supply of nitrogen compounds
hist(nut.PCA$scores[ ,"Comp.1"]) # look at distribution

wind.soil.df <- wind.soil.df %>% 
  mutate(nut.PC1 = nut.PCA$scores[ ,"Comp.1"])

#### Create Principal Component Scores for all soil properties
props.all <- colnames(select(wind.soil.df, NO3.N:percent.TOM, avg.moisture.vwc))

soil.mat <- as.matrix(wind.soil.df[ ,props.all])
#soil.nb <- estim_ncpPCA(soil.mat, method.cv = "loo") # suggests using 2 PCs for imputing data
w.soil.imp <- imputePCA(soil.mat, ncp = 2, scale = TRUE)

props.PCA <- princomp(w.soil.imp$completeObs, cor = TRUE)
biplot(props.PCA)
summary(props.PCA) # 52% of variance explained by first 2 components
props.PCA$loadings[ ,c(1,2)] # positive values of PC1 indicate moist soil with lots of organic matter and higher supply rates of micronutrients, but lower supply of nitrogen compounds. Variation in PC2 describes variation in ionic elements.


wind.soil.df <- wind.soil.df %>% mutate(soil.PC1 = props.PCA$scores[ ,"Comp.1"], soil.PC2 = props.PCA$scores[ ,"Comp.2"], Wind.Exposure = ifelse(Treatment == "exposed","Exposed","Unexposed")) %>% select(Block, Wind.Exposure, Plot_code, Total.N:percent.TOM, avg.moisture.vwc, nut.PC1:soil.PC2)

#### Save .csv files ----
write.csv(aa.soil.2013, "final_data/ant_aphid_soil_2013.csv")
write.csv(wind.soil.df, "final_data/wind_soil_df.csv")
