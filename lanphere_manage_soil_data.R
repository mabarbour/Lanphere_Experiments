## manage soil characteristics for lanphere experiments (mostly wind)

## load required libraries
library(dplyr)
library(tidyr)

## manage soil nutrient data
prs <- read.csv('~/Documents/Lanphere_Experiments/wind_experiment/data/PRS Probes/PRS_probe_data_2012.csv') %>%
  tbl_df() %>%
  select(Sample.Label:Cd) %>%
  mutate(wind.treat = ifelse(Treatment == "e", "exposed", "unexposed"),
         plot.id = paste(Wind.Site, wind.treat, sep = "_")) %>%
  select(Block = Wind.Site, Treatment = wind.treat, plot.id, Total.N:Cd)
glimpse(prs)

## manage total organic matter
windTOM <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Wind - Soil Total Organic Matter.csv") %>%
  tbl_df() %>%
  filter(BAD.Sample == "n") %>%
  mutate(OvenWt = Cruc....Oven.Dry.Wt. - Crucible.Wt.,
         IgniteWt = Cruc....Ignited.Wt. - Crucible.Wt.,
         PercTOM = (OvenWt - IgniteWt)/OvenWt,
         plot.id = paste(Wind.Site.No., Wind.Treatment, sep = "_")) %>%
  select(plot.id, percent.TOM = PercTOM)
glimpse(windTOM)

## soil moisture, temperature, and electrical conductivity for Sep 18
wind.sep.18.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Lanphere Dunes - Soil Moisture-Temp-EC - Sep 18 2012.csv") %>%
  tbl_df() %>%
  filter(Treatment %in% c("unexposed","exposed")) %>%
  select(Measurement.Time:EC.2) %>%
  group_by(Site, Treatment) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE))) %>%
  mutate(moisture.vwc.Sep18_2012 = (moisture + moisture.1 + moisture.2)/3,
         temp.C.Sep18_2012 = (temp + temp.1 + temp.2)/3,
         EC.Sep18_2012 = (EC + EC.1 + EC.2)/3,
         plot.id = paste(Site, Treatment, sep = "_")) %>%
  ungroup() %>%
  select(plot.id, moisture.vwc.Sep18_2012, temp.C.Sep18_2012, EC.Sep18_2012)
glimpse(wind.sep.18.2012)

## soil moisture, temperature, and electrical conductivity for Sep 22
wind.sep.22.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Lanphere Soil Measurements/Lanphere Dunes - Soil Moisture-Temp-EC - Sep 22 2012.csv") %>%
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

## soil moisture, temperature, electrical conductivity for 2013 in both wind and ant-aphid experiments
soil.2013 <- read.csv('~/Documents/Lanphere_Experiments/wind_experiment/data/wind_soil_2013.csv', skip = 1) %>%
  tbl_df() %>%
  select(-survey.time) %>%
  group_by(experiment, date, block, treatment) %>%
  summarise_each(funs(mean_na.rm = mean(., na.rm = TRUE))) 

aa.soil.2013 <- soil.2013 %>%
  filter(experiment == "ant_aphid") %>%
  mutate(plot.id = paste(block, treatment, sep = "_")) %>%
  ungroup() %>%
  select(block, treatment, plot.id, moisture.vwc.July9_13 = moisture.vwc,
         temp.C.July9_13 = temp.C, EC.July9_13 = electrical.conductivity) 

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

## join all of the wind soil data together
wind.soil.df <- left_join(prs, windTOM) %>%
  left_join(., wind.sep.18.2012) %>%
  left_join(., wind.sep.22.2012) %>%
  left_join(., wind.soil.earlyJuly_13) %>%
  left_join(., wind.soil.July7_13) %>%
  left_join(., wind.soil.July9_13)

## write .csv files
write.csv(aa.soil.2013, "~/Documents/Lanphere_Experiments/final_data/ant_aphid_soil_2013.csv")
write.csv(wind.soil.df, "~/Documents/Lanphere_Experiments/final_data/wind_soil_df.csv")
