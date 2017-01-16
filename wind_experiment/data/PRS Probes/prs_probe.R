## load required libraries
library(dplyr)
library(lme4)
library(car)

## upload data
prs <- read.csv('~/Documents/Lanphere_Experiments/wind_experiment/data/PRS Probes/PRS_probe_data_2012.csv') %>%
  tbl_df() %>%
  select(Sample.Label:Cd) %>%
  mutate(Block = as.factor(Wind.Site), Sample.Label = as.factor(Sample.Label)) %>%
  select(Block, Treatment, Sample.Label, Total.N:Cd)

## exploratory data analysis
with(prs, interaction.plot(Treatment, Block, NH4.N, col = 1:10)) 

# Total.N appears to have the strongest variation among treatment. None of the other micronutrients show much variation.
Total.N.lmer <- lmer(Total.N ~ Treatment + (1|Block), prs) 
summary(Total.N.lmer)
Anova(Total.N.lmer, test.statistic = "F") # marginally significant
plot(Total.N.lmer)

