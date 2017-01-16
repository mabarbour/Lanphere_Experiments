## upload data

arth.2012 <-read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/arthropod_wind_visual_data_2012.csv") %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  filter(Dead. < 1)

arth.2012$spiders <- with(arth.2012, Spiders.black.with.yellow.legs + Spiders.Other)

with(arth.2012, interaction.plot(Wind.Exposure, Genotype, spiders, mean, col= 1:10 ))

spiders.lmer <- lmer(log(spiders+1) ~ Wind.Exposure + (Wind.Exposure|Genotype) + (1|Block), arth.2012)
summary(spiders.lmer)
Anova(spiders.lmer, test.statistic = "F")
