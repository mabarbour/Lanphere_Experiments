## load libraries
library(dplyr)
library(lme4)
library(car)

## upload data

# plant info
wind_plant_info <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_info.csv") %>%
  tbl_df() %>%
  mutate(block = as.factor(block)) %>%
  select(block, treatment, genotype, plant.position, plant.id)

# spodoptera experiment #1
spod.1 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/spodoptera/spodoptera_herbivory_wind_experiment_1.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp1 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.weight.exp1 = larva.wet.weight, survival.exp1)

# spodoptera experiment #2
spod.2 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/spodoptera/spodoptera_herbivory_wind_experiment_2.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp2 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.weight.exp2 = larvae.wet.weight, survival.exp2)

# join data sets
spod.df.tmp <- left_join(wind_plant_info, spod.1, by = "plant.id")
spod.df <- left_join(spod.df.tmp, spod.2, by = "plant.id")

# exploratory analyses
with(spod.df, table(treatment, genotype, survival.exp1))
with(spod.df, table(treatment, genotype, survival.exp2))

# Genotype effect in first experiment.
with(spod.df, interaction.plot(treatment, genotype, larva.wet.weight.exp1, 
                               fun = function(x) mean(x, na.rm = TRUE), col = 1:10))
wt.1.lmer <- lmer(larva.wet.weight.exp1 ~ treatment + (treatment|genotype) + (1|block), data = spod.df)
summary(wt.1.lmer)
Anova(wt.1.lmer, test.statistic = "F")

surv.1.glmer <- glmer(survival.exp1 ~ genotype*treatment + (1|block), spod.df, family = binomial) # model not converging
surv.1.glm <- glm(survival.exp1 ~ genotype*treatment, spod.df, family = binomial)
summary(surv.1.glm)
anova(surv.1.glm, test = "Chi")

# GxE interaction for experiment 2 Spodoptera weights. Too few data to analyze glmer
with(spod.df, interaction.plot(treatment, genotype, larva.wet.weight.exp2, 
                               fun = function(x) mean(x, na.rm = TRUE), col = 1:10))
wt.2.lmer <- lmer(larva.wet.weight.exp2 ~ treatment + (treatment|genotype) + (1|block), data = spod.df)
summary(wt.2.lmer)
Anova(wt.2.lmer, test.statistic = "F")

surv.2.glmer <- glmer(survival.exp2 ~ genotype*treatment + (1|block), spod.df, family = binomial) # model not converging
surv.2.glm <- glm(survival.exp2 ~ genotype*treatment, spod.df, family = binomial)
summary(surv.2.glm)
anova(surv.2.glm, test = "Chi")
