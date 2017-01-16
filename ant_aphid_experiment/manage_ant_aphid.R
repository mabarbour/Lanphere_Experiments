## load libraries
library(dplyr)
library(tidyr)
library(lme4) # mixed-effect modelling
library(car) # for Anova()
library(RLRsim) # for random effects testing

## upload data

# experimental design info
setup <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Ant-Aphid_Experiment_Setup.csv") %>%
  tbl_df() %>% 
  select(Block, Aphid.Treat = Aphids.or.No.Aphids, Ant.Mound.Dist = Distant.to.Ant.Mound,
         Genotype, Plant_Position = Plant.Position) %>%
  mutate(Block.Plant_Position = as.numeric(paste(Block, Plant_Position, sep=".")),
         Block = as.factor(Block))

# visual arthropod surveys. Appears to be missing last survey of the year
visual <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Ant-Aphid_Data_v2.csv", stringsAsFactors=F) %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block),
         total.aphid.abund = Green.Aphids + Orange.Aphids,
         total.ant.abund = Red.Ants + Black.Ants,
         total.frog.hoppers = Leaf.Hoppers.Camo + Spittlebugs,
         total.tortricid = Leaf.Bundlers + Leaf.Edge.Rollers) %>%
  select(relative.Date:Aphid.Treatment, Unique.Number, Dead, Aphids_Added:Spiders, LTF.abund, Leaf.Hoppers, 
         total.frog.hoppers, total.tortricid, Giant.Willow.Aphids, total.aphid.abund, Fly.Predators,
         Chalcid.Parasitoids, Beetles, Other.Predators:Other.Herbivores, 
         total.ant.abund) %>% # excluding Tent.Mines because of detectability issues during surveys, they also should be indicative of LTF of abundance. I also excluded Old.Herbivores and Other.Arthropods, because I don't know what these categories mean right now. Leaf.Hoppers accounts for blue, brown and black/yellow leaf hopper nymphs. Spiders includes Spider.blk.ylw and Spiders.Other. Only a couple instances whether this underestimates spider richness (should be 2 instead of 1).
  group_by(Block, Distance.to.Ant.Mound, Genotype, Aphid.Treatment, Unique.Number) %>%
  summarise_each(funs(max.narm = max(., na.rm = TRUE))) %>% # take the maximum number of individuals observed on a plant over the entire survey period as a conservative estimate of abundance over the entire growing season. is this doing the appropriate thing for total predator and non_aphid.herb.abund?
  mutate(non_aphid.herb.abund = LTF.abund + Leaf.Hoppers + total.frog.hoppers + total.tortricid + Giant.Willow.Aphids + Other.Herbivores,
         non_ant.predator.abund = Chalcid.Parasitoids + Other.Predators + Spiders) %>%
  filter(Dead < 1) # exclude plants that were classified as dead on final survey of the year.

# strong genotype and aphid effect on ant abundance. No distance to mound effect. 
# for red ants only, there is a significant Aphid, and marginal Distance to ant mound effect. Possible GxE for aphid treatment...
# marginal genotype effect on spider abundance.
# strong genotype effect and distance to ant mound on LTF abundance
# strong genotype effect and marginal Ant.Mound x Aphid interaction on non-aphid herb abund
# genotype effect on non-ant predator abundance
# distance to ant.mound, Genotype, and marginal Distance x Aphid effect on leaf hopper abundance
# close to marginal genotype effect on aphid abundance...consider modelling weekly population growth instead of maximum abundance over the season.
# Genotype and distance to ant mound effect on total tortricid abundance
with(visual, interaction.plot(Distance.to.Ant.Mound, Genotype, total.tortricid, mean))
plot(LTF.abund ~ Distance.to.Ant.Mound, visual)
explore.lmer <- lmer(log(total.tortricid+1) ~ Distance.to.Ant.Mound*Aphid.Treatment*Genotype + (1|Block),
                 visual)
summary(explore.lmer)
Anova(explore.lmer, test.statistic = "F")
plot(explore.lmer)


# dissection data
#LTFdata <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Arthropod Dissections/Leaf Miner Dissections - Ant-Aphid Experiment 2012.csv", stringsAsFactors=F)

# plant architecture data
arch <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/ant_aphid_architecture_data.csv") %>% 
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  select(Block, Aphid.Treatment = Aphids.or.No.Aphids, Ant.Dist = Distant.to.Ant.Mound, Genotype, Survived, Mature.Green.Leaves, Shoots.Sprouting, Max.Height, mature.shoot.1:mature.shoot.18)
arch$total.mature.shoot.length = rowSums(select(arch, mature.shoot.1:mature.shoot.18))
arch$mature.shoot.count = rowSums(select(arch, mature.shoot.1:mature.shoot.18) > 0)
arch$max.mature.shoot.length = apply(X = select(arch, mature.shoot.1:mature.shoot.18),
                                     MARGIN = 1, FUN = max)

# exploratory archicture analysis
# strong genotype effect on plant height and marginal effect of ant.distance
# strong genotype effect on number of shoots and marginal effect of aphid treatment
# strong genotype effect and Aphid x Ant.Dist interaction on sprouting shoots and maximum mature shoot length
# strong genotype effect on total mature shoot length
plot(Max.Height ~ Genotype, arch)
plot(Max.Height ~ Ant.Dist, arch)
with(arch, interaction.plot(Ant.Dist, Aphid.Treatment, Shoots.Sprouting, 
                            fun = function(x) mean(x, na.rm = TRUE), col = 1:10))
exp.arch.lmer <- lmer(log(Max.Height) ~ Aphid.Treatment + log(Ant.Dist) + Genotype + (1|Block),
                      arch)
summary(exp.arch.lmer)
Anova(exp.arch.lmer, test.statistic = "F")

with(arch, table(Genotype, Aphid.Treatment, Survived))
exp.arch.glmer <- glmer(Survived ~ Aphid.Treatment*Ant.Dist*Genotype + (1|Block), arch, family = binomial) # model failing to converge

exp.arch.lm <- lm(log(Max.Height) ~ log(Ant.Dist) + Genotype, arch)
summary(exp.arch.lm)
library(visreg)
visreg(exp.arch.lm, trans = exp)


# plant trait data
trich <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Trichome Density - Ant-Aphid Experiment 2012.csv") %>%
  tbl_df() %>%
  select(Block.Plant_Position = Collection.No., Trichome.Density)
DryLeafWt <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Leaf Weights - Dry - Ant-Aphid Experiment 2012.csv") %>% 
  tbl_df() %>%
  select(Block.Plant_Position = Collection.No., Dry.Leaf.Weight.g = Dry.Leaf.Weight)
WetLeafWt <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Ant-aphid Plant Trait Data.csv") %>% 
  tbl_df() %>%
  select(Block.Plant_Position = Collection.Number, Wet.Leaf.Weight.g = Wet.Leaf.Weight) # didnt' retain percent browning, because I don't know how informative it will be since we didn't take leaves from plants that had leaves that were substantially browned.
WetLeafWt[3,"Wet.Leaf.Weight.g"] <- 0.103 # typo, double-checked this with real data.

leafwts <- full_join(DryLeafWt, WetLeafWt, by = "Block.Plant_Position") %>%
  mutate(LDMC = Dry.Leaf.Weight.g/Wet.Leaf.Weight.g,
         WC = (Wet.Leaf.Weight.g - Dry.Leaf.Weight.g)/Dry.Leaf.Weight.g)

leaftraits <- left_join(leafwts, trich, by = "Block.Plant_Position")

leaftraits.df <- left_join(setup, leaftraits, by = "Block.Plant_Position")


## Exploratory data analysis
plot(WC ~ Genotype, leaftraits.df)
hist(leaftraits.df$Wet.Leaf.Weight.g)
which(leaftraits.df$Block.Plant_Position == 5.47)
leaftraits.df[293, ] # may be a typo

# no effect on water content
WC.lmer <- lmer(WC ~ Ant.Mound.Dist*Aphid.Treat*Genotype + (1|Block), data = leaftraits.df[-293, ])
summary(WC.lmer)
Anova(WC.lmer, test = "F")

# no effect on LDMC
LDMC.lmer <- lmer(LDMC ~ Ant.Mound.Dist*Aphid.Treat*Genotype + (1|Block), data = leaftraits.df[-293, ])
summary(LDMC.lmer)
Anova(LDMC.lmer, test = "F")

cor.test(leaftraits.df$LDMC[-293], leaftraits.df$WC[-293]) # strong negative correlation, suggesting they may be redundant. Maybe just use water content for manuscript

# strong Genotype effect on trichome density
trich.lmer <- lmer(Trichome.Density ~ Ant.Mound.Dist*Aphid.Treat*Genotype + (1|Block), data = leaftraits.df)
summary(trich.lmer)
Anova(trich.lmer, test = "F")
