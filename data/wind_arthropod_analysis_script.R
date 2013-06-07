#### This script analyzes arthropod data collected from my wind experiment

# set working directory and load required libraries
setwd("~/Documents/Lanphere_Experiments/data")
library("vegan")
library("BiodiversityR")
library("mvabund")
library("MASS")
library("bipartite")
library("lme4")

# source in managed data sets
source('wind_manage_arthropod_dissection_data.R')
source('wind_manage_visual_arthropod_surveys.R')
source('wind_manage_plant_trait_data.R')

###### Merge and manage data sets
# merge all of the arthropod and plant info data together
all_arthropods <- merge(arthropod_wind_visual_data_2012, compact_wind_miners, by.x = "plant_code", by.y = "ID_plant_miners", all.x=TRUE)
all_arthropods <- merge(all_arthropods, compact_wind_ties, by.x = "plant_code", by.y = "ID_plant_ties", all.x = TRUE)

# change all of the NA's in the columns to zeros.  These actually reflect real zeros but the reason for the NA's is because of an artifact from merging the data sets together
all_arthropods[ ,7:31][is.na(all_arthropods[ ,7:31])] <- 0

# focus on herbivore columns. Includes "ALL" instead of just living Caloptilia and tortricids
all_herbivore_data_frame <- with(all_arthropods, cbind.data.frame(leaf_hoppers_blue_brown, Leaf.Hoppers.Camo, Light.Green.Aphids, stem_gall, volcano_gall, mite_gall, yellow_knife_leaf_hopper, twisty_gall, giant_willow_aphid, psyllidae, sawflies, LH.nymph.other, aphid_species_used_in_ant_aphid_experiment, all_Caloptilia, all_tortricid)) 


all_arthropods$herbivore_rich <- specnumber(all_herbivore_data_frame)
all_arthropods$herbivore_abund <- rowSums(all_herbivore_data_frame)

# focus on predator/parasitoid columns
all_predator_data_frame <- with(all_arthropods, cbind.data.frame(Black.Ants,Red.ants,Spiders.black.with.yellow.legs, Spiders.Other, sp_327_LT.1, LTF_329, LEF_eulophid_cocoons, sp_323_all))

all_arthropods$enemy_rich <- specnumber(all_predator_data_frame)
all_arthropods$enemy_abund <- rowSums(all_predator_data_frame)

# arthropod richness and abundance in general
all_arthropods$all_rich <- all_arthropods$herbivore_rich + all_arthropods$enemy_rich
all_arthropods$all_abund <- all_arthropods$herbivore_abund + all_arthropods$enemy_abund

# predator/prey richness. added one to all values to deal with instances of zero species for herbivores.
all_arthropods$predator_prey_richness_ratio <- (all_arthropods$enemy_rich+1)/(all_arthropods$herbivore_rich + 1)

# all arthropods but no dead plants included.  Weird, but it is screwing up my ANOVA analysis, so right now I'm not comfortable using it.
no_dead_all_arthropods <- subset(all_arthropods, Dead. == 0)

########## Exploratory Plots
# richness
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, herbivore_rich, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, enemy_rich, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, all_rich, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, predator_prey_richness_ratio, col=1:10, lwd=2))


# abundance
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, herbivore_abund, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, enemy_abund, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, all_abund, col=1:10, lwd=2))


# species specific abundances
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, all_Caloptilia, col=1:10, lwd=2))
with(all_arthropods, interaction.plot(Wind.Exposure, Genotype, Light.Green.Aphids, col=1:10, lwd=2))

# multivariate abundances, note that the order is reversed along the y-axis (don't know why). Also, doesn't currently account for dead plants.
#herb_mvabund <- mvabund(all_herbivore_data_frame)
#boxplot(herb_mvabund)

#pred_mvabund <- mvabund(all_predator_data_frame)
#boxplot(pred_mvabund)

########## Data analysis. log and sqrt transformations aren't working with the data...DON'T FORGET THAT QUANTIFYING LINKS AND NETWORK SIZE, ETC. MAY BE A WAY TO EXTRACT MORE INFORMATION FROM SIMPLY RICHNESS OF DIFFERENT TROPHIC LEVELS.
# richness data
aov_herb_rich <- aov(log(herbivore_rich+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_herb_rich)

aov_enemy_rich <- aov(log(enemy_rich+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_enemy_rich)

aov_all_rich <- aov(log(all_rich+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_all_rich)

# genearalized mixed effect model for richness
all_arthropods$block <- as.factor(all_arthropods$Block)
glmer_all_rich_null <- glmer(all_rich ~ 1 + (1 | Wind.Exposure:block), data=all_arthropods, family=poisson)

glmer_all_rich_genotype <- glmer(all_rich ~ Genotype + (1 | Wind.Exposure:block), data=all_arthropods, family=poisson)

glmer_all_rich_exposure <- glmer(all_rich ~ Wind.Exposure + (1 | Wind.Exposure:block), data=all_arthropods, family=poisson)

anova(glmer_all_rich_null, glmer_all_rich_genotype)

#

aov_predator_prey_rich_ratio <- aov(log(predator_prey_richness_ratio + 1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_predator_prey_rich_ratio)

# abundance data
aov_herb_abund <- aov(log(herbivore_abund+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_herb_abund)

aov_enemy_abund <- aov(log(enemy_abund +1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_enemy_abund)

aov_all_abund <- aov(log(all_abund + 1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_all_abund)

aov_Caloptilia_abund <- aov(log(all_Caloptilia+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_Caloptilia_abund)
aov_Light.Green.Aphid_abund <- aov(Light.Green.Aphids ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= all_arthropods, qr=TRUE)
summary(aov_Light.Green.Aphid_abund)

# dissimilarity. NOTHING is showing through.  Possibly because I had to do the weird transformation of all zero values to create highly similar or dissimilar communities.
# transform all zero values
#bray_all_arthropods <- vegdist(no_dead_all_arthropods[ ,7:31], method="bray")
#no_zero_arthropod_dissimilarity <- dist.zeroes(no_dead_all_arthropods[ ,7:31], bray_all_arthropods)
#anosim_genotype_arthropods <- anosim(no_zero_arthropod_dissimilarity, all_arthropods$Genotype)
#anosim_exposure_arthropods <- anosim(no_zero_arthropod_dissimilarity, all_arthropods$Wind.Exposure)


###################### Caloptilia data set only
Caloptilia_only <- merge(arthropod_wind_visual_data_2012, sub_Caloptilia, by.x ="plant_code", by.y = "ID_plant_miners", all.x=TRUE)
Caloptilia_only[ ,24:28][is.na(Caloptilia_only[ ,24:28])] <- 0
omit_Caloptilia <- Caloptilia_only[-c(32,55), ]
Caloptilia_only <- merge(Caloptilia_only, sub_plant_traits, by = "plant_code")

###### plots of abundance of different stages. Abundance appears to switch between early and late instars
# summary data for fresh_gallery_mines
gallery_treatment_means <- aggregate(fresh_gallery_mine ~ Genotype * Wind.Exposure, data=Caloptilia_only, mean)
gallery_exposure_means <- aggregate(fresh_gallery_mine ~ Wind.Exposure, data=gallery_treatment_means, mean)
gallery_exposure_sd <- aggregate(fresh_gallery_mine ~ Wind.Exposure, data=gallery_treatment_means, sd)
gallery_exposure_means$standard_error <- gallery_exposure_sd$fresh_gallery_mine/sqrt(10)

# interaction plot of fresh gallery mines
with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, fresh_gallery_mine, col=1:10, lwd=4, legend=F, ylab="Leaf mine Abundance", xlab="", las=1))

# add summary points and error bars
points(x=0.95, gallery_exposure_means[1,2], pch=19)
arrows(x0=0.95, y0=gallery_exposure_means[1,2], x1=0.95, y1=gallery_exposure_means[1,2]+gallery_exposure_means[1,3], angle=90, length=0.1)
arrows(x0=0.95, y0=gallery_exposure_means[1,2], x1=0.95, y1=gallery_exposure_means[1,2]-gallery_exposure_means[1,3], angle=90, length=0.1)
points(x=2.05, gallery_exposure_means[2,2], pch=19)
arrows(x0=2.05, y0=gallery_exposure_means[2,2], x1=2.05, y1=gallery_exposure_means[2,2]+gallery_exposure_means[2,3], angle=90, length=0.1)
arrows(x0=2.05, y0=gallery_exposure_means[2,2], x1=2.05, y1=gallery_exposure_means[2,2]-gallery_exposure_means[2,3], angle=90, length=0.1)

# summary data for leaf tip folders
LTF_treatment_means <- aggregate(LTF_Caloptilia ~ Genotype * Wind.Exposure, data=Caloptilia_only, mean)
LTF_exposure_means <- aggregate(LTF_Caloptilia ~ Wind.Exposure, data=LTF_treatment_means, mean)
LTF_exposure_sd <- aggregate(LTF_Caloptilia ~ Wind.Exposure, data=LTF_treatment_means, sd)
LTF_exposure_means$standard_error <- LTF_exposure_sd$LTF_Caloptilia/sqrt(10)

# interaction plot of leaf tip folders
with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, LTF_Caloptilia, col=1:10, lwd=3, legend=F, ylab="Leaf tip fold Abundance", xlab="", las=1))

# add summary points and error bars
points(x=0.95, LTF_exposure_means[1,2], pch=19)
arrows(x0=0.95, y0=LTF_exposure_means[1,2], x1=0.95, y1=LTF_exposure_means[1,2]+LTF_exposure_means[1,3], angle=90, length=0.1)
arrows(x0=0.95, y0=LTF_exposure_means[1,2], x1=0.95, y1=LTF_exposure_means[1,2]-LTF_exposure_means[1,3], angle=90, length=0.1)
points(x=2.05, LTF_exposure_means[2,2], pch=19)
arrows(x0=2.05, y0=LTF_exposure_means[2,2], x1=2.05, y1=LTF_exposure_means[2,2]+LTF_exposure_means[2,3], angle=90, length=0.1)
arrows(x0=2.05, y0=LTF_exposure_means[2,2], x1=2.05, y1=LTF_exposure_means[2,2]-LTF_exposure_means[2,3], angle=90, length=0.1)

with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, LEF_Caloptilia, col=1:10, lwd=2))
with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, tent_mine_Caloptilia, col=1:10, lwd=2))

# Data analysis
aov_fresh_gallery <- aov(log(fresh_gallery_mine + 1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= Caloptilia_only)
summary(aov_fresh_gallery)
plot(fitted(aov_fresh_gallery[[4]]), studres(aov_fresh_gallery[[4]]))
abline(h=0, lty=2)

aov_LTF <- aov(log(LTF_Caloptilia +1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data=Caloptilia_only, qr=TRUE)
summary(aov_LTF)
plot(fitted(aov_LTF[[4]]), studres(aov_LTF[[4]]))
abline(h=0, lty=2)

###### plot plant traits
# summary data for mature shoots
mature_shoot_treatment_means <- aggregate(mature_shoot_total ~ Genotype * Wind.Exposure, data=Caloptilia_only, mean)
mature_shoot_exposure_means <- aggregate(mature_shoot_total ~ Wind.Exposure, data=mature_shoot_treatment_means, mean)
mature_shoot_exposure_sd <- aggregate(mature_shoot_total ~ Wind.Exposure, data=mature_shoot_treatment_means, sd)
mature_shoot_exposure_means$standard_error <- mature_shoot_exposure_sd$mature_shoot_total/sqrt(10)

# mature shoot total
with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, mature_shoot_total, col=1:10, lwd=3, legend=F, ylab="Mature shoot growth (cm)", xlab="", las=1))

# add summary points and error bars
points(x=0.95, mature_shoot_exposure_means[1,2], pch=19)
arrows(x0=0.95, y0=mature_shoot_exposure_means[1,2], x1=0.95, y1=mature_shoot_exposure_means[1,2]+mature_shoot_exposure_means[1,3], angle=90, length=0.1)
arrows(x0=0.95, y0=mature_shoot_exposure_means[1,2], x1=0.95, y1=mature_shoot_exposure_means[1,2]-mature_shoot_exposure_means[1,3], angle=90, length=0.1)
points(x=2.05, mature_shoot_exposure_means[2,2], pch=19)
arrows(x0=2.05, y0=mature_shoot_exposure_means[2,2], x1=2.05, y1=mature_shoot_exposure_means[2,2]+mature_shoot_exposure_means[2,3], angle=90, length=0.1)
arrows(x0=2.05, y0=mature_shoot_exposure_means[2,2], x1=2.05, y1=mature_shoot_exposure_means[2,2]-mature_shoot_exposure_means[2,3], angle=90, length=0.1)

# summary data for basal shoots
immature_shoot_treatment_means <- aggregate(immature_shoot_total ~ Genotype * Wind.Exposure, data=Caloptilia_only, mean)
immature_shoot_exposure_means <- aggregate(immature_shoot_total ~ Wind.Exposure, data=immature_shoot_treatment_means, mean)
immature_shoot_exposure_sd <- aggregate(immature_shoot_total ~ Wind.Exposure, data=immature_shoot_treatment_means, sd)
immature_shoot_exposure_means$standard_error <- immature_shoot_exposure_sd$immature_shoot_total/sqrt(10)

# shoots sprouting at base
with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, immature_shoot_total, col=1:10, lwd=2, legend=F, ylab="Basal shoot growth (cm)", xlab="", las=1))

# add summary points and error bars
points(x=0.95, immature_shoot_exposure_means[1,2], pch=19)
arrows(x0=0.95, y0=immature_shoot_exposure_means[1,2], x1=0.95, y1=immature_shoot_exposure_means[1,2]+immature_shoot_exposure_means[1,3], angle=90, length=0.1)
arrows(x0=0.95, y0=immature_shoot_exposure_means[1,2], x1=0.95, y1=immature_shoot_exposure_means[1,2]-immature_shoot_exposure_means[1,3], angle=90, length=0.1)
points(x=2.05, immature_shoot_exposure_means[2,2], pch=19)
arrows(x0=2.05, y0=immature_shoot_exposure_means[2,2], x1=2.05, y1=immature_shoot_exposure_means[2,2]+immature_shoot_exposure_means[2,3], angle=90, length=0.1)
arrows(x0=2.05, y0=immature_shoot_exposure_means[2,2], x1=2.05, y1=immature_shoot_exposure_means[2,2]-immature_shoot_exposure_means[2,3], angle=90, length=0.1)

###### analysis on plant traits
aov_mature_shoots <- aov(log(mature_shoot_total+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data=Caloptilia_only, qr=TRUE)
summary(aov_mature_shoots) # no genotype by environment interaction with mature shoot totals
plot(fitted(aov_mature_shoots[[4]]), studres(aov_mature_shoots[[4]]))
abline(h=0, lty=2)

aov_sprouts <- aov(log(immature_shoot_total+1) ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data=Caloptilia_only, qr=TRUE)
summary(aov_sprouts)
plot(fitted(aov_sprouts[[4]]), studres(aov_sprouts[[4]]))
abline(h=0, lty=2)

###### look for correlations between plant traits and Caloptilia abundances
traits_Caloptilia_data <- aggregate(cbind(immature_shoot_total,mature_shoot_total, Mature.Green.Leaves, LTF_Caloptilia, fresh_gallery_mine) ~ Genotype * Wind.Exposure, data=Caloptilia_only, mean)

# plot and look at correlations between leaf tip fold abundance versus mature green leaves
plot(LTF_Caloptilia ~ Mature.Green.Leaves, traits_Caloptilia_data)
cor.test(log(traits_Caloptilia_data$LTF_Caloptilia+1), log(traits_Caloptilia_data$Mature.Green.Leaves+1))
lm_LTF_leaves <- lm(log(traits_Caloptilia_data$LTF_Caloptilia+1) ~ log(traits_Caloptilia_data$Mature.Green.Leaves+1))
summary(lm_LTF_leaves) # mature leaf abundance explains about 20% of the variation.

# look at correlation with mature shoot total
plot(LTF_Caloptilia ~ mature_shoot_total, traits_Caloptilia_data)
cor.test(log(traits_Caloptilia_data$LTF_Caloptilia+1), log(traits_Caloptilia_data$mature_shoot_total+1)) 
lm_LTF_shoots <- lm(log(traits_Caloptilia_data$LTF_Caloptilia+1) ~ log(traits_Caloptilia_data$mature_shoot_total+1))
summary(lm_LTF_shoots) # mature shoot total explains about 23% of the variation

# look at correlation between mature shoot and immature shoots
plot(immature_shoot_total ~ mature_shoot_total, traits_Caloptilia_data)
cor.test(log(traits_Caloptilia_data$immature_shoot_total+1), log(traits_Caloptilia_data$mature_shoot_total+1)) # negative correlation between mature and immature shoot totals
plot(fresh_gallery_mine ~ LTF_Caloptilia, traits_Caloptilia_data)
cor.test(log(traits_Caloptilia_data$fresh_gallery_mine+1), log(traits_Caloptilia_data$LTF_Caloptilia+1)) # no negative correlation between abundances of early versus late larva instars.

# plot and look at correlations between gallery mine and sprouting shoots. Not looking great...
plot(fresh_gallery_mine ~ immature_shoot_total, traits_Caloptilia_data)
cor.test(log(traits_Caloptilia_data$fresh_gallery_mine+1), log(traits_Caloptilia_data$immature_shoot_total+1))
lm_gallery_sprouts <- lm(log(traits_Caloptilia_data$fresh_gallery_mine+1) ~ log(traits_Caloptilia_data$immature_shoot_total+1))
summary(lm_gallery_sprouts) # marginally significant effect.  only explains about 11% of the variance (adjusted R-squared)

# create new variables: sum early instars and later instars
#Caloptilia_only$early_instars <- Caloptilia_only$gallery_mine_Caloptilia + Caloptilia_only$tent_mine_Caloptilia
#Caloptilia_only$late_instars <- Caloptilia_only$LTF_Caloptilia + Caloptilia_only$LEF_Caloptilia

# plots of data
#with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, early_instars, col=1:10, lwd=2))
#with(Caloptilia_only, interaction.plot(Wind.Exposure, Genotype, late_instars, col=1:10, lwd=2))


