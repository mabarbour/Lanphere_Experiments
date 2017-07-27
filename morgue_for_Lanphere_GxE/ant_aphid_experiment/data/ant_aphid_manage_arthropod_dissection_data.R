
############ upload and manage data
ant_aphid_leaf_miner_dissection_summer_2012 <- read.csv("~/Documents/Lanphere_Experiments/data/ant_aphid_arthropod_dissection_summer_2012.csv", stringsAsFactor=F, skip=1)

# create variable unique ID for each plant in both data sets
ant_aphid_arthropod_dissection_summer_2012$plant_code <- with(ant_aphid_arthropod_dissection_summer_2012, paste(Block, Plant, sep="_"))

# composite variable for living or freshly dead gallery mines (late season cohort abundance)
ant_aphid_arthropod_dissection_summer_2012$fresh_gallery_mine <- ant_aphid_arthropod_dissection_summer_2012$gallery_mine_living + ant_aphid_arthropod_dissection_summer_2012$gallery_mine_fresh_but_stiff + ant_aphid_arthropod_dissection_summer_2012$gallery_mine_freshly_dead

# composite variable for any evidence of a leaf tip fold (early season cohort abundance)
ant_aphid_arthropod_dissection_summer_2012$LTF_abund <- ant_aphid_arthropod_dissection_summer_2012$LTF_329 + ant_aphid_arthropod_dissection_summer_2012$LTF_empty_with_frass + ant_aphid_arthropod_dissection_summer_2012$LTF_dead + ant_aphid_arthropod_dissection_summer_2012$LTF_other + ant_aphid_arthropod_dissection_summer_2012$LTF_empty

# create a more manageable data set for statistical analyses
compact_ant_aphid_dissect <- subset(ant_aphid_arthropod_dissection_summer_2012, select = c("fresh_gallery_mine", "LTF_abund", "plant_code"))
sum_compact_ant_aphid_dissect <- rowsum(compact_ant_aphid_dissect[ ,1:2], compact_ant_aphid_dissect$plant_code)
sum_compact_ant_aphid_dissect$plant_code <- rownames(sum_compact_ant_aphid_dissect)

# upload plant info data
Ant.Aphid_Experiment_Setup <- read.csv("~/Documents/Lanphere_Experiments/data/Ant-Aphid_Experiment_Setup.csv")
plant_info_ant_aphid_experiment <- subset(Ant.Aphid_Experiment_Setup, select=c("Block","Aphids.or.No.Aphids", "Distant.to.Ant.Mound", "Plant.Position", "Genotype"))
plant_info_ant_aphid_experiment$plant_code <- with(plant_info_ant_aphid_experiment, paste(Block, Plant.Position, sep="_"))

# exploratory stats
ant_aphid_dissect_plant_info <- merge(plant_info_ant_aphid_experiment, sum_compact_ant_aphid_dissect, by="plant_code", all.x=TRUE)
ant_aphid_dissect_plant_info[ ,7:8][is.na(ant_aphid_dissect_plant_info[ ,7:8])] <- 0

# focus on LTF abundance
with(ant_aphid_dissect_plant_info, interaction.plot(Aphids.or.No.Aphids, Genotype, LTF_abund, col=1:10, lwd=2))
with(ant_aphid_dissect_plant_info, interaction.plot(Distant.to.Ant.Mound, Genotype, LTF_abund, col=1:10, lwd=2))
with(ant_aphid_dissect_plant_info, interaction.plot(Distant.to.Ant.Mound, Aphids.or.No.Aphids, LTF_abund, col=1:10, lwd=2))
plot(ant_aphid_dissect_plant_info$LTF_abund ~ ant_aphid_dissect_plant_info$Genotype)


aov_antaphid_LTF <- aov(log(LTF_abund + 1) ~ Genotype * Aphids.or.No.Aphids * Distant.to.Ant.Mound + Error(as.factor(Block)/Distant.to.Ant.Mound), data=ant_aphid_dissect_plant_info, qr=TRUE)
summary(aov_antaphid_LTF)
plot(fitted(aov_antaphid_LTF[[4]]), studres(aov_antaphid_LTF[[4]])) 
abline(h=0, lty=2)

# focus on gallery mine abundance
with(ant_aphid_dissect_plant_info, interaction.plot(Aphids.or.No.Aphids, Genotype, fresh_gallery_mine, col=1:10, lwd=2))
with(ant_aphid_dissect_plant_info, interaction.plot(Distant.to.Ant.Mound, Genotype, fresh_gallery_mine, col=1:10, lwd=2))
with(ant_aphid_dissect_plant_info, interaction.plot(Distant.to.Ant.Mound, Aphids.or.No.Aphids, fresh_gallery_mine, col=1:10, lwd=2))

aov_antaphid_gallery <- aov(log(fresh_gallery_mine + 1) ~ Genotype * Aphids.or.No.Aphids * Distant.to.Ant.Mound + Error(as.factor(Block)/Distant.to.Ant.Mound), data=ant_aphid_dissect_plant_info, qr=TRUE)
summary(aov_antaphid_gallery)
plot(fitted(aov_antaphid_gallery[[4]]), studres(aov_antaphid_gallery[[4]])) # residuals are "fanned" and look horrible.
abline(h=0, lty=2)