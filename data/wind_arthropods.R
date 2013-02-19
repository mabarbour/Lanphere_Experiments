# This code analyzes arthropod community data collected from my willow genotype-by-wind exposure experiment at Lanphere Dunes, California (May 2012 - ______)

######## Setup working directory and load required libraries
setwd("~/Documents/Lanphere_Experiments/data")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
library("bipartite")

######## Upload arthropod community data and do some data management
# upload data
wind_arthropod_visual_survey_data_summer_2012 <- read.csv("~/Documents/Lanphere_Experiments/data/wind_arthropod_data_summer_2012.csv", skip = 1, header=TRUE)

# change "parasitoid cocoon on aphids" column to "lacewing eggs" because that is what they actually area
wind_arthropod_visual_survey_data_summer_2012$lacewing_eggs <- wind_arthropod_visual_survey_data_summer_2012$Parisatoid.cocoons.on.Aphids

# collapse previously separate morphospecies
wind_arthropod_visual_survey_data_summer_2012$psyllidae <- wind_arthropod_visual_survey_data_summer_2012$psyllid + wind_arthropod_visual_survey_data_summer_2012$psyllid_nymph # psyllid adult and nymphs together
wind_arthropod_visual_survey_data_summer_2012$sawflies <- wind_arthropod_visual_survey_data_summer_2012$striped_sawfly + wind_arthropod_visual_survey_data_summer_2012$sawfly_larva # all sawfly larva together
wind_arthropod_visual_survey_data_summer_2012$aphid_species_used_in_ant_aphid_experiment <- wind_arthropod_visual_survey_data_summer_2012$Green.Aphids + wind_arthropod_visual_survey_data_summer_2012$Orange.Aphids

# gather up known predator morphospecies (currently excluding lacewing_eggs)
predator_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(Black.Ants,Red.ants,Spiders.black.with.yellow.legs, Spiders.Other))

# gather up known herbivore morphospecies
herbivore_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(Leaf.Hoppers.Blue.Brown, Leaf.Hoppers.Camo, Light.Green.Aphids, stem_gall, volcano_gall, mite_gall, yellow_knife_leaf_hopper, twisty_gall, giant_willow_aphid, psyllidae, sawflies, LH.nymph.blk.ylw, LH.nymph.other, aphid_species_used_in_ant_aphid_experiment)) # currently leaf hopper nymphs are not separated

# data frame for relevant info on each plant
plant_info_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(Genotype, Block, Wind.Exposure))

# put all of the relevant data back together (currently excluding any counts of eggs, other than lacewings...)
arthropod_visual_data_2012 <- cbind.data.frame(plant_info_data_frame, predator_data_frame, herbivore_data_frame)

avg_arthropod_by_treatment <- aggregate(arthropod_visual_data_2012[ ,4:21], by = list(arthropod_visual_data_2012$Genotype, arthropod_visual_data_2012$Wind.Exposure), mean)
avg_arthropod_by_treatment$Group_combine <- paste(avg_arthropod_by_treatment$Group.1, avg_arthropod_by_treatment$Group.2, sep="_")
rownames(avg_arthropod_by_treatment) <- avg_arthropod_by_treatment$Group_combine

# some new summary variables for each plant
#arthropod_visual_data_2012$species_richness_arthropods <- specnumber(community_matrix_arthropods)
#arthropod_visual_data_2012$abundance_arthropods <- rowSums(community_matrix_arthropods)
#dissimilarity_arthropod <- vegdist(community_matrix_arthropods, method = "bray")

################# Exploratory Data Analysis
visweb(avg_arthropod_by_treatment[ ,3:20], preynames=TRUE) # light green aphids are the most abundant arthropods and appeared to be associated with wind exposed sites.

visweb(avg_arthropod_by_treatment[ ,c(3:8, 10:20)], type="diagonal", preynames=TRUE) # exclude ligh green aphids
