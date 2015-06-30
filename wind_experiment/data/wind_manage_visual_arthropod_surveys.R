# This code analyzes arthropod community data collected from my willow genotype-by-wind exposure experiment at Lanphere Dunes, California (May 2012 - ______)

######## Setup working directory and load required libraries
#setwd("~/Documents/Lanphere_Experiments/data")

######## Upload arthropod community data and do some data management
# upload data
wind_arthropod_visual_survey_data_summer_2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_arthropod_data_summer_2012.csv", skip = 1, header=TRUE)

# change "parasitoid cocoon on aphids" column to "lacewing eggs" because that is what they actually area
wind_arthropod_visual_survey_data_summer_2012$lacewing_eggs <- wind_arthropod_visual_survey_data_summer_2012$Parisatoid.cocoons.on.Aphids

# collapse previously separate morphospecies
wind_arthropod_visual_survey_data_summer_2012$psyllidae <- wind_arthropod_visual_survey_data_summer_2012$psyllid + wind_arthropod_visual_survey_data_summer_2012$psyllid_nymph # psyllid adult and nymphs together
wind_arthropod_visual_survey_data_summer_2012$sawflies <- wind_arthropod_visual_survey_data_summer_2012$striped_sawfly + wind_arthropod_visual_survey_data_summer_2012$sawfly_larva # all sawfly larva together
wind_arthropod_visual_survey_data_summer_2012$aphid_species_used_in_ant_aphid_experiment <- wind_arthropod_visual_survey_data_summer_2012$Green.Aphids + wind_arthropod_visual_survey_data_summer_2012$Orange.Aphids
wind_arthropod_visual_survey_data_summer_2012$leaf_hoppers_blue_brown <- wind_arthropod_visual_survey_data_summer_2012$Leaf.Hoppers.Blue.Brown + wind_arthropod_visual_survey_data_summer_2012$LH.nymph.blk.ylw # big assumption right now that black and yellow leaf hopper nymphs become blue/brown leaf hoppers (based on temporal observations in the field).

# gather up known predator morphospecies (currently excluding lacewing_eggs)
visual_predator_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(Black.Ants,Red.ants,Spiders.black.with.yellow.legs, Spiders.Other))

# gather up known herbivore morphospecies
visual_herbivore_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(leaf_hoppers_blue_brown, Leaf.Hoppers.Camo, Light.Green.Aphids, stem_gall, volcano_gall, mite_gall, yellow_knife_leaf_hopper, twisty_gall, giant_willow_aphid, psyllidae, sawflies, LH.nymph.other, aphid_species_used_in_ant_aphid_experiment)) # currently leaf hopper nymphs are not separated

# create unique variable ID for each plant that matches arthropod dissection data set. Also, create separate data frame with just relevant plant information
wind_arthropod_visual_survey_data_summer_2012$plant_code <- with(wind_arthropod_visual_survey_data_summer_2012, paste(Wind.Exposure, Block, Plant.Position, sep="_"))

plant_info_data_frame <- with(wind_arthropod_visual_survey_data_summer_2012, cbind.data.frame(Genotype, Block, Wind.Exposure, Plant.Position, plant_code, Dead.))

# put all of the relevant data back together (currently excluding any counts of eggs)
arthropod_wind_visual_data_2012 <- cbind.data.frame(plant_info_data_frame, visual_predator_data_frame, visual_herbivore_data_frame)


write.csv(arthropod_wind_visual_data_2012, "~/Documents/Lanphere_Experiments/wind_experiment/data/arthropod_wind_visual_data_2012.csv")
