# This script organizes my leaf miner and leaf bundler dissection data collected in July of 2012 for my wind experiment at Lanphere Dunes, California.  

######### set working directory and load any required libraries
setwd("~/Documents/Lanphere_Experiments/data")

######### upload and manage data
# load leaf miner data
wind_leaf_miners <- read.csv("~/Documents/Lanphere_Experiments/data/wind_dissections_leaf_miner_summer_2012.csv", stringsAsFactors=F, skip = 1) # this data set focuses on the leaf mining moth, Caloptilia pruniella

# load leaf tie data
wind_leaf_ties <- read.csv("~/Documents/Lanphere_Experiments/data/wind_dissections_leaf_tie_summer_2012.csv", skip = 1, stringsAsFactors=F)

###### manage data and create new variables
# create variable unique ID for each plant in both data sets
wind_leaf_miners$plant_code <- with(wind_leaf_miners, paste(Treatment, Block, Plant, sep="_")) 

wind_leaf_ties$plant_code <- with(wind_leaf_ties, paste(Treatment, Block, Plant, sep="_"))

# dead Caloptilia at all stages of life cycle
wind_leaf_miners$dead_Caloptilia <- with(wind_leaf_miners, gallery_mine_dead + gallery_mine_fresh_but_stiff + gallery_freshly_dead + tent_mine_stiff + tent_mine_freshly_dead + tent_mine_dead + tent_mine_fresh_but_stiff + LTF_dead + LEF_dead_larva)

# live Caloptilia at all stages of life cycle
wind_leaf_miners$live_Caloptilia <- with(wind_leaf_miners, gallery_mine_living + tent_mine_living_larva + LEF_live_adult)

# live and dead Caloptilia at all stages of life cycle
wind_leaf_miners$all_Caloptilia <- wind_leaf_miners$live_Caloptilia + wind_leaf_miners$dead_Caloptilia

# live tortricid moths (sp. 319 - 322 and sp. 328). Assuming they are all the same species.  Evidence suggests this, because I have seen some "intermediates" mophotypes between life stages
wind_leaf_ties$live_tortricid <- with(wind_leaf_ties, LT_sp_328 + live_320_LT + LB_moth_pupa_LT + living_319_LT + living_320_LT + living_320_321_LT + living_321_LT + living_322 + letter_T_moth_LT)

# dead larva of tortricid species abova
wind_leaf_ties$dead_tortricid <- with(wind_leaf_ties, LT_dead_319 + LT_dead_320 + dead_fresh_320_321)

# live and dead tortricids at all stages
wind_leaf_ties$all_tortricid <- wind_leaf_ties$live_tortricid + wind_leaf_ties$dead_tortricid

## Parasitoids. May need to add more after GRC once I have a chance to double check everything. Will collpase together ones that appear on both tortricid moths and Caloptilia.

# sp_323 (larva, unknown what adult is)
wind_leaf_miners$sp_323_all <- with(wind_leaf_miners, gallery_mine_323 + tent_mine_323) # pics of parasitoid sp 323 in iPhoto database (unknown type though)

# sp_329
wind_leaf_miners$LTF_329 # no additions necesary

# eulophid cocoons found on leaf edge fold stage of Caloptilia
wind_leaf_miners$LEF_eulophid_cocoons # no additions necessary

# ichneumonid sp. 327 associated with totricid moth leaf bundles.
wind_leaf_ties$sp_327_LT.1

# wind_leaf_ties$white_parasitoid # assuming this is the one that makes the white papery cocoons...omitting for now

######### For merging into the arthropod survey data, I will make each data set of a more manageable size for examining patterns of arthropod abundance and richness
sub_wind_miners <- subset(wind_leaf_miners, select = c("live_Caloptilia", "all_Caloptilia", "LTF_329", "LEF_eulophid_cocoons", "sp_323_all", "plant_code"))
plant_sum_wind_miner_abund <- rowsum(sub_wind_miners[ ,1:5], sub_wind_miners$plant_code)
ID_plant_miners <- rownames(plant_sum_wind_miner_abund)
compact_wind_miners <- cbind(plant_sum_wind_miner_abund, ID_plant_miners)

sub_wind_leaf_ties <- subset(wind_leaf_ties, select = c("live_tortricid", "all_tortricid", "sp_327_LT.1", "plant_code"))
plant_sum_wind_ties_abund <- rowsum(sub_wind_leaf_ties[ ,1:3], sub_wind_leaf_ties$plant_code)
ID_plant_ties <- rownames(plant_sum_wind_ties_abund)
compact_wind_ties <- cbind(plant_sum_wind_ties_abund, ID_plant_ties)

##### Caloptilia by life stage dataset
wind_leaf_miners$gallery_mine_Caloptilia <- with(wind_leaf_miners, gallery_mine_missing + gallery_mine_dead + gallery_mine_fresh_but_stiff + gallery_freshly_dead + gallery_mine_living + gallery_mine_323 + gallery_mine_parasitized_unknown + gallery_mine_other_category)

wind_leaf_miners$fresh_gallery_mine <- with(wind_leaf_miners, gallery_mine_fresh_but_stiff + gallery_freshly_dead + gallery_mine_living)

wind_leaf_miners$tent_mine_Caloptilia <- with(wind_leaf_miners, tent_mine_353 + tent_mine_parasitized_unknown + tent_mine_missing + tent_mine_other + tent_mine_stiff + tent_mine_empty_with_frass + tent_mine_freshly_dead + tent_mine_dead + tent_mine_empty + tent_mine_fresh_but_stiff + tent_mine_living_larva + tent_mine_323)

wind_leaf_miners$LTF_Caloptilia <- with(wind_leaf_miners, LTF_empty + LTF_empty_with_frass + LTF_327 + LTF_329 + LTF_other + LTF_dead + LTF_347)

wind_leaf_miners$LEF_Caloptilia <- with(wind_leaf_miners, LEF_intact + LEF_eulophid_cocoons + LEF_other_parasitoid_. + LEF_dead_larva + LEF_emerged + LEF_empty + LEF_live_adult)

# bring together data set
sub_Caloptilia <- subset(wind_leaf_miners, select = c("gallery_mine_Caloptilia", "fresh_gallery_mine", "tent_mine_Caloptilia", "LTF_Caloptilia", "LEF_Caloptilia", "plant_code"))
plant_sum_Caloptilia <- rowsum(sub_Caloptilia[ ,1:5], sub_Caloptilia$plant_code)
ID_plant_miners <- rownames(plant_sum_Caloptilia)
sub_Caloptilia <- cbind(ID_plant_miners, plant_sum_Caloptilia)