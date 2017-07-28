## load required libraries ----
library(dplyr)
library(tidyr)
library(vegan) # for rarefaction

## useful functions ----
left_join_NA_to_0 <- function(x, y, ...) {
  left_join(x = x, y = y, by = ...) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0)))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

## upload and manage wind experiment data ----

## wind visual data 2012
wind.2012.vis <- read.csv("manage_raw_data/raw_data/arthropod_wind_visual_data_2012.csv") %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block),
         Year = 2012, Date = "July22") %>%
  select(Year, Date, Block, Wind.Exposure, Genotype, plant_code, 
         Dead = Dead., 
         ant_F_obscuripes = Red.ants,
         ant_black = Black.Ants,
         spider_BY = Spiders.black.with.yellow.legs,
         spider_unk = Spiders.Other, # don't have good morphospecies data on these spiders
         leafhopper_C_reductus = leaf_hoppers_blue_brown,
         leafhopper_camo = Leaf.Hoppers.Camo,
         leafhopper_YK = yellow_knife_leaf_hopper,
         leafhopper_nymph_unk = LH.nymph.other,
         aphid_LG = Light.Green.Aphids,
         # omitting twisty_gall because I was never able to confirm whether this was a galling insect.
         gall_R_salicisbattatus = stem_gall,
         gall_Iteomyia = volcano_gall,
         gall_Aculus = mite_gall,
         psyllid = psyllidae,
         sawfly_larva = sawflies,
         aphid_Tuberolachnus = giant_willow_aphid,
         aphid_Aphis = aphid_species_used_in_ant_aphid_experiment) %>%
  filter(Dead < 1) %>% # only retain living plants
  select(-Dead, -Date)

glimpse(wind.2012.vis)

## load leaf miner data for wind 2012 
wind_leaf_miners <- read.csv("manage_raw_data/raw_data/wind_dissections_leaf_miner_summer_2012.csv", stringsAsFactors=F, skip = 1) %>% # this data set focuses on the leaf mining moth, Caloptilia pruniella
  tbl_df() %>%
  mutate(plant_code = paste(Treatment, Block, Plant, sep="_"),
         dead_Caloptilia = gallery_mine_dead + gallery_mine_fresh_but_stiff + gallery_freshly_dead + tent_mine_stiff + tent_mine_freshly_dead + tent_mine_dead + tent_mine_fresh_but_stiff + LTF_dead + LEF_dead_larva,
         live_Caloptilia = gallery_mine_living + tent_mine_living_larva + LEF_live_adult,
         all_Caloptilia = live_Caloptilia + dead_Caloptilia,
         ptoid_sp_323_all = gallery_mine_323 + tent_mine_323, # sp_323 (larva, unknown what adult is); pics of parasitoid sp 323 in iPhoto database (unknown type though)
         ptoid_sp_329 = LTF_329,
         ptoid_LEF_eulophid_cocoons = LEF_eulophid_cocoons, # eulophid cocoons found on leaf edge fold stage of Caloptilia
         
         # Caloptilia by life-stage data
         gallery_mine_Caloptilia = gallery_mine_missing + gallery_mine_dead + gallery_mine_fresh_but_stiff + gallery_freshly_dead + gallery_mine_living + gallery_mine_323 + gallery_mine_parasitized_unknown + gallery_mine_other_category,
         fresh_gallery_mine = gallery_mine_fresh_but_stiff + gallery_freshly_dead + gallery_mine_living,
         tent_mine_Caloptilia = tent_mine_353 + tent_mine_parasitized_unknown + tent_mine_missing + tent_mine_other + tent_mine_stiff + tent_mine_empty_with_frass + tent_mine_freshly_dead + tent_mine_dead + tent_mine_empty + tent_mine_fresh_but_stiff + tent_mine_living_larva + tent_mine_323,
         LTF_Caloptilia = LTF_empty + LTF_empty_with_frass + LTF_327 + LTF_329 + LTF_other + LTF_dead + LTF_347,
         LEF_Caloptilia = LEF_intact + LEF_eulophid_cocoons + LEF_other_parasitoid_. + LEF_dead_larva + LEF_emerged + LEF_empty + LEF_live_adult) %>%
  
      # subset of data I'm interested in. I decided to only retain LTF_Caloptilia because I think is the most reliably quantified and good indicator of damage cause by Caloptilia across years and experiments. I also retain ptoid_sp_329 since it was the parasitoid associated with this life-stage. I may not end up using this data though.
  select(plant_code, 
         LTF_Caloptilia) %>% # indicator of amount of damage and this will be a consistent measure of Caloptilia herbivory across years and experiments.
         # ptoid_LTF_Braconid = ptoid_sp_329) not including because this came from dissection data and not visual surveys.
  group_by(plant_code) %>%
  summarise_each(funs(sum))

## load leaf tie data for wind 2012 
wind_leaf_ties <- read.csv("manage_raw_data/raw_data/wind_dissections_leaf_tie_summer_2012.csv", skip = 1, stringsAsFactors=F) %>%
  tbl_df() %>%
  mutate(plant_code = paste(Treatment, Block, Plant, sep="_"),
         live_tortricid = LT_sp_328 + live_320_LT + LB_moth_pupa_LT + living_319_LT + living_320_LT + living_320_321_LT + living_321_LT + living_322 + letter_T_moth_LT, # live tortricid moths (sp. 319 - 322 and sp. 328). Assuming they are all the same species.  Evidence suggests this, because I have seen some "intermediates" mophotypes between life stages
         dead_tortricid = LT_dead_319 + LT_dead_320 + dead_fresh_320_321,
         all_tortricid = live_tortricid + dead_tortricid,
         ptoid_sp_327_LT.1 = sp_327_LT.1) %>% # ichneumonid sp. 327 associated with totricid moth leaf bundles.
  
  # didn't select any spider data to avoid double-counting with visual surveys. Also, didn't count the other morphospecies because I'm unsure about them.
  # to stay consistent with the LTF data, I'm only retaining 'leaf_dam_LER' since this is a consistent indicator of damage caused by this tortricid moth across years and experiments.
  select(plant_code, 
         leaftier_Tortricid = leaf_dam_LER,
         tentmine_Phyllonorycter = other_tent_mine_pupa) %>% # damage caused by mid-level instars
  group_by(plant_code) %>%
  summarise_each(funs(sum))

## merge wind 2012 arthropod datasets together 
wind.2012.arth.df <- left_join_NA_to_0(wind.2012.vis, wind_leaf_miners) %>%
  left_join_NA_to_0(wind_leaf_ties)

## upload and manage wind data for 2013 (visual surveys only)
May10 <- read.csv('manage_raw_data/raw_data/Wind_Arthropod_data_survey_1_May_10_2013.csv', skip = 4, stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(plant_code = paste(Wind.Exposure, Block, Plant.Position, sep = "_"), Year = 2013, Date = "May10")

June27 <- read.csv('manage_raw_data/raw_data/Wind_Arthropod_data_survey_2_June_27_2013.csv', skip = 1, stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(plant_code = paste(Wind.Exposure, Block, Plant.Position, sep = "_"), Year = 2013, Date = "June27")

July10 <- read.csv('manage_raw_data/raw_data/wind_arthropod_data_survey_3_July_10_2013.csv', skip = 1, stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(plant_code = paste(Wind.Exposure, Block, Plant.Position, sep = "_"), Year = 2013, Date = "July10")

wind.2013.vis <- bind_rows(May10, June27, July10) 
wind.2013.vis[is.na(wind.2013.vis)] <- 0

wind.2013.vis.df <- wind.2013.vis %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure),
         Genotype = as.factor(Genotype),
         plant_code = as.factor(plant_code),
         leafhopper_C_reductus = Colladonus_reductus + LH_nymph_black_banded_abdomen,
         tentmine_Phyllonorycter = phyllonorycter_tent_mine + phyllonoricter_tent_mine,
         gall_Pontania = gall_Pontania_leaf + petiole_gall_or_Iteomyia_nobumps + leaf_gall_without_points,
         spider_Larionoides = spider_larionoides_sp + larinoides_spider + unknown_spider_sp_hairy, # added unknown hairy spider sp. because this matches Larionoides as a potential sp.
         spider_CS = clear_small_spider + clear_spider,
         red_scale = red_scale_insect + red_scale_insects,
         caterpillar_unk = unknown_tiny_moth_larva_silking_topofleaf + striped.moth.larva + unknown_moth_silking_underside_of_leaf,
         LTF_Caloptilia = Calop_LTF + Calop_LTF_old) %>% # this is the most reliable indicator of Caloptilia damage. Omitting Calop_tent_mine, Calop_pre_tentmine, Calop_cocoon, and Calop_tent_mine_old.
  select(Year, Date, Block, Wind.Exposure, Genotype, plant_code, Dead,
         LTF_Caloptilia,
         tentmine_Phyllonorycter,
         leaftier_Tortricid = leaf_edge_roll_and_silk, # this is the most reliable indicator of Tortricid moth damage. Omitting Leaf_bundle.
         gall_R_rigidae = gall_shoot_top,
         gall_R_salicisbrassicoides = gall_bud,
         gall_Pontania,
         gall_Aculus = mite_leaf_gall,
         leafhopper_C_reductus,
         leafhopper_green = green_leafhopper,
         leafhopper_unk = unknown_leafhopper_adult,
         sawfly_larva,
         caterpillar_looper = looper_rattlesnake,
         caterpillar_LB = new_leaf_bundler, # not the same as the Tortricid moth that typically causes leaf bundles.
         caterpillar_unk,
         red_scale,
         psyllid = adult_psyllid_or_aphid,
         grasshopper,
         stinkbug,
         ant_F_obscuripes = red_ants,
         spider_Theridion = spider_Theridion_sp.,
         spider_BY = spider_black_with_yellow_legs,
         spider_NW = nursery_web_spider,
         spider_Tetragnathid,
         spider_CS,
         spider_Larionoides) %>%
   filter(Dead < 1) # only retain living plants

wind.2013.vis.df.max <- wind.2013.vis.df %>%
  group_by(Block, Wind.Exposure, Genotype, plant_code) %>%
  summarise_each(funs(max)) %>%
  select(Year, Block:plant_code, LTF_Caloptilia:spider_Larionoides)

## merge wind arthropod 2012 and 2013 datasets
wind.arth.df <- bind_rows(wind.2012.arth.df, wind.2013.vis.df.max) %>%
  mutate(plant_code = as.factor(plant_code))
wind.arth.df[is.na(wind.arth.df)] <- 0
glimpse(wind.arth.df)

wind.arth.df <- wind.arth.df %>% separate(plant_code, into = c("Block.tmp","Wind.Exposure.tmp","Plant_Position")) %>% mutate(treat.tmp = ifelse(Wind.Exposure == "Exposed", "E", "U"), plant_ID = paste(Block, treat.tmp, Genotype, Plant_Position, sep = "")) %>% select(Year, Block, Wind.Exposure, Genotype, plant_ID, ant_F_obscuripes:spider_Larionoides)

# create data on total abundance, richness, and rarefied richness
wind.arth.df$total.abund <- rowSums(
  select(wind.arth.df, ant_F_obscuripes:spider_Larionoides))
wind.arth.df$total.rich <- rowSums(
  select(wind.arth.df, ant_F_obscuripes:spider_Larionoides) > 0)
wind.arth.df$total.rarerich <- rarefy(select(wind.arth.df, ant_F_obscuripes:spider_Larionoides), 2) - 1

na.rarerich <- which(wind.arth.df$total.abund < 2) # can only test rarefied richness when at least 2 individuals were sampled.
wind.arth.df$total.rarerich[na.rarerich] <- NA

wind.arth.df.families <- wind.arth.df %>%
  mutate(Gracilliaridae_miner = LTF_Caloptilia + tentmine_Phyllonorycter,
         Tortricidiae_leaftier = leaftier_Tortricid,
         Cecidomyiidae_gall = gall_R_rigidae + gall_R_salicisbrassicoides + gall_R_salicisbattatus + gall_Iteomyia,
         Tenthredinidae_sawfly = gall_Pontania + sawfly_larva,
         Cicadellidae = leafhopper_C_reductus + leafhopper_green + leafhopper_unk + leafhopper_camo + leafhopper_YK + leafhopper_nymph_unk,
         Psyllidae = psyllid,
         Eriophyidae = gall_Aculus,
         caterpillar_other = caterpillar_looper + caterpillar_LB + caterpillar_unk,
         Coccoidea = red_scale,
         Orthoptera = grasshopper,
         Aphididae = aphid_Tuberolachnus + aphid_Aphis + aphid_LG,
         Pentatomidae = stinkbug,
         Formica_ant = ant_F_obscuripes + ant_black,
         Spider = spider_Theridion + spider_BY + spider_NW + spider_Tetragnathid + spider_CS + spider_Larionoides + spider_unk)

write.csv(wind.arth.df.families, 'final_data/wind_arthropod_df.csv')

## upload and manage ant-aphid data for 2012 ----

# ant-aphid experimental info
#aa.info <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Ant-Aphid_Experiment_Setup.csv") %>%
 # tbl_df() %>%
  #mutate(plant_code = paste(Block, Plant.Position, sep = "_")) %>%
  #select(Block, Aphids.or.No.Aphids, Distant.to.Ant.Mound, Plant.Position, Genotype, plant_code)

# ant-aphid visual survey data 2012. Starting with just ants and focal aphids
aa.2012.vis <- read.csv("manage_raw_data/raw_data/Ant-Aphid_Data_v2.csv", stringsAsFactors=F, skip = 1) %>%
  tbl_df() %>%
  mutate(X = as.factor(1:1800), # unique ID for each data point
         #plant_code = paste(Block, Plant.Position, sep = "_"),
         Plot_code = paste(Block, Distance.to.Ant.Mound, sep = "_"),
         plant_ID = paste(Block, Distance.to.Ant.Mound, Aphid.Treatment, Genotype, Plant.Position, sep = "_"),
         Block = as.factor(Block),
         aphid_Aphis = Green.Aphids + Orange.Aphids,
         leafhopper = Leaf.Hoppers.Brown.Blue + LH.nymph.blk.ylw + Leaf.Hoppers.Camo + leafhopper.brownish + leafhopper.grey + leafhopper.garden.sp7 + scuttler, # note that this may underestimate richness estimates slightly since these ambiguous categories are being combined, but it will give us a conservative estimate.
         froghopper = Frog.Hopper.brown + brown.megatron + Spittlebugs,
         aphid_Tuberolachnus = Giant.Willow.Aphids + giant.willow.aphids,
         psyllid = psyillid_sp13 + psyllid) %>%
  select(Date_rel = relative.Date,
         Block, 
         Genotype,
         Ant.mound.dist = Distance.to.Ant.Mound,
         Aphid.treatment = Aphid.Treatment, 
         #plant_code, 
         plant_ID, X, 
         Dead, Plot_code,
         
         # only selected species/morphospecies in which I'm confident in their identification, trophic position, and ones that are not transients (e.g adult flies, Chalcid.Parasitoids)
         aphid_Aphis,
         aphid_Tuberolachnus,
         aphid_LG = Light.Green.Aphids, # note that we only had data available on this herbivore from the last survey, although the notes suggest it started to show up as early as July 14.
         psyllid,
         leafhopper,
         froghopper,
         spiders = Spiders,
         sawfly_larva = sawfly.larva,
         ant_F_obscuripes = Red.Ants, 
         ant_black = Black.Ants,
         syrphid_larva = syriphid.larva,
         gall_R_salicisbattatus = stem.galls,
         grasshopper,
         leaftier_Tortricid = Leaf.Edge.Rollers, # omitting Leaf.Bundlers so I have one consistent indicator of Tortricid moth damage/abundance. Note that no data was taking on leaf edge rollers for last survey, July 20. I think this was because I was planning to take data on them from my leaf dissections, but I never did this (apparently)
         LTF_Caloptilia = LTF.abund) %>% # Note that no data was taking on LTF.abund for last survey, July 20. I think this was because I was taking data on them from my leaf dissections. 
  filter(Dead < 1) # only retain living plants

aa.2012.vis.max <- aa.2012.vis %>%
  select(-X, -Date_rel, -Dead) %>%
  group_by(Block, Genotype, Ant.mound.dist, Aphid.treatment, Plot_code, plant_ID) %>%
  summarise_each(funs(max.narm = max(., na.rm = TRUE))) %>% # take the maximum number of individuals observed on a plant over the entire survey period as a conservative estimate of abundance over the entire growing season. 
  ungroup()

aphid.growth.df <- aa.2012.vis %>%
  select(Date_rel:plant_ID, aphid_Aphis) %>%
  spread(Date_rel, aphid_Aphis) %>% 
  filter(Aphid.treatment == "aphid") %>%
  mutate(start.7 = 5,
         start.14 = ifelse(.$'7' > 5, .$'7', 5),
         start.25 = ifelse(.$'14' > 5, .$'14', 5),
         start.35 = ifelse(.$'25' > 5, .$'25', 5),
         start.45 = ifelse(.$'35' > 5, .$'35', 5),
         start.51 = ifelse(.$'45' > 5, .$'45', 5),
         week1.growth = (log(.$'7' + 1) - log(start.7))/(7-0),
         week2.growth = (log(.$'14' + 1) - log(start.14))/(14-7),
         week3.growth = (log(.$'25' + 1) - log(start.25))/(25-14),
         week4.growth = (log(.$'35' + 1) - log(start.35))/(35-25),
         week5.growth = (log(.$'45' + 1) - log(start.45))/(45-35),
         week6.growth = (log(.$'51' + 1) - log(start.51))/(51-45)) %>%
  select(plant_ID, week1.growth:week6.growth) %>%
  gather(Week, Aphid.growth.rate, week1.growth:week6.growth) %>%
  mutate(Date_rel = ifelse(Week == "week1.growth", 7, 
                                ifelse(Week == "week2.growth", 14,
                                       ifelse(Week == "week3.growth", 25, 
                                              ifelse(Week == "week4.growth", 35,
                                                     ifelse(Week == "week5.growth", 45, 
                                                            ifelse(Week == "week6.growth", 51, NA))))))) %>%
  select(Date_rel, plant_ID, 
         Aphis.growth.rate = Aphid.growth.rate) %>%
  left_join(., select(aa.2012.vis, Date_rel:Dead)) %>%
  select(Date_rel, Block:Aphid.treatment, plant_ID, X, Aphis.growth.rate)

## ant-aphid leaf miner dissection data 2012 
# note that there is no leaf tier dissection data.
aa.leafminers.2012 <- read.csv("manage_raw_data/raw_data/ant_aphid_arthropod_dissection_summer_2012.csv", stringsAsFactor=F, skip=1) %>%
  tbl_df() %>%
  mutate(plant_code = paste(Block, Plant, sep="_"),
         fresh_gallery_mine = gallery_mine_living + gallery_mine_fresh_but_stiff + gallery_mine_freshly_dead, # composite variable for living or freshly dead gallery mines (late season cohort abundance)
         LTF_Caloptilia_dissect = LTF_329 + LTF_empty_with_frass + LTF_dead + LTF_other + LTF_empty) %>% # composite variable for any evidence of a leaf tip fold (early season cohort abundance)
  
  # only retained LTF_Caloptilia because this is the one reliable and consistent indicator of Caloptilia damage across years and experiments
  select(plant_code, LTF_Caloptilia_dissect) %>%
  group_by(plant_code) %>%
  summarise_each(funs(sum)) 

## merge ant-aphid datasets together 
aa.arth.df <- aa.2012.vis.max #left_join_NA_to_0(aa.2012.vis.max, aa.leafminers.2012) # decided to not include leaf miner dissection data because I feel the counts of LTF from the field were more accurate.
glimpse(aa.arth.df)

aa.arth.names <- colnames(select(aa.arth.df, aphid_Tuberolachnus:sawfly_larva, ant_black:LTF_Caloptilia))

# generate new columns for total abundance and richness that do not include aphid_Aphis or ant_F_obscuripes
aa.arth.df$total.abund <- rowSums(aa.arth.df[ ,aa.arth.names])
aa.arth.df$total.rich <- rowSums(aa.arth.df[ ,aa.arth.names] > 0)
aa.arth.df$total.rarerich <- rarefy(aa.arth.df[ ,aa.arth.names], 2) - 1
aa.na.rarerich <- which(aa.arth.df$total.abund < 2) # can only test rarefied richness when at least 2 individuals were sampled.
aa.arth.df$total.rarerich[aa.na.rarerich] <- NA
#aa.herb.abund <- rowSums(

aa.arth.df.families <- aa.arth.df %>%
  mutate(Gracilliaridae_miner = LTF_Caloptilia, # no tentmine_Phyllonorycter
         Tortricidiae_leaftier = leaftier_Tortricid,
         Cecidomyiidae_gall = gall_R_salicisbattatus, # no gall_R_rigidae or gall_R_salicisbrassicoides 
         Tenthredinidae_sawfly = sawfly_larva, # no gall_Pontania 
         Cicadellidae = leafhopper,
         Cercopidae = froghopper,
         Psyllidae = psyllid,
         #Eriophyidae = gall_Aculus,
         #caterpillar_other = caterpillar_looper + caterpillar_LB + caterpillar_unk,
         #Coccoidea = red_scale,
         Orthoptera = grasshopper,
         Aphididae = aphid_Tuberolachnus + aphid_LG, # not-including aphid_Aphis since they were part of the experiment
         Syrphidae = syrphid_larva,
         #Pentatomidae = stinkbug,
         Formica_ant = ant_black, # not including ant_F_obscuripes since they were part of the experiment
         Spider = spiders)
  
write.csv(aa.arth.df.families, 'final_data/ant_aphid_arthropod_df.csv')
write.csv(aphid.growth.df, 'final_data/ant_aphid_Aphis_popgrowth_df.csv')
