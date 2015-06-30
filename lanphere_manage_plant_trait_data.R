## Description: This script organizes plant architecture and leaf quality data collected in the summers of 2012 and 2013 as part of the Lanphere dune Experiments.

## load required packages
library(dplyr)
library(tidyr)

#### upload and manage wind and ant-aphid data

## upload wind plant info
wind_plant_info <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_info.csv") %>%
  tbl_df() %>%
  mutate(block = as.factor(block)) %>%
  select(block, treatment, genotype, plant.position, plant.id)

## wind 2012 - plant traits

# plant architecture 
wind.arch2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_traits_field_summer_2012.csv") %>%
  tbl_df() %>%
  mutate(treatment = ifelse(Wind.Exposure == "Exposed", "E", "U"),
         plant.id = paste0(Block,treatment,Plant.Position))
wind.arch2012$mature.shoot.total.length <- rowSums(select(wind.arch2012, X1.shoot.M:X10.shoot.M))
wind.arch2012$all.shoot.total.length <- rowSums(select(wind.arch2012, X1.shoot.M:X6.shoot.I))
wind.arch2012$mature.shoot.count <- rowSums(select(wind.arch2012, X1.shoot.M:X10.shoot.M) > 0)
wind.arch2012$mature.shoot.avg.length <- rowMeans(select(wind.arch2012, X1.shoot.M:X10.shoot.M))
wind.arch2012$all.shoot.avg.length <- rowMeans(select(wind.arch2012, X1.shoot.M:X6.shoot.I))
wind.arch2012 <- wind.arch2012 %>%
  select(plant.id, Dead = Dead., Height = Max.Height, 
         leaf.count = Mature.Green.Leaves, 
         mature.shoot.total.length, all.shoot.total.length, 
         mature.shoot.avg.length, 
         all.shoot.avg.length, mature.shoot.count, 
         immature.shoot.count = Shoots.sprouting.at.base) %>%
  mutate(all.shoot.count = mature.shoot.count + immature.shoot.count)

# wind leaf quality data
wind.WetLeafWts.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Wind leaf trait data.csv") %>%
  tbl_df() %>%
  .[-which(.$Collection.Number %in% 
             names(which(table(.$Collection.Number) > 1))), ] %>% # removes all duplicates
  mutate(Collection.Number = as.character(Collection.Number)) %>%
  separate(Collection.Number, c("treatment", "Block.position"), sep = 1) %>% 
  filter(treatment != "w", Block.position != " 7.9(2)") # remove leaf with unknown treatment, and U7.9(2) which is another duplicate.
wind.WetLeafWts.2012$Block.position[which(wind.WetLeafWts.2012$Block.position == "7.9(1)")] <- "7.9" # remove (1)
wind.WetLeafWts.2012 <- wind.WetLeafWts.2012 %>%
  separate(Block.position, c("Block", "position")) %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         plant.id = paste0(Block, treatment, position)) %>%
  select(plant.id, leaf_wet.wt.g = Wet.Leaf.Weight, 
         leaf_percent.browned = Percent.Browned)

wind.trichomes.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Trichome Density - Wind Experiment 2012.csv", skip = 1) %>%
  tbl_df() %>% 
  .[-c(9,89), ] %>% # remove NA6.10 (unknown treatment) and U7.9-2 (duplicate sample)
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  separate(Collection.No., c("treatment","Block.position"), sep = 1) %>%
  separate(Block.position, c("Block", "position")) %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         plant.id = paste0(Block, treatment, position)) %>%
  select(plant.id, leaf_trichome.density = Trichome.Density) 

wind.DryLeafWt.2012 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Leaf Weights - Dry - Wind Experiment 2012.csv", skip = 1) %>%
  tbl_df() %>%
  .[-c(38,112), ] %>% # remove U7.9-2 (duplicate sample) and NA6.10 (unknown treatment)
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  separate(Collection.No., c("treatment","Block.position"), sep = 1) %>%
  separate(Block.position, c("Block", "position")) %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         plant.id = paste0(Block, treatment, position)) %>%
  select(plant.id, leaf_dry.wt.g = Dry.Leaf.Wt)

# combine wind 2012 trait data
wind.traits.2012 <- left_join(wind_plant_info, wind.arch2012) %>%
  left_join(., wind.trichomes.2012) %>%
  left_join(., wind.WetLeafWts.2012) %>%
  left_join(., wind.DryLeafWt.2012) %>%
  mutate(leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2012) %>%
  select(Year, block:leaf_trichome.density, leaf_WC, leaf_percent.browned)
glimpse(wind.traits.2012)

## wind 2013 - plant traits

# plant architecture 
wind.arch2013 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_growth_data_July10_2013.csv", skip = 1) %>%
  tbl_df() %>%
  mutate(treatment = ifelse(Wind.Exposure == "Exposed", "E", "U"),
         plant.id = paste0(Block,treatment,Plant.Position))
wind.arch2013$all.shoot.total.length <- rowSums(select(wind.arch2013, 
                                                       X1_shoot_length:X21_shoot_length))
wind.arch2013$all.shoot.count <- rowSums(select(wind.arch2013, 
                                                X1_shoot_length:X21_shoot_length) > 0)
wind.arch2013$all.shoot.avg.length <- rowMeans(select(wind.arch2013, 
                                                      X1_shoot_length:X21_shoot_length))
wind.arch2013 <- wind.arch2013 %>%
  select(plant.id, Dead = Dead_for_plant_survey, Height, 
         all.shoot.total.length, all.shoot.avg.length, all.shoot.count)

# leaf dry weight
wind.drywt.2013 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/SLA_leaf_wts_2013.csv") %>%
  .[-which(.$plant.code %in% names(which(table(.$plant.code) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.code.tmp = gsub("u","U",plant.code)) %>%
  mutate(plant.code.tmp2 = gsub("e","E",plant.code.tmp)) %>%
  select(plant.id = plant.code.tmp2, leaf_dry.wt.g = dry.wt.g)
wind.drywt.2013[which(wind.drywt.2013$plant.id == "3E6"), "leaf_dry.wt.g"] <- 0.0446 # type double checked with data.

# leaf wet weight
wind.wetwt.2013 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wet_leaf_weights_wind_experiment_2013.csv") %>%
  tbl_df() %>%
  mutate(Treat.tmp = ifelse(treatment == "u", "U", "E"),
         plant.id = paste0(block,Treat.tmp,position)) %>%
  select(plant.id, leaf_wet.wt.g = collection.weight..after.label.) %>%
  .[-which(.$plant.id %in% names(which(table(.$plant.id) > 1))), ] # removes all duplicates

# leaf area
wind.leafarea.2013 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Wind_Leaf_Scans_July_2013.csv") %>%
  .[-which(.$plant.code %in% names(which(table(.$plant.code) > 1))), ] %>% # removes all duplicates
  .[-103, ] %>% # innaccurate leaf measurement so I removed it.
  tbl_df() %>%
  mutate(plant.id = gsub(" ","",plant.code)) %>%
  select(plant.id, leaf_area.mm2 = leaf.area.mm2) 

# manage C:N data
batch1 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_Wind_run_1_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch1) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch1 <- batch1[-1, ] # remove first line which has no data
batch1$sample.id[which(batch1$sample.amount.mg == 5.15)] <- "3E01" # typo that I double checked with the original data.

batch2 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_Wind_run_2_results_CN.csv", stringsAsFactors = FALSE) # unknown why there is a warning error, loaded data looks okay
colnames(batch2) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch2 <- batch2[-1, ]
batch2$sample.id[2:3] <- c("6E1","2E5") # changed them to their appropriate ID

batch3 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_Wind_run_3_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch3) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent","possible.switch","notes")
batch3 <- batch3[-c(1,15,16), ] # also removed the ones that may have been possibly switched

batch4 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_wind_run_4_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch4) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch4$sample.id[c(5,8,11,14,17,19,20,23,25,27,30,31,32)] <- c("1E4","1E9","7E1","2E2","5E1","4E8","3E3","5E5","1E10","3E8","7E3","2E9","6E6") # changed them to their appropriate ID

batch7 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_wind_run_7_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch7) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch7 <- batch7[-c(1:42), ] # removed samples that were not part of the wind experiment

combined <- rbind.data.frame(batch1, batch2, batch3[ ,1:10], batch4, batch7) # didn't include notes from batch3
combined[ ,5:10] <- sapply(combined[ ,5:10], as.numeric) # convert necessary variables to numeric ones
combined$C_N <- combined$C.weight.percent/combined$N.weight.percent

duplicate.samples <- which(combined$sample.id %in% c("1E10")) # double checked original data and confirmed that I won't be able to resolve this duplicate, so I'm removing the entire sample since I don't know which one is erroneous.
row.remove <- which(combined$sample.id %in% c("blind standard","blind std", "Blind Standard", "standard 1", "standard 2", "standard 3", "bypass", "std2", "std 3", "std1")) # identify data from standards and bypasses

C_N_combined <- combined[-c(row.remove, duplicate.samples), ] %>% # removes unecessary data from standards and bypasses as as well as the duplicated samples
  tbl_df() %>%
  select(plant.id = sample.id, leaf_C_N = C_N, 
         leaf_C = C.weight.percent, leaf_N = N.weight.percent)


# spodoptera experiment #1 
spod.1 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/spodoptera/spodoptera_herbivory_wind_experiment_1.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp1 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.wt.exp1 = larva.wet.weight, 
         larva.survival.exp1 = survival.exp1)

# spodoptera experiment #2 
spod.2 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/spodoptera/spodoptera_herbivory_wind_experiment_2.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp2 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.wt.exp2 = larvae.wet.weight, 
         larva.survival.exp2 = survival.exp2)

## join together datasets for wind 2013 plant traits
wind.traits.2013 <- left_join(wind_plant_info, wind.arch2013) %>%
  left_join(., wind.wetwt.2013) %>%
  left_join(., wind.drywt.2013) %>%
  left_join(., wind.leafarea.2013) %>%
  left_join(., C_N_combined) %>%
  left_join(., spod.1) %>%
  left_join(., spod.2) %>%
  mutate(SLA = leaf_dry.wt.g/leaf_area.mm2,
         leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2013) %>%
  select(Year, block:all.shoot.count, SLA, leaf_WC, leaf_C_N:larva.survival.exp2)
glimpse(wind.traits.2013)

## bind together wind 2012 and 2013 plant trait data
wind.traits.df <- bind_rows(wind.traits.2012, wind.traits.2013)
glimpse(wind.traits.df)
write.csv(wind.traits.df, "~/Documents/Lanphere_Experiments/final_data/wind_trait_df.csv")

## manage ant-aphid data

# ant-aphid setup (plant info)
aa.setup <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Ant-Aphid_Experiment_Setup.csv") %>%
  tbl_df() %>% 
  select(Block, Aphid.Treatment = Aphids.or.No.Aphids, 
         Ant.Mound.Dist = Distant.to.Ant.Mound,
         Genotype, Plant_Position = Plant.Position) %>%
  mutate(plant.id = as.character(paste(Block, Plant_Position, sep=".")),
         Block = as.factor(Block))

# ant-aphid plant architecture data 2012
aa.arch2012 <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/ant_aphid_architecture_data_2012.csv") %>% 
  tbl_df() %>%
  mutate(plant.id = as.character(paste(Block, Plant.Position, sep = ".")),
         Dead = ifelse(Survived == 1, 0, 1)) %>%
  select(plant.id, Dead, Height = Max.Height, leaf_count = Mature.Green.Leaves, 
         immature.shoot.count = Shoots.Sprouting, mature.shoot.1:mature.shoot.18)

aa.arch2012$mature.shoot.total.length = rowSums(select(aa.arch2012, 
                                                       mature.shoot.1:mature.shoot.18))
aa.arch2012$mature.shoot.count = rowSums(select(aa.arch2012,
                                                mature.shoot.1:mature.shoot.18) > 0)
aa.arch2012$mature.shoot.avg.length = rowMeans(select(aa.arch2012, 
                                                       mature.shoot.1:mature.shoot.18))

aa.arch2012 <- aa.arch2012 %>%
  mutate(all.shoot.count = mature.shoot.count + immature.shoot.count) %>%
  select(plant.id, Dead, leaf_count, Height, mature.shoot.count, immature.shoot.count,
         all.shoot.count, mature.shoot.total.length, mature.shoot.avg.length)

# ant-aphid plant architecture data 2013
aa.arch2013 <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/ant_aphid_plant_growth_2013.csv", skip = 1) %>%
  tbl_df() %>%
  mutate(plant.id = as.character(paste(Block, Plant.Position, sep = ".")))  %>%
  select(plant.id, Dead:X22.shoot.length)

aa.arch2013$all.shoot.total.length = rowSums(select(aa.arch2013, X1.shoot.length:X22.shoot.length),
                                             na.rm = TRUE)
aa.arch2013$all.shoot.avg.length = rowMeans(select(aa.arch2013, X1.shoot.length:X22.shoot.length),
                                             na.rm = TRUE)
aa.arch2013$all.shoot.count = rowSums(select(aa.arch2013, X1.shoot.length:X22.shoot.length) > 0,
                                             na.rm = TRUE)

aa.arch2013 <- aa.arch2013 %>%
  select(plant.id, Dead, Height, all.shoot.total.length, all.shoot.avg.length, all.shoot.count)

# ant-aphid leaf quality data 2012
aa.trich2012 <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Trichome Density - Ant-Aphid Experiment 2012.csv") %>%
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.No.)) %>%
  select(plant.id, leaf_trichome.density = Trichome.Density)

aa.DryLeafWt2012 <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Leaf Weights - Dry - Ant-Aphid Experiment 2012.csv") %>%
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.No.)) %>%
  select(plant.id, leaf_dry.wt.g = Dry.Leaf.Weight)

aa.WetLeafWt2012 <- read.csv("~/Documents/Lanphere_Experiments/ant_aphid_experiment/data/Plant Traits/Ant-aphid Plant Trait Data.csv") %>% # checked for duplicates and there were none
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.Number)) %>%
  select(plant.id, leaf_wet.wt.g = Wet.Leaf.Weight,
         leaf_percent.browned = Percent.Browned) # note that I can only use percent browning as a covariate for the ant-aphid data, because I didn't sample plants if all of their leaves were substantially browned. 
aa.WetLeafWt2012[3,"leaf_wet.wt.g"] <- 0.103 # typo, double-checked this with real data.

# join together ant-aphid trait data 2012
aa.traits2012 <- left_join(aa.setup, aa.arch2012) %>%
  left_join(., aa.WetLeafWt2012) %>%
  left_join(., aa.DryLeafWt2012) %>%
  left_join(., aa.trich2012) %>%
  mutate(leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2012) %>%
  select(Year, Block:mature.shoot.avg.length, leaf_WC, leaf_percent.browned, leaf_trichome.density)

# join together ant-aphid trait data 2013
aa.traits2013 <- left_join(aa.setup, aa.arch2013) %>%
  mutate(Year = 2013) %>%
  select(Year, Block:all.shoot.count)

# bind together ant-aphid trait data from 2012 and 2013
aa.traits.df <- bind_rows(aa.traits2012, aa.traits2013)
glimpse(aa.traits.df)

write.csv(aa.traits.df, "~/Documents/Lanphere_Experiments/final_data/ant_aphid_trait_df.csv")
