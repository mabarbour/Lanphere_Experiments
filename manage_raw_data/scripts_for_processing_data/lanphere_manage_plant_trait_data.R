########################################
## Description: This script organizes plant architecture and leaf quality data collected in the summers of 2012 and 2013 as part of the Lanphere dune Experiments.
## Code author: Matt Barbour
## Email: barbour@zoology.ubc.ca
########################################

#### load required packages ----
library(dplyr)
library(tidyr)
library(missMDA) # for imputing missing values in PCA

#### Manage Wind data ----

## upload wind plant info ----
wind_plant_info <- read.csv("manage_raw_data/raw_data/wind_plant_info.csv") %>%
  tbl_df() %>%
  mutate(block = as.factor(block)) %>%
  select(block, treatment, genotype, plant.position, plant.id)

## plant architecture 2012 ----
wind.arch2012 <- read.csv("manage_raw_data/raw_data/wind_plant_traits_field_summer_2012.csv") %>%
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

## leaf wet weights 2012 ----
wind.WetLeafWts.2012 <- read.csv("manage_raw_data/raw_data/Wind_leaf_trait_data.csv") %>%
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

## leaf trichome densities 2012 ----
wind.trichomes.2012 <- read.csv("manage_raw_data/raw_data/Trichome_Density_Wind_Experiment_2012.csv", skip = 1) %>%
  tbl_df() %>% 
  .[-c(9,89), ] %>% # remove NA6.10 (unknown treatment) and U7.9-2 (duplicate sample)
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  separate(Collection.No., c("treatment","Block.position"), sep = 1) %>%
  separate(Block.position, c("Block", "position")) %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         plant.id = paste0(Block, treatment, position)) %>%
  select(plant.id, leaf_trichome.density = Trichome.Density) 

## leaf dry weights 2012 ----
wind.DryLeafWt.2012 <- read.csv("manage_raw_data/raw_data/Leaf_Weights_Dry_Wind_Experiment_2012.csv", skip = 1) %>%
  tbl_df() %>%
  .[-c(38,112), ] %>% # remove U7.9-2 (duplicate sample) and NA6.10 (unknown treatment)
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  separate(Collection.No., c("treatment","Block.position"), sep = 1) %>%
  separate(Block.position, c("Block", "position")) %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         plant.id = paste0(Block, treatment, position)) %>%
  select(plant.id, leaf_dry.wt.g = Dry.Leaf.Wt)

## combine wind 2012 trait data ----
wind.traits.2012 <- left_join(wind_plant_info, wind.arch2012) %>%
  left_join(., wind.trichomes.2012) %>%
  left_join(., wind.WetLeafWts.2012) %>%
  left_join(., wind.DryLeafWt.2012) %>%
  mutate(leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2012) %>%
  select(Year, block:leaf_trichome.density, leaf_WC, leaf_percent.browned)
glimpse(wind.traits.2012)

wind.traits.2012$Height[which(wind.traits.2012$Height > 80)] <- NA # replacing biologically unreasonable values with NA
wind.traits.2012$leaf_WC[which(wind.traits.2012$leaf_WC < 0)] <- NA # replacing biologically unreasonable values with NA

## perform PCA with plant traits for 2012
trait.mat.2012 <- as.matrix(wind.traits.2012[ ,c("Height","all.shoot.count","all.shoot.avg.length","leaf_WC","leaf_trichome.density")])
(trait.nb.2012 <- estim_ncpPCA(trait.mat.2012, method.cv = "loo", scale = TRUE)) # 4 PCs
w.trait.imp.2012 <- imputePCA(trait.mat.2012, ncp = 4, scale = TRUE)

trait.PCA.12 <- prcomp(w.trait.imp.2012$completeObs, scale = TRUE)
biplot(trait.PCA.12)
summary(trait.PCA.12) # first 2 PCs explain 60% of the variance
trait.PCA.12$rotation <- trait.PCA.12$rotation*-1
trait.PCA.12$x <- trait.PCA.12$x*-1
biplot(trait.PCA.12) 
trait.PCA.12$rotation # high values of PC1 indicate larger plants, high values of PC2 indicate plants with greater leaf WC and trichome density.

write.csv(data.frame(trait.PCA.12$rotation), "final_data/wind.trait.PCA.2012.csv")

wind.traits.2012 <- mutate(wind.traits.2012, trait.PC1 = trait.PCA.12$x[ ,"PC1"], trait.PC2 = trait.PCA.12$x[ ,"PC2"])

## plant architecture 2013 ----
wind.arch2013 <- read.csv("manage_raw_data/raw_data/wind_plant_growth_data_July10_2013.csv", skip = 1) %>%
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

## leaf dry weight 2013 ----
wind.drywt.2013 <- read.csv("manage_raw_data/raw_data/SLA_leaf_wts_2013.csv") %>%
  .[-which(.$plant.code %in% names(which(table(.$plant.code) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.code.tmp = gsub("u","U",plant.code)) %>%
  mutate(plant.code.tmp2 = gsub("e","E",plant.code.tmp)) %>%
  select(plant.id = plant.code.tmp2, leaf_dry.wt.g = dry.wt.g)
wind.drywt.2013[which(wind.drywt.2013$plant.id == "3E6"), "leaf_dry.wt.g"] <- 0.0446 # type double checked with data.

## leaf wet weight 2013 ----
wind.wetwt.2013 <- read.csv("manage_raw_data/raw_data/wet_leaf_weights_wind_experiment_2013.csv") %>%
  tbl_df() %>%
  mutate(Treat.tmp = ifelse(treatment == "u", "U", "E"),
         plant.id = paste0(block,Treat.tmp,position)) %>%
  select(plant.id, leaf_wet.wt.g = collection.weight..after.label.) %>%
  .[-which(.$plant.id %in% names(which(table(.$plant.id) > 1))), ] # removes all duplicates

## leaf area 2013 ----
wind.leafarea.2013 <- read.csv("manage_raw_data/raw_data/Wind_Leaf_Scans_July_2013.csv") %>%
  .[-which(.$plant.code %in% names(which(table(.$plant.code) > 1))), ] %>% # removes all duplicates
  .[-103, ] %>% # innaccurate leaf measurement so I removed it.
  tbl_df() %>%
  mutate(plant.id = gsub(" ","",plant.code)) %>%
  select(plant.id, leaf_area.mm2 = leaf.area.mm2) 

## manage C:N data 2013 ----
batch1 <- read.csv("manage_raw_data/raw_data/Lanphere_Wind_run_1_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch1) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch1 <- batch1[-1, ] # remove first line which has no data
batch1$sample.id[which(batch1$sample.amount.mg == 5.15)] <- "3E01" # typo that I double checked with the original data.

batch2 <- read.csv("manage_raw_data/raw_data/Lanphere_Wind_run_2_results_CN.csv", stringsAsFactors = FALSE) # unknown why there is a warning error, loaded data looks okay
colnames(batch2) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch2 <- batch2[-1, ]
batch2$sample.id[2:3] <- c("6E1","2E5") # changed them to their appropriate ID

batch3 <- read.csv("manage_raw_data/raw_data/Lanphere_Wind_run_3_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch3) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent","possible.switch","notes")
batch3 <- batch3[-c(1,15,16), ] # also removed the ones that may have been possibly switched

batch4 <- read.csv("manage_raw_data/raw_data/Lanphere_wind_run_4_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch4) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch4$sample.id[c(5,8,11,14,17,19,20,23,25,27,30,31,32)] <- c("1E4","1E9","7E1","2E2","5E1","4E8","3E3","5E5","1E10","3E8","7E3","2E9","6E6") # changed them to their appropriate ID

batch7 <- read.csv("manage_raw_data/raw_data/Lanphere_wind_run_7_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch7) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch7 <- batch7[-c(1:42), ] # removed samples that were not part of the wind experiment

combined <- rbind.data.frame(batch1, batch2, batch3[ ,1:10], batch4, batch7) # didn't include notes from batch3
combined[ ,5:10] <- sapply(combined[ ,5:10], as.numeric) # convert necessary variables to numeric ones
combined$C_N <- with(combined, C.weight.mg/N.weight.mg) #combined$C.weight.percent/combined$N.weight.percent # incorrect, should be C.weight.mg/N.weight.mg

duplicate.samples <- which(combined$sample.id %in% c("1E10")) # double checked original data and confirmed that I won't be able to resolve this duplicate, so I'm removing the entire sample since I don't know which one is erroneous.
row.remove <- which(combined$sample.id %in% c("blind standard","blind std", "Blind Standard", "standard 1", "standard 2", "standard 3", "bypass", "std2", "std 3", "std1")) # identify data from standards and bypasses

C_N_combined <- combined[-c(row.remove, duplicate.samples), ] %>% # removes unecessary data from standards and bypasses as as well as the duplicated samples
  tbl_df() %>%
  select(plant.id = sample.id, leaf_C_N = C_N, 
         leaf_C.mg = C.weight.mg, leaf_N.mg = N.weight.mg,
         leaf_C.perc = C.weight.percent, leaf_N.perc = N.weight.percent)

## spodoptera experiment #1 2013 ----
spod.1 <- read.csv("manage_raw_data/raw_data/spodoptera_herbivory_wind_experiment_1.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp1 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.wt.exp1 = larva.wet.weight, 
         larva.survival.exp1 = survival.exp1)

## spodoptera experiment #2 2013 ----
spod.2 <- read.csv("manage_raw_data/raw_data/spodoptera_herbivory_wind_experiment_2.csv", skip = 1) %>%
  tbl_df() %>%
  filter(missing.larva != "missing_larva") %>%
  mutate(treatment = ifelse(treatment == "e", "E", "U"),
         survival.exp2 = ifelse(dead == "d", 0, 1),
         plant.id = paste0(block, treatment, position)) %>%
  select(plant.id, larva.wet.wt.exp2 = larvae.wet.weight, 
         larva.survival.exp2 = survival.exp2)

## root C:N 2013 ----
rootCN <- read.csv('manage_raw_data/raw_data/RootCNdataLanphereDunes.csv', stringsAsFactors = FALSE) 
rootCN[2,"Sample"] <- "1EF3" # double-checked with plant info and this is the correct code
rootCN[76,"Sample"] <- "10UI6"

rootCN.df <- rootCN %>%
  tbl_df() %>%
  # need to split up plant ID into block, treatment, genotype and position
  mutate(TreatGeno = gsub("[[:digit:]]", "", Sample),
         Block__Position = gsub("[^[:digit:]]", "_", Sample)) %>%
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  separate(Block__Position, into = c("Block","Plant.Position")) %>%
  mutate(plant.id = paste(Block, Wind.Exposure, Plant.Position, sep = ""),
         root_CN = Carbon_mg/Nitrogen_mg) %>%
  filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) %>% 
  filter(Nitrogen_percent > 0 & Carbon_percent > 0 & Carbon_percent < 100 & root_CN < 700) %>% # restricting to biologically reasonable values
  select(plant.id, root_N.perc = Nitrogen_percent, root_N.mg = Nitrogen_mg, root_C.perc = Carbon_percent, root_C.mg = Carbon_mg, root_CN)

## join together datasets for wind 2013 plant traits ----
wind.traits.2013 <- left_join(wind_plant_info, wind.arch2013) %>%
  left_join(., wind.wetwt.2013) %>%
  left_join(., wind.drywt.2013) %>%
  left_join(., wind.leafarea.2013) %>%
  left_join(., C_N_combined) %>%
  left_join(., spod.1) %>%
  left_join(., spod.2) %>%
  left_join(., rootCN.df) %>%
  mutate(SLA = leaf_area.mm2/leaf_dry.wt.g*1000, # mm2/mg
         leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2013) #%>%
  #select(Year, block:all.shoot.count, SLA, leaf_WC, leaf_C_N:larva.survival.exp2)
glimpse(wind.traits.2013)

wind.traits.2013$Height[which(wind.traits.2013$Height > 80)] <- NA # replacing biologically unreasonable values with NA
wind.traits.2013$leaf_WC[which(wind.traits.2013$leaf_WC < 0)] <- NA # replacing biologically unreasonable values with NA

## perform PCA with plant traits for 2013 ----
trait.mat.2013 <- as.matrix(wind.traits.2013[ ,c("Height","all.shoot.count","all.shoot.avg.length","SLA","leaf_WC","leaf_C_N")])
(trait.nb.2013 <- estim_ncpPCA(trait.mat.2013, method.cv = "loo", scale = TRUE)) # 5 PCs
w.trait.imp.2013 <- imputePCA(trait.mat.2013, ncp = 5, scale = TRUE)

trait.PCA.13 <- prcomp(w.trait.imp.2013$completeObs, scale = TRUE)
biplot(trait.PCA.13 )
summary(trait.PCA.13 ) # first 2 PCs explain 72% of the variance
#trait.PCA$rotation <- trait.PCA$rotation*-1
#trait.PCA$x <- trait.PCA$x*-1
#biplot(trait.PCA) 
trait.PCA.13$rotation # High values of PC1 indicate larger plants, but with lower leaf WC and SLA. High values of PC2 indicate higher values of leaf C:N, but lower values of SLA and and leaf WC

write.csv(data.frame(trait.PCA.13$rotation), "final_data/wind.trait.PCA.2013.csv")

wind.traits.2013 <- mutate(wind.traits.2013, trait.PC1 = trait.PCA.13$x[ ,"PC1"], trait.PC2 = trait.PCA.13$x[ ,"PC2"])

## bind together wind 2012 and 2013 plant trait data and save to a new file ----
wind.traits.df <- bind_rows(wind.traits.2012, wind.traits.2013)
glimpse(wind.traits.df)

# subset data for manuscript
wind.traits.df <- wind.traits.df %>% mutate(treat.tmp = ifelse(treatment == "Exposed","E","U"), plant_ID = paste(block, treat.tmp, genotype, plant.position, sep = "")) %>% select(Year, Block = block, Wind.Exposure = treatment, Genotype = genotype, plant_ID, Dead, Height, all.shoot.count, all.shoot.avg.length, leaf_trichome.density, leaf_WC, SLA, leaf_area.mm2, leaf_C_N:leaf_N.perc, trait.PC1, trait.PC2, root_N.perc:root_CN)

#wind.traits.df$Height[which(wind.traits.df$Height > 80)] <- NA # replacing biologically unreasonable values with NA
#wind.traits.df$leaf_WC[which(wind.traits.df$leaf_WC < 0)] <- NA # replacing biologically unreasonable values with NA

write.csv(wind.traits.df, "final_data/wind_trait_df.csv")

#### Manage ant-aphid data ----

## ant-aphid setup (plant info) ----
aa.setup <- read.csv("manage_raw_data/raw_data/Ant_Aphid_Experiment_Setup.csv") %>%
  tbl_df() %>% 
  select(Block, Aphid.Treatment = Aphids.or.No.Aphids, 
         Ant.Mound.Dist = Distant.to.Ant.Mound,
         Genotype, Plant_Position = Plant.Position) %>%
  mutate(plant.id = as.character(paste(Block, Plant_Position, sep=".")),
         Block = as.factor(Block))

## ant-aphid plant architecture data 2012 ----
aa.arch2012 <- read.csv("manage_raw_data/raw_data/ant_aphid_architecture_data_2012.csv") %>% 
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

## ant-aphid plant architecture data 2013 ----
aa.arch2013 <- read.csv("manage_raw_data/raw_data/ant_aphid_plant_growth_2013.csv", skip = 1) %>%
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

## ant-aphid leaf quality data 2012 ----
aa.trich2012 <- read.csv("manage_raw_data/raw_data/Trichome_Density_Ant-Aphid_Experiment_2012.csv") %>%
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.No.)) %>%
  select(plant.id, leaf_trichome.density = Trichome.Density)

aa.DryLeafWt2012 <- read.csv("manage_raw_data/raw_data/Leaf_Weights_Dry_Ant-Aphid_Experiment_2012.csv") %>%
  .[-which(.$Collection.No. %in% names(which(table(.$Collection.No.) > 1))), ] %>% # removes all duplicates
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.No.)) %>%
  select(plant.id, leaf_dry.wt.g = Dry.Leaf.Weight)

aa.WetLeafWt2012 <- read.csv("manage_raw_data/raw_data/Ant_aphid_Plant_Trait_Data.csv") %>% # checked for duplicates and there were none
  tbl_df() %>%
  mutate(plant.id = as.character(Collection.Number)) %>%
  select(plant.id, leaf_wet.wt.g = Wet.Leaf.Weight,
         leaf_percent.browned = Percent.Browned) # note that I can only use percent browning as a covariate for the ant-aphid data, because I didn't sample plants if all of their leaves were substantially browned. 
aa.WetLeafWt2012[3,"leaf_wet.wt.g"] <- 0.103 # typo, double-checked this with real data.

## join together ant-aphid trait data 2012 ----
aa.traits2012 <- left_join(aa.setup, aa.arch2012) %>%
  left_join(., aa.WetLeafWt2012) %>%
  left_join(., aa.DryLeafWt2012) %>%
  left_join(., aa.trich2012) %>%
  mutate(leaf_WC = (leaf_wet.wt.g - leaf_dry.wt.g)/leaf_dry.wt.g,
         Year = 2012) %>%
  select(Year, Block:mature.shoot.avg.length, leaf_WC, leaf_percent.browned, leaf_trichome.density)

## Get principal component scores for plant trait in 2012 ----
library(missMDA) # for imputing PCA for missing trait values

aa.trait.mat.2012 <- as.matrix(aa.traits2012[ ,c("Height","all.shoot.count","mature.shoot.avg.length","leaf_WC","leaf_trichome.density")])

(aa.trait.nb.2012 <- estim_ncpPCA(aa.trait.mat.2012, method.cv = "loo", scale = TRUE)) #  4 PCs
aa.trait.imp.2012 <- imputePCA(aa.trait.mat.2012, ncp = 4, scale = TRUE)

aa.trait.PCA <- prcomp(aa.trait.imp.2012$completeObs, scale = TRUE)
biplot(aa.trait.PCA)
summary(aa.trait.PCA) # first 2 components explain 68% of the variance

aa.trait.PCA$rotation <- aa.trait.PCA$rotation*-1
aa.trait.PCA$x <- aa.trait.PCA$x*-1
biplot(aa.trait.PCA) 
aa.trait.PCA$rotation # positive values of PC1 indicate larger plants (taller, more shoots, and longer shoots). Positive values of PC2 indicate plants with more shoots, but smaller in height, less leaf WC and few trichomes.

write.csv(data.frame(aa.trait.PCA$rotation), "final_data/aa.trait.PCA.2012.csv")

aa.traits2012 <- mutate(aa.traits2012, trait.PC1 = aa.trait.PCA$x[ ,"PC1"], trait.PC2 = aa.trait.PCA$x[ ,"PC2"])

## join together ant-aphid trait data 2013 ----
aa.traits2013 <- left_join(aa.setup, aa.arch2013) %>%
  mutate(Year = 2013) %>%
  select(Year, Block:all.shoot.count)

## bind together ant-aphid trait data from 2012 and 2013 and save to a file -----
aa.traits.df <- bind_rows(aa.traits2012, aa.traits2013)
glimpse(aa.traits.df)

aa.traits.df <- aa.traits.df %>% mutate(plant_ID = paste(Block, Ant.Mound.Dist, Aphid.Treatment, Genotype, Plant_Position, sep = "_")) %>% select(Year, Block, Aphid.treatment = Aphid.Treatment, Ant.mound.dist = Ant.Mound.Dist, Genotype, plant_ID, Height, all.shoot.count, mature.shoot.avg.length, leaf_WC, leaf_trichome.density, trait.PC1, trait.PC2)

write.csv(aa.traits.df, "final_data/ant_aphid_trait_df.csv")
