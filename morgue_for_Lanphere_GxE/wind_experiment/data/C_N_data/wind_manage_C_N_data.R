# load libraries
library(dplyr)
library(lme4)
library(car)
#library(reshape2)
#library(ggplot2)
#library(GGally)

# source in useful functions
#source("~/Documents/miscellaneous_R/summarySE_function.R")

# set working directory
#setwd("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/")

######## Upload datasets and clean them up. 
# need to evaluate whether there are duplicates in the dataset or not. especially because excel may have messed up some of the sample IDs

batch1 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/Lanphere_Wind_run_1_results_CN.csv", stringsAsFactors = FALSE)
colnames(batch1) <- c("sample.time","unknown.variable", "sample.id", "sample.amount.mg", "N.reten.time.min","N.weight.mg","N.weight.percent","C.reten.time.min","C.weight.mg","C.weight.percent")
batch1 <- batch1[-1, ] # remove first line which has no data

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

# combine datasets together
combined <- rbind.data.frame(batch1, batch2, batch3[ ,1:10], batch4, batch7) # didn't include notes from batch3
combined[ ,5:10] <- sapply(combined[ ,5:10], as.numeric) # convert necessary variables to numeric ones
combined$C_N <- combined$C.weight.percent/combined$N.weight.percent

row.remove <- which(combined$sample.id %in% c("blind standard","blind std", "Blind Standard", "standard 1", "standard 2", "standard 3", "bypass", "std2", "std 3", "std1")) # identify data from standards and bypasses

#combined <- combined[-row.remove, ]

which(combined$sample.id %in% c("3E10","1E10")) # duplicated samples.

#combined <- combined[-c(7,36), ] # removing 2 of the duplicated samples. Both of these appeared to be outlying C_N values, but I should 

#combined[which(combined$sample.id %in% c("3E10","1E10")), ]

combined <- combined[-row.remove, ] %>% # removes unecessary data from standards and bypasses
  tbl_df() %>%
  select(plant.id = sample.id, C_N)

wind_plant_info <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_info.csv") %>%
  tbl_df() 

## Join datasets. 
wind_CN_data <- left_join(wind_plant_info, combined, by = "plant.id")

# identify duplicates
with(wind_CN_data, table(genotype, block, treatment)) # duplicates for Genotype G in Exposed treatments in Block 1 and 3
duplicated(wind_CN_data$plant.id) # there are a few duplicates
which(duplicated(wind_CN_data$plant.id))
wind_CN_data$plant.id[which(duplicated(wind_CN_data$plant.id))]
which(wind_CN_data$plant.id %in% c("3E10","1E10")) # duplicated samples.


## exploratory analyses
plot(log(C_N) ~ genotype, wind_CN_data) # clearly a strong genotype effect

# strong genotype effect, but no GxE or E effect. Doesn't matter whether I log transform it or not and doesn't matter which duplicate samples I remove.
CN.lmer <- lmer(log(C_N) ~ treatment + (treatment|genotype) + (1|block), data = wind_CN_data[-c(51,52), ])
summary(CN.lmer)
Anova(CN.lmer) # very strong genotype effect, no GxE or environmental effect.
plot(CN.lmer)

