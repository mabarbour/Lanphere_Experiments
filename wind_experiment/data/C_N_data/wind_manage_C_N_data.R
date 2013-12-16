# load libraries
library(reshape2)
library(ggplot2)
library(GGally)

# source in useful functions
source("~/Documents/miscellaneous_R/summarySE_function.R")

# set working directory
setwd("~/Documents/Lanphere_Experiments/wind_experiment/data/C_N_data/")

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

combined <- combined[-row.remove, ] # removes unecessary data from standards and bypasses

wind_plant_info <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/wind_plant_info.csv")

##### Merge datasets. NEED TO INVESTIGATE DUPLICATES
wind_CN_data <- merge(wind_plant_info, combined, by.x = "plant.id", by.y = "sample.id", all.y = TRUE) # retain all of C_N data
dim(wind_CN_data)
duplicated(wind_CN_data$plant.id) # there are a few duplicates
which(duplicated(wind_CN_data$plant.id))
wind_CN_data$plant.id[which(duplicated(wind_CN_data$plant.id))]

#### Explore data
hist(wind_CN_data$C_N) # pretty broad distribution
hist(log(wind_CN_data$C_N)) # log transformation normalizes data distribution more
ggpairs(wind_CN_data[ ,c("C.weight.percent", "N.weight.percent", "C_N")])

plot(C_N ~ genotype, wind_CN_data)
plot(C_N ~ treatment, wind_CN_data)
plot(C_N ~ as.factor(block), wind_CN_data)

with(wind_CN_data, table(genotype, treatment))

summarizeCN <- summarySE(data = wind_CN_data, measurevar = "C_N", groupvars = c("treatment","genotype")) # note that sample sizes are quite low for some of the genotypes
ggplot(data = summarizeCN, aes(x = treatment, y = C_N, group = genotype, color = genotype)) + geom_line()

# preliminary analyses for split plot design. need to triple check that this is the correct coding or whether a mixed effects model would be better
aov_C_N <- aov(log(C_N) ~ genotype * treatment + Error(as.factor(block)/treatment), data= wind_CN_data, qr=TRUE)
summary(aov_C_N) # does this not allow me examine the effect of treatment? Right now, genotype appears to have a strong effect. 
