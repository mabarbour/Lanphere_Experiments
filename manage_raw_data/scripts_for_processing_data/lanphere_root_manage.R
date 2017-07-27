## Description: manages raw data for root-associated communities

## load required libraries ----
#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
#library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(dplyr)
library(tidyr)
library(vegan)
#library(lme4)
#library(ggplot2)
#library(cowplot) # throwing an error for unknown reason
#library(effects)

## Upload and manage fungal community data ----

fungal <- read.csv('~/Documents/Lanphere_Experiments/Output/hfNormalizedCounts.csv', check.names = FALSE) 
f.OTUs <- fungal[ ,1]
t.fungal <- as.data.frame(t(fungal[ ,-1]))
colnames(t.fungal) <- as.character(f.OTUs)#fungal$X)
t.fungal$plant_ID <- rownames(t.fungal)

# need to split up plant ID into block, treatment, genotype and position
t.fungal$TreatGeno <- gsub("[[:digit:]]", "", t.fungal$plant_ID)
t.fungal$Block__Position <- gsub("[^[:digit:]]", "_", t.fungal$plant_ID)

# split up columns again into final data frame
fungal.df <- t.fungal %>%
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) %>% 
  separate(Block__Position, into = c("Block","Plant.Position")) %>%
  mutate(#Genotype = as.factor(Genotype),
         #Wind.Exposure = as.factor(Wind.Exposure),
         #Block = as.factor(Block),
         Plot_code = interaction(Block, Wind.Exposure))#,
         #Sample = interaction(Block, Wind.Exposure, Genotype, Plant.Position),
         #GxE = C(interaction(Genotype, Wind.Exposure), contr = "contr.sum", how.many = 9))

#comm <- colnames(select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position, -Plot_code, -Sample, -GxE)) # OTU names

#sort(colSums(fungal.df[ ,comm])/sum(colSums(fungal.df[ ,comm]))) # each OTU only makes up at most 1% of data, OTU_1 and OTU_14 are the most abundant

fungal.df <- data.frame(Block = fungal.df$Block, 
                        Wind.Exposure = ifelse(fungal.df$Wind.Exposure == "E", "Exposed","Unexposed"), 
                        Plot_code = fungal.df$Plot_code, 
                        Genotype = fungal.df$Genotype,
                        plant_ID = fungal.df$plant_ID,
                        #fungal.abund = rowSums(fungal.df[ ,f.OTUs]),
                        #fungal.rich = rowSums(fungal.df[ ,f.OTUs] > 0), 
                        fungal.rarerich = rarefy(round(fungal.df[ ,f.OTUs],0), min(rowSums(round(fungal.df[ ,f.OTUs],0)))),
                        #fungal.rarerich = rarefy(round(fungal.df[ ,f.OTUs],0), 2) - 1, # round abundance since rarefy only works on integers
                        fungal.df[ ,f.OTUs]) 

write.csv(x = fungal.df, file = "fungal.df.csv")

## Upload and manage bacterial community data ----

bacteria <- read.csv('~/Documents/Lanphere_Experiments/Output/hbNormalizedCounts.csv', check.names = FALSE) 
#bacteria <- read.csv("~/Documents/Lanphere_Experiments/bacteria_phyla_otutable.csv", check.names = FALSE, stringsAsFactors = FALSE)
b.OTUs <- bacteria[ ,1]
t.bacteria <- as.data.frame(t(bacteria[ ,-1]))
colnames(t.bacteria) <- as.character(b.OTUs)#bacteria$X)
t.bacteria$plant_ID <- rownames(t.bacteria)

# need to split up plant ID into block, treatment, genotype and position
t.bacteria$TreatGeno <- gsub("[[:digit:]]", "", t.bacteria$plant_ID)
t.bacteria$Block__Position <- gsub("[^[:digit:]]", "_", t.bacteria$plant_ID)

# split up columns again into final data frame
bacteria.df <- t.bacteria[-53, ] %>% # unknown Genotype and plant position
  #separate(TreatGeno, into = c("X","Wind.Exposure","Genotype"), sep = c(1,2)) %>%
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  #filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) %>% 
  filter(Genotype %in% c("Fb","Gb","Ib","Jb","Lb","Sb","Tb","Ub","Wb","Xb")) %>% 
  #separate(Block__Position, into = c("blank","Block","Plant.Position")) %>%
  separate(Block__Position, into = c("Block","Plant.Position","blank")) %>%
  mutate(#Genotype = as.factor(Genotype),
         #Wind.Exposure = as.factor(Wind.Exposure),
         #Wind.Expos.alt = ifelse(Wind.Exposure == "E","Exposed","Unexposed"),
         #Block = as.factor(Block),
         Plot_code = interaction(Block, Wind.Exposure))

bacteria.df$Genotype <- as.factor(unlist(strsplit(as.character(bacteria.df$Genotype), split='b', fixed=TRUE)))

bacteria.df$plant_ID <- with(bacteria.df, paste(Block, Wind.Exposure, Genotype, Plant.Position, sep = ""))

bacteria.df <- data.frame(Block = bacteria.df$Block, 
                          Wind.Exposure = ifelse(bacteria.df$Wind.Exposure == "E", "Exposed","Unexposed"), 
                          Plot_code = bacteria.df$Plot_code,
                          Genotype = bacteria.df$Genotype,
                          plant_ID = bacteria.df$plant_ID, 
                          #bacteria.abund = rowSums(bacteria.df[ ,b.OTUs]),
                          #bacteria.rich = rowSums(bacteria.df[ ,b.OTUs] > 0),
                          bacteria.rarerich = rarefy(round(bacteria.df[ ,b.OTUs],0), min(rowSums(round(bacteria.df[ ,b.OTUs])))), # rarefying down to sample with lowest abundance
                          #bacteria.rarerich = rarefy(round(bacteria.df[ ,b.OTUs],0), 2) - 1,
                          bacteria.df[ ,b.OTUs]) 

write.csv(x = bacteria.df, file = "bacteria.df.csv")


#bacteria.df$plant_code <- with(bacteria.df, interaction(Wind.Expos.alt, Block, Plant.Position))

#bacteria.df$GxE <- with(bacteria.df, C(interaction(Genotype, Wind.Exposure), contr = "contr.sum", how.many = 9))

#bact.comm <- colnames(select(bacteria.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position, -Plot_code, -Wind.Expos.alt, -blank, -GxE)) # OTU names

#sort(colSums(bacteria.df[ ,bact.comm])/sum(colSums(bacteria.df[ ,bact.comm]))) # each OTU only makes up less than 1% of data.



