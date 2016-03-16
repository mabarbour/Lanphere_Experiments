# lanphere root community analyses (fungus and bacteria)

# load required libraries
library(dplyr)
library(tidyr)
library(vegan)

## Fungal community data ----

# fungal
fungal <- read.csv('~/Documents/Lanphere_Experiments/Output/hfNormalizedCounts.csv', check.names = FALSE) 
OTUs <- fungal[ ,1]
t.fungal <- as.data.frame(t(fungal[ ,-1]))
colnames(t.fungal) <- as.character(OTUs)#fungal$X)
t.fungal$plant_ID <- rownames(t.fungal)

# need to split up plant ID into block, treatment, genotype and position
t.fungal$TreatGeno <- gsub("[[:digit:]]", "", t.fungal$plant_ID)
t.fungal$Block__Position <- gsub("[^[:digit:]]", "_", t.fungal$plant_ID)

# split up columns again into final data frame
fungal.df <- t.fungal %>%
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  separate(Block__Position, into = c("Block","Plant.Position")) %>%
  mutate(Genotype = as.factor(Genotype),
         Wind.Exposure = as.factor(Wind.Exposure),
         Block = as.factor(Block),
         Plot_code = interaction(Block, Wind.Exposure)) %>%
  filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) # need to figure out where "M" came from...

wind.effect <- fungal.df %>%
  select(-Genotype, -Plant.Position, -plant_ID) %>%
  group_by(Block, Wind.Exposure, Plot_code) %>%
  summarise_each(funs(mean))

## analysis
comm <- colnames(select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position, -Plot_code)) # OTU names

sort(colSums(fungal.df[ ,comm])/sum(colSums(fungal.df[ ,comm]))) # each OTU only makes up at most 1% of data, OTU_1 and OTU_14 are the most abundant

# Interaction effect
adonis(fungal.df[ ,comm] ~ Wind.Exposure*Genotype, 
       data = fungal.df, 
       strata = fungal.df$Block, 
       method = "horn")

# Genotype effect
geno.adonis <- adonis(fungal.df[ ,comm] ~ Genotype, 
       data = fungal.df, 
       strata = fungal.df$Plot_code, 
       method = "horn")
geno.adonis
geno.adonis$coefficients[ ,1:10]

fungal.permdisp <- betadisper(vegdist(fungal.df[ ,comm],
                                      method = "horn"), 
                              group = fungal.df$Genotype, 
                              bias.adjust = TRUE)
anova(fungal.permdisp) # okay
plot(fungal.permdisp) # doesn't look like too much differentiation among communities...
boxplot(fungal.permdisp)

meandist(vegdist(fungal.df[ ,comm], 
                 method = "horn"), 
         group = fungal.df$Genotype)

# Wind exposure effect
adonis(wind.effect[ ,comm] ~ Wind.Exposure, 
       data = wind.effect, 
       strata = wind.effect$Block, 
       method = "horn")

## Bacterial community data ----
bacteria <- read.csv('~/Documents/Lanphere_Experiments/Output/hbNormalizedCounts.csv', check.names = FALSE) 
OTUs <- bacteria[ ,1]
t.bacteria <- as.data.frame(t(bacteria[ ,-1]))
colnames(t.bacteria) <- as.character(OTUs)#bacteria$X)
t.bacteria$plant_ID <- rownames(t.bacteria)

# need to split up plant ID into block, treatment, genotype and position
t.bacteria$TreatGeno <- gsub("[[:digit:]]", "", t.bacteria$plant_ID)
t.bacteria$Block__Position <- gsub("[^[:digit:]]", "_", t.bacteria$plant_ID)

## need to remove "b" at end of plant_ID and codes

# split up columns again into final data frame
bacteria.df <- t.bacteria[-53, ] %>% # unknown Genotype and plant position
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  separate(Block__Position, into = c("Block","Plant.Position","blank")) %>%
  mutate(Genotype = as.factor(Genotype),
         Wind.Exposure = as.factor(Wind.Exposure),
         Block = as.factor(Block),
         Plot_code = interaction(Block, Wind.Exposure)) %>%
  filter(Genotype %in% c("Fb","Gb","Ib","Jb","Lb","Sb","Tb","Ub","Wb","Xb")) # need to figure out where "M" came from...

wind.effect.bact <- bacteria.df %>%
  select(-Genotype, -Plant.Position, -plant_ID, -blank) %>%
  group_by(Block, Wind.Exposure, Plot_code) %>%
  summarise_each(funs(mean))

## analysis
bact.comm <- colnames(select(bacteria.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position, -Plot_code, -blank)) # OTU names

sort(colSums(bacteria.df[ ,bact.comm])/sum(colSums(bacteria.df[ ,bact.comm]))) # each OTU only makes up less than 1% of data.

# Interaction effect
adonis(bacteria.df[ ,bact.comm] ~ Wind.Exposure*Genotype, 
       data = bacteria.df, 
       strata = bacteria.df$Block, 
       method = "horn")

# Genotype effect
adonis(bacteria.df[ ,bact.comm] ~ Genotype, 
                      data = bacteria.df, 
                      strata = bacteria.df$Plot_code, 
                      method = "horn")

bacteria.permdisp <- betadisper(vegdist(bacteria.df[ ,bact.comm],
                                      method = "horn"), 
                              group = bacteria.df$Genotype, 
                              bias.adjust = TRUE)
anova(bacteria.permdisp) # okay
plot(bacteria.permdisp) # doesn't look like too much differentiation among communities...
boxplot(bacteria.permdisp)

meandist(vegdist(bacteria.df[ ,bact.comm], 
                 method = "horn"), 
         group = bacteria.df$Genotype)

# Wind exposure effect. No effect once I aggregate.
adonis(wind.effect.bact[ ,bact.comm] ~ Wind.Exposure, 
       data = wind.effect.bact, 
       strata = wind.effect.bact$Block, 
       method = "horn")

## OLD ----
# lmer
library(pbkrtest)
library(lmerTest)

otu14 <- lmer(OTU_1606 ~ Wind.Exposure*Genotype +
               (1|Block) + (1|Block:Wind.Exposure),
             data = fungal.df)
plot(otu14)
summary(otu14)
anova(otu14, ddf = "Kenward-Roger")


## load library for analysis
library(mvabund)
library(vegan)

fungal.mv <- mvabund(select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position))
sort(colSums(fungal.mv))
OTU.common <- which(colSums(fungal.mv) > 800)
fungal.com.mv <- fungal.mv[ ,OTU.common]

meanvar.plot(fungal.mv ~ fungal.df$Wind.Exposure)

plot.mvabund(fungal.com.mv, type = "bx")
plot.mvformula(fungal.com.mv ~ fungal.df$Block, legend = FALSE)

fungal.many <- manyany("lm", fungal.com.mv, 
                       data = fungal.df, 
                       log(fungal.com.mv+1) ~ Wind.Exposure*Genotype + Block, 
                       composition = TRUE)
fungal.null <- manyany("glm", fungal.com.mv, data = fungal.df, fungal.com.mv ~ 1, composition = TRUE, family = poisson(link = "log"))
anova.manyany(fungal.many, fungal.null, test = "LR")
summary(fungal.many, p.uni = "unadjusted")
plot(fungal.many)

fungal.community <- select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position)
OTU.common <- which(colSums(fungal.community) > 800)
fungal.adonis <- adonis(fungal.community[ ,OTU.common] ~ Wind.Exposure*Genotype + Block, fungal.df, 
                        #strata = fungal.df$Block, 
                        method = "horn")
fungal.adonis

fungal.permdisp <- betadisper(vegdist(fungal.community[ ,OTU.common], method = "horn"), group = fungal.df$Wind.Exposure, bias.adjust = TRUE)
anova(fungal.permdisp)
plot(fungal.permdisp)
boxplot(fungal.permdisp)

meandist(vegdist(fungal.community[ ,OTU.common], method = "horn"), group = fungal.df$Block)
