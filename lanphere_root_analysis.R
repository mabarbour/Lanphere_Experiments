# lanphere root community analyses (fungus and bacteria)

# load required libraries
library(dplyr)
library(tidyr)

## upload and organize community data ----

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
  mutate(Genotype = factor(Genotype)) %>%
  filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) # need to figure out where "M" came from...

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
