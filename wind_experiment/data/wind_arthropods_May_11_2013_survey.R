
###### upload required libraries
library(ggplot2)
library(BiodiversityR)
library(plyr)
library(lme4)

###### upload and manage the data the data
wind_arthropods_May_11_2013 <- read.csv("~/Documents/Lanphere_Experiments/wind_experiment/data/Wind_Arthropod_data_survey_1_May_10_2013.csv", skip=4)

str(wind_arthropods_May_11_2013)
wind_arthropods_May_11_2013$Block <- as.factor(wind_arthropods_May_11_2013$Block)

wind_arthropods_May_11_2013_subset <- subset(wind_arthropods_May_11_2013, select=-c(Calop_pre_tentmine, Leaf_bundle, Calop_LTF, Calop_cocoon, miscellaneous_arthropods_present_see_notes, Notes)) # remove Calop_pre_tentmine because of its difficulty in reliably sampling.  Also remove leaf_bundle temporarily because I don't know yet how I want to treat the "old" mines.  Also, remove Calop_LTF because I think the same larva may be responsible for multiple leaf folds on tiny, crappy looking plants.  Also remove Calop_cocoon for now although it may yield important information on survival.  Also remove miscellaneous arthropods and notes

wind_arthropods_May_11_2013_subset <- wind_arthropods_May_11_2013_subset[-which(wind_arthropods_May_11_2013_subset$Dead == 1), ] # remove all plants labelled as dead
wind_arthropods_May_11_2013_subset <- wind_arthropods_May_11_2013_subset[ ,-21]


###### Exploratory data analysis

# how many plants in each treatment
with(wind_arthropods_May_11_2013_subset, table(Genotype, Wind.Exposure)) # Genotypes G and S are creating an unbalanced data set and seem particularly susceptible to wind exposure.

# create arthropod abundance and richness variables
wind_arthropods_May_11_2013_subset$arthropod_abundance <- rowSums(wind_arthropods_May_11_2013_subset[ ,5:20])
wind_arthropods_May_11_2013_subset$arthropod_richness <- rowSums(wind_arthropods_May_11_2013_subset[ ,5:20] > 0)

## updated data analysis
arth.abund.lmer <- lmer(arthropod_abundance ~ Wind.Exposure*Genotype + (1|Block), data = wind_arthropods_May_11_2013_subset)
summary(arth.abund.lmer)
Anova(arth.abund.lmer) # clear effect of wind exposure on arthropod abundance



# ggplot graphs
abundance <- ggplot(wind_arthropods_May_11_2013_subset, aes(x=Wind.Exposure, y=arthropod_abundance, color=Genotype))
abundance + geom_boxplot() # unexposed plots appear to have higher abundances, a couple of outlying datapoints though.

richness <- ggplot(wind_arthropods_May_11_2013_subset, aes(x=Wind.Exposure, y=arthropod_richness, color=Genotype))
richness + geom_boxplot() # same pattern as abundance

qplot(wind_arthropods_May_11_2013_subset$arthropod_abundance) # highly right skewed distributions
qplot(wind_arthropods_May_11_2013_subset$arthropod_richness) # highly right skewed distributions

abundance_means <- ddply(wind_arthropods_May_11_2013_subset, .(Wind.Exposure,Genotype), summarise, val=mean(arthropod_abundance))
abundance_medians <- ddply(wind_arthropods_May_11_2013_subset, .(Wind.Exposure,Genotype), summarise, val=median(arthropod_abundance)) # Not enough variation to detect a strong pattern

abundance + geom_point(data=abundance_means, aes(y=val)) + geom_line(data=abundance_means, aes(y=val, group=Genotype)) # One of the outlying datapoints comes from the 8 Calop_tent_mines with the 4 leaf_bumps. The second outlying data point came from the same plot with also a very high Calop_tent_mine count (6).

richness_means <- ddply(wind_arthropods_May_11_2013_subset, .(Wind.Exposure,Genotype), summarise, val=mean(arthropod_richness))

richness + geom_point(data=richness_means, aes(y=val)) + geom_line(data=richness_means, aes(y=val, group=Genotype)) # suggests a possible interactive effect of plant genotype on arthropod richness. However, the mixed effects model suggest that wind exposure dominates plant genotype.

### Model the data

# problem with the anova model is that mortality has created an imbalance in my data.
aov_richness <- aov(arthropod_richness ~ Genotype * Wind.Exposure + Error(as.factor(Block)/Wind.Exposure), data= wind_arthropods_May_11_2013_subset, qr=TRUE)
summary(aov_richness)

# use linear mixed effects model. Create nested models and fit with likelihood ration test. Code taken from presentation at: http://www.stat.wisc.edu/courses/st572-larget/Spring2007/handouts16-4.pdf

lmer_richness_null <- lmer(arthropod_richness ~ 1 + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset)

lmer_richness_genotype <- lmer(arthropod_richness ~ Genotype + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset)

lmer_richness_exposure <- lmer(arthropod_richness ~ Wind.Exposure + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset)

lmer_richness_genotype_exposure <- lmer(arthropod_richness ~ Genotype + Wind.Exposure + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset)

lmer_richness_full <- lmer(arthropod_richness ~ Genotype * Wind.Exposure + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset)

anova(lmer_richness_null, lmer_richness_genotype) # suggests no signficant effect of Genotype
anova(lmer_richness_null, lmer_richness_exposure) # SUGGESTS A SIGNIFICANT EFFECT OF WIND EXPOSURE.
anova(lmer_richness_exposure, lmer_richness_genotype_exposure) # adding Genotype as a main effect doesn't improve model fit
anova(lmer_richness_exposure, lmer_richness_genotype_exposure, lmer_richness_full) # adding genotype as main effect or its interaction doesn't improve model fit.

# Linear mixed effects model suggests that wind exposure dominates the effects of plant genotype on arthropod richness. This effect is persistent with log transformed richness data.  Now lets try and fit a genearlized mixed effects model

glmer_richness_null <- glmer(arthropod_richness ~ 1 + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset, family=poisson)

glmer_richness_exposure <- glmer(arthropod_richness ~ Wind.Exposure + (1 | Wind.Exposure:Block), data=wind_arthropods_May_11_2013_subset, family=poisson)

anova(glmer_richness_null, glmer_richness_exposure)

# The Generalized mixed effects model assuming a poisson distribution shows even stronger support for a wind effect. 







