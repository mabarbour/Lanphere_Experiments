
###### upload required libraries
library(ggplot2)
library(BiodiversityR)
library(plyr)
library(lme4)

##### upload and manage data

plant_growth_May12_2013 <- read.csv("~/Documents/Lanphere_Experiments/data/wind_plant_growth_data_May12-13_2013.csv", skip=1)

str(plant_growth_May12_2013)

plant_growth_May12_2013$Block <- as.factor(plant_growth_May12_2013$Block)
plant_growth_May12_2013$Total_shoot_growth <- rowSums(plant_growth_May12_2013[ ,9:29])
plant_growth_May12_2013$Height_Difference <- plant_growth_May12_2013$Height - plant_growth_May12_2013$Original_Height 

plant_growth_May12_2013 <- plant_growth_May12_2013[ ,-30] # remove Notes section

plant_growth_May12_2013_not_dead <- plant_growth_May12_2013[-which(plant_growth_May12_2013$Dead_for_plant_survey == 1), ]

##### Exploratory Data analysis

total_shoot_growth <- ggplot(plant_growth_May12_2013_not_dead, aes(x=Wind.Exposure, y=Total_shoot_growth, color=Genotype))
total_shoot_growth + geom_boxplot()

total_shoot_growth_means <- ddply(plant_growth_May12_2013_not_dead, .(Wind.Exposure,Genotype), summarise, val=mean(Total_shoot_growth))

total_shoot_growth + geom_point(data=total_shoot_growth_means, aes(y=val)) + geom_line(data=total_shoot_growth_means, aes(y=val, group=Genotype)) # appears to be a strong effect of genotype on shoot growth

height_difference <- ggplot(plant_growth_May12_2013_not_dead, aes(x=Wind.Exposure, y=Height_Difference, color=Genotype))
height_difference + geom_boxplot()

height_difference_means <- ddply(plant_growth_May12_2013_not_dead, .(Wind.Exposure,Genotype), summarise, val=mean(Height_Difference))

height_difference + geom_point(data=height_difference_means, aes(y=val)) + geom_line(data=height_difference_means, aes(y=val, group=Genotype)) # most of the Genotypes have dropped in their average height

### Seems like wind pruning is having a dramatic effect on plant height which is influencing arthropod richness over plant biomass.

#####  Model the data

lmer_total_shoot_growth_null <- lmer(Total_shoot_growth ~ 1 + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_total_shoot_growth_genotype <- lmer(Total_shoot_growth ~ Genotype + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_total_shoot_growth_expsoure <- lmer(Total_shoot_growth ~ Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_total_shoot_growth_genotype_expsoure <- lmer(Total_shoot_growth ~ Genotype + Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_total_shoot_growth_genotype_exposure_interaction <- lmer(Total_shoot_growth ~ Genotype*Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)

anova(lmer_total_shoot_growth_null, lmer_total_shoot_growth_genotype) # signficant effect of Genotype
anova(lmer_total_shoot_growth_null, lmer_total_shoot_growth_expsoure) # no signifcant effect of Wind exposure
anova(lmer_total_shoot_growth_genotype, lmer_total_shoot_growth_genotype_expsoure) # no signficant effect compared to wind + genotype model
anova(lmer_total_shoot_growth_genotype, lmer_total_shoot_growth_genotype_exposure_interaction) # no interactive effect.

### Shoot growth depends on plant genotype and not wind expsoure

lmer_height_difference_null <- lmer(Height_Difference ~ 1 + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_height_difference_genotype <- lmer(Height_Difference ~ Genotype + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_height_difference_exposure <- lmer(Height_Difference ~ Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_height_difference_genotype_exposure <- lmer(Height_Difference ~ Genotype + Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)
lmer_height_difference_genotype_exposure_interaction <- lmer(Height_Difference ~ Genotype*Wind.Exposure + (1 | Wind.Exposure:Block), data=plant_growth_May12_2013_not_dead)

anova(lmer_height_difference_null, lmer_height_difference_genotype) # significant effect of Genotype on height difference
anova(lmer_height_difference_null, lmer_height_difference_exposure) # signficant effect of wind.exposure
anova(lmer_height_difference_genotype, lmer_height_difference_exposure) # exposure model isn't significantly better than genotype
anova(lmer_height_difference_genotype, lmer_height_difference_genotype_exposure) # main effect of Genotype and wind exposure is so far the best
anova(lmer_height_difference_genotype_exposure, lmer_height_difference_genotype_exposure_interaction) # interaction model isn't great.

####  Appears to be a main effect of both Genotype and Wind exposure, although wind exposure appears to have a larger effect on the degree of wind pruning (change in original height)

qplot(plant_growth_May12_2013_not_dead$Height_Difference) #, geom="density") # distribution appears bimodal...may be an issue for the models


qplot(plant_growth_May12_2013_not_dead$Total_shoot_growth) # heavily right skewed distribution
qplot(log(plant_growth_May12_2013_not_dead$Total_shoot_growth + 1)) # log transformatin make the data more symmetrical

qplot(plant_growth_May12_2013_not_dead$Height)

### Will need to evaluate whether I need to use any special distributions for modelling these data.

#######  In summary, Genotype appears to have the strongest effect on biomass, however, the wind pruning alters the distribution of that biomass to lower on the plant.  This may actually mean that those plants have lower herbivore loads...


