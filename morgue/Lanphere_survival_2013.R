# set working directory
setwd("~/Documents/Lanphere_Experiments/Lanphere_R_scripts_&_Results")

# Upload survival data from Lanphere Experiments as of Sunday April 28, 2013
survival_ant_aphid <- read.csv("~/Documents/Lanphere_Experiments/data/Lanphere_Phenology_2013_ant_aphid_experiment.csv", skip=1)
survival_ant_aphid <- subset(survival_ant_aphid, select=Block:Dead)
survival_ant_aphid$Block <- as.factor(survival_ant_aphid$Block)

survival_wind <- read.csv("~/Documents/Lanphere_Experiments/data/Lanphere_Phenology_2013_wind_experiment.csv", skip=1)
survival_wind <- subset(survival_wind, select=Block:Dead)
survival_wind$Block <- as.factor(survival_wind$Block)

###### Tabulate Survival data

# ant-aphid
sum(survival_ant_aphid$Dead)/300 # over 50% mortality. From my field observations, many of these plants were sprouting shoots fromt their base.

block_survival <- with(survival_ant_aphid, table(Dead, Block)) # rather large mortality in block 4 and relatively less mortality in block 5.
plot(block_survival)

aphid_survival <- with(survival_ant_aphid, table(Dead, Aphid.Treatment))
plot(aphid_survival) # slightly higher survival in aphid treatment

distance_survival <- with(survival_ant_aphid, table(Dead, Ant.Mound.Dist))
plot(distance_survival) # little difference across distance treatments

genotype_survival <- with(survival_ant_aphid, table(Dead, Genotype))
plot(genotype_survival) # appears to be a lot of variation in survival across genotypes

aphiddistI_survival <- with(survival_ant_aphid, table(Dead, Aphid.Treatment, Ant.Mound.Dist)) # no apparently strong interaction...But it appears to show up in the model if genotype is included as a main effect...

# all possible models.  Need to double check that I have specified the model correctly with this complicated split plot design.

glmer_ant_aphid_full <- glmer(Dead ~ Genotype * Aphid.Treatment * Ant.Mound.Dist + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid) # note that right now I'm treating Ant.Mound.Dist as a factor in the random effect portion...It doesn't appear to be able to handle it as a continuous factor.

glmer_null <- glmer(Dead ~ 1 + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_geno <- glmer(Dead ~ Genotype + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_dist <- glmer(Dead ~ Ant.Mound.Dist + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_aphid <- glmer(Dead ~ Aphid.Treatment + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_genodistI <- glmer(Dead ~ Genotype*Ant.Mound.Dist + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_genodistI_aphid <- glmer(Dead ~ Genotype*Ant.Mound.Dist + Aphid.Treatment + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_genoaphidI <- glmer(Dead ~ Genotype*Aphid.Treatment + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_genoaphidI_dist <- glmer(Dead ~ Genotype*Aphid.Treatment + Ant.Mound.Dist + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_aphiddistI <- glmer(Dead ~ Aphid.Treatment*Ant.Mound.Dist + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

glmer_aphiddistI_geno <- glmer(Dead ~ Aphid.Treatment*Ant.Mound.Dist + Genotype + (1 | Block/as.factor(survival_ant_aphid$Ant.Mound.Dist)), family=binomial, data=survival_ant_aphid)

anova(glmer_null, glmer_geno, glmer_aphid, glmer_dist, glmer_genoaphidI, glmer_genoaphidI_dist, glmer_genodistI, glmer_genodistI_aphid, glmer_aphiddistI, glmer_aphiddistI_geno, glmer_ant_aphid_full) # according to AIC, the two best competing models have genotype as an independent effect and included aphid and distance interaction.

anova(glmer_geno, glmer_aphiddistI_geno) # these two models don't appear to be siginficantly different from each other

### Wind experiment
sum(survival_wind$Dead)/200 # ~ 20% mortality.

wind_block_survival <- with(survival_wind, table(Dead, Block)) # rather large mortality in block 4 and relatively less mortality in block 5.
plot(wind_block_survival) # little heterogeneity across blocks

wind_exposure_survival <- with(survival_wind, table(Dead, Wind.Exposure))
plot(wind_exposure_survival) # margingally higher survival in exposed treatments

wind_genotype_survival <- with(survival_wind, table(Dead, Genotype))
plot(wind_genotype_survival) # heterogeneity across genotypes in survival. Genotype S had very high mortality, but that may be because many of their roots were stripped during planting...May need to explore if an effect is still there when this genotype is removed.

# model the data
glmer_wind_full <- glmer(Dead ~ Genotype * Wind.Exposure + (1 | Wind.Exposure:Block), family=binomial, data=survival_wind) # I think this is the appropriate way to specify the random effect.

glmer_null_wind <- glmer(Dead ~ 1 + (1 | Wind.Exposure:Block), family=binomial, data=survival_wind)

glmer_geno_wind <- glmer(Dead ~  Genotype + (1 | Wind.Exposure:Block), family=binomial, data=survival_wind)

glmer_exposure <- glmer(Dead ~  Wind.Exposure + (1 | Wind.Exposure:Block), family=binomial, data=survival_wind)

glmer_geno_exposure <- glmer(Dead ~  Genotype + Wind.Exposure + (1 | Wind.Exposure:Block), family=binomial, data=survival_wind)

anova(glmer_null_wind, glmer_exposure, glmer_geno_wind, glmer_geno_exposure, glmer_wind_full)

# let's see if the genotype effect remains after Genotype S is removed.  May need to go back to notes on whether roots were stripped or not during planting...
noS_wind <- subset(survival_wind, Genotype != "S")

noS_exposure <- glmer(Dead ~ Wind.Exposure + (1 | Wind.Exposure:Block), family=binomial, data=noS_wind)
noS_geno <- glmer(Dead ~ Genotype + (1 | Wind.Exposure:Block), family=binomial, data=noS_wind)

noS_null <- glmer(Dead ~ 1 + (1 | Wind.Exposure:Block), family=binomial, data=noS_wind)

noS_geno_exposure <- glmer(Dead ~ Genotype + Wind.Exposure + (1 | Wind.Exposure:Block), family=binomial, data=noS_wind)

anova(noS_null, noS_exposure, noS_geno, noS_geno_exposure) # when Genotype S is removed, the results are qualitatively the same.





