########################################################################################
##      Ant-Aphid by Genotype Experiment data analysis
########################################################################################

### Play around with "date" as a random effect (can I do this or am I violating some assumption?)
### Need to model aphid population growth, following Agrawal 2004
### Need to get herbivores organized.  Look at leaf bundlers and leaf tip folders, as well as homopterans.


## set working directory
setwd("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment")

## load libraries
library("lattice")
library("lme4")

### upload data
# experimental design setup
setup <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Ant-Aphid_Experiment_Setup.csv")

# arthropod data. NEED TO DOUBLE CHECK!!!
raw <- read.csv("~/Desktop/WillowProject/Ant-Aphid Experiment/Ant-Aphid_Data_v2.csv", stringsAsFactors=F)

# dissection data
LTFdata <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Arthropod Dissections/Leaf Miner Dissections - Ant-Aphid Experiment 2012.csv", stringsAsFactors=F)

# plant trait data
trich <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Plant Traits/Trichome Density - Ant-Aphid Experiment 2012.csv")
DryLeafWt <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Plant Traits/Leaf Weights - Dry - Ant-Aphid Experiment 2012.csv")
WetLeafWt <- read.csv("~/Desktop/Lanphere Experiments/Ant-Aphid Experiment/Plant Traits/Ant-aphid Plant Trait Data.csv") # also has percent browning, but I don't know how informative it will be since we didn't take leaves from plants that had leaves that were substantially browned.

# merge and manage plant trait data
setup$Collection.No. <- with(setup, paste(Block,Plant.Position, sep="."))
setup$Zeroes <- seq(0,0,0) # create column with only zeroes
aaTrait <- merge(setup, trich, by="Collection.No.")
aaTrait <- merge(aaTrait, DryLeafWt, by="Collection.No.")
aaTrait <- merge(aaTrait, WetLeafWt, by.x = "Collection.No.", by.y="Collection.Number")
aaTrait <- subset(aaTrait, select = c(Genotype, Block, Distant.to.Ant.Mound, Plant.Position, Collection.No., Aphids.or.No.Aphids, Wet.Leaf.Weight, Percent.Browned, Trichome.Density.x, Dry.Leaf.Weight.x))
aaTrait$LDMC <- with(aaTrait, Dry.Leaf.Weight.x/Wet.Leaf.Weight)
aaTrait$WC <- with(aaTrait, (Wet.Leaf.Weight-Dry.Leaf.Weight.x)/Dry.Leaf.Weight.x)

# merge and manage Caloptilia dissection data
LTFdata$Collection.No. <- with(LTFdata, paste(Block,Plant, sep="."))
aaLTF <- merge(setup, LTFdata, by="Collection.No.")
emptyLTF <- grep("empty", aaLTF$Leaf.Tip.Fold)
aaLTF$Leaf.Tip.Fold <- replace(aaLTF$Leaf.Tip.Fold, emptyLTF,"LTFlarv") # changed all entries with "empty" to "LTFlarv" to indicate they survived.
aaLTF$LTFsurvive <- with(aaLTF, ifelse(aaLTF$Leaf.Tip.Fold == "LTFlarv", 1, 0))

############################ Exploratory Data Analysis ##############################
### Plant Traits
# leaf dry matter content (LDMC). Appears to be a major outlying point.
xyplot(LDMC~Distant.to.Ant.Mound+Aphids.or.No.Aphids, group = Genotype, aaTrait, type="l")
with(aaTrait,interaction.plot(Distant.to.Ant.Mound,Genotype,LDMC))
with(aaTrait, interaction.plot(Distant.to.Ant.Mound,Aphids.or.No.Aphids,LDMC))
with(aaTrait, interaction.plot(Aphids.or.No.Aphids,Genotype,LDMC))
LDMC.glm <- glm(LDMC~Genotype*Distant.to.Ant.Mound*Aphids.or.No.Aphids,aaTrait,family="gaussian")
summary(LDMC.glm)
anova(LDMC.glm,test="F") # appears to be a marginal main effect of distance to ant mound

# water content (WC).  Appears to be a major outlying point.
xyplot(WC~Distant.to.Ant.Mound+Aphids.or.No.Aphids,group=Genotype,aaTrait, type ="l")
WC.glm <- glm(WC~Genotype*Distant.to.Ant.Mound*Aphids.or.No.Aphids,aaTrait,family="gaussian")
summary(WC.glm)
anova(WC.glm,test="F") # wow...apparently there is a significant 3-way interaction

# trichome density. Lines are going crazy...
xyplot(Trichome.Density.x~Distant.to.Ant.Mound+Aphids.or.No.Aphids,group=Genotype,aaTrait, type = "l")
tri.glm <- glm(Trichome.Density.x~Genotype*Distant.to.Ant.Mound*Aphids.or.No.Aphids,aaTrait,family="poisson")
summary(tri.glm)
anova(tri.glm,test="Chisq") # wow...apparently there is ANOTHER significant 3-way interaction.  This seems a bit fishy to me...

########## Note the sampling bias in this data set (use table function below)!  Not equal samples for all genotypes, by any means.  It still may be okay though...
with(aaTrait, table(Genotype))
with(aaTrait, table(Genotype,Distant.to.Ant.Mound,Aphids.or.No.Aphids)) 

### Caloptilia LTF larva abundance. Remember that on July 4 I surveyed LTFs and for a number of them I found parasitoid cocoons, which I believe I removed.  In other words, I think I have two sets up data on LTF abundance and parasitism (will be interesting to see whether the genotype patterns hold up)
######### Still need to remove dead plants from analysis...Also need to include Mature.Green.Leaves as a covariate to measure "density" rather than simply "abundance".
aggLTF <- aggregate(LTFsurvive~Collection.No.+Genotype+Distant.to.Ant.Mound+Aphids.or.No.Aphids, FUN=sum,aaLTF) 
aggLTF <- merge(setup,aggLTF,all.x=T)
aggLTF$LTFsurvive[is.na(aggLTF$LTFsurvive)] <- 0 # change all of the NA's in the column to zero
with(aggLTF, table(Genotype,LTFsurvive) )
with(aggLTF, table(Distant.to.Ant.Mound,LTFsurvive))
with(aggLTF, table(Aphids.or.No.Aphids,LTFsurvive))
aaLTF.glm <- glm(LTFsurvive~Genotype*Aphids.or.No.Aphids*Distant.to.Ant.Mound + Block,family="poisson",aggLTF) # if I try to include interactions, then I receive the warning "fitted rates numerically 0 occurred".
summary(aaLTF.glm)
anova(aaLTF.glm,test="Chisq") # genotype, but no other factor, has a significant effect on LTF larva abundance.  However, if I include ALL interactions, I get a highly signficant genotype by distance to ant mound interaction (as well as a warning message).  Don't know what to do with this yet.  Maybe follow what Mooney did and only examine 2-way interactions.  This may be appropriate if I don't find a 3-way interaction.

# Look for mechanism influencing Leaf tip fold abundance. 
AAtraitLTF <- merge(aggLTF,aaTrait)
subAAtraitLTF <- subset(AAtraitLTF, select = c(LTFsurvive, Trichome.Density.x,LDMC,WC,Percent.Browned))
splom(subAAtraitLTF)
cor(subAAtraitLTF) # unusual apparent correlation with LDMC (but doesn't appear to be one in the graph)
with(subAAtraitLTF, cor.test(LTFsurvive,LDMC)) # highly significant positive correlation.  Remember though that my data is NON-NORMAL.  May be suspect, wind data suggests there is a slight NEGATIVE correlation.  Just goes to show you the potential problem with fishing for data.  Also, LDMC does not appear to be genetically determined.
with(subAAtraitLTF, cor.test(LTFsurvive,Percent.Browned)) # no correlation, but remember that I avoided substantiall browned leaves for the ant-aphid experiment.


## check the data
str(raw)

## omit dead trees
raw$dead <- ifelse(raw$Notes == "dead", 1,0)

## new data columns
herbivores <- with(raw, Spittlebugs+Leaf.Tip.Folds..pre.July.14.category.+Leaf.Hoppers) # new data column with abundance of all herbivores added together.
ants <- with(raw, Red.Ants+Black.Ants) # new data column with abundance of all ants added together

## aggregated data
SumDead <- aggregate(dead ~ Unique.Number+Genotype+Distance.to.Ant.Mound+Aphid.Treatment+Block, FUN=sum,data=raw) #method has not reliably identified all dead trees.

strAvgHerb <- aggregate(herbivores ~ Unique.Number + Genotype + Distance.to.Ant.Mound + Aphid.Treatment + Block, FUN=mean, data=raw)

SumSpid <- aggregate(Spiders ~ Unique.Number + Genotype + Distance.to.Ant.Mound + Aphid.Treatment + Block, FUN=sum, data=raw)

SumRedAnts <- aggregate(Red.Ants ~ Unique.Number + Genotype + Distance.to.Ant.Mound + Aphid.Treatment + Block, FUN=sum, data=raw)

SumAphids <- aggregate(Green.Aphids ~ Unique.Number + Genotype + Distance.to.Ant.Mound + Aphid.Treatment + Block, FUN=sum, data=subset(raw, Aphid.Treatment = "aphid"))

MaxLTF <- aggregate(Leaf.Tip.Folds..pre.July.14.category.~Unique.Number + Genotype + Distance.to.Ant.Mound + Aphid.Treatment + Block, FUN=max, data=raw)

## data currently in use for analysis
data <- cbind(SumAphids,SumDead$dead)
data <- subset(data, SumDead$dead < 1) 
data <- subset(data, Aphid.Treatment == "aphid")

## Effect on Aphid Abundance
bwplot(Green.Aphids ~ as.factor(Genotype), data=data)
xyplot(Green.Aphids ~ as.numeric(Distance.to.Ant.Mound), type = "a", data=data)
with(data, interaction.plot(Distance.to.Ant.Mound,Genotype,Green.Aphids, col=1:10, ylab="Green Aphid Abundance"))
glm.aphid <- glm(Green.Aphids ~ Genotype*Distance.to.Ant.Mound, data, family=quasipoisson)
summary(glm.aphid) # note that the dispersion factor is incredibly high (likely due to high aphid mortality.  Maybe I should convert this into aphid population growth rates...)
anova(glm.aphid, test="Chisq")

## Effect on Red Ant abundance
bwplot(Red.Ants ~ as.factor(Genotype), data=data)
bwplot(Red.Ants ~ as.factor(Aphid.Treatment), data=data)
xyplot(Red.Ants ~ as.numeric(Distance.to.Ant.Mound), type = "a", data=data)
xyplot(Red.Ants ~ as.numeric(Distance.to.Ant.Mound), group = Aphid.Treatment, type = "a", data=data, auto.key=T)
xyplot(Red.Ants ~ as.numeric(Distance.to.Ant.Mound), group = Genotype, type = "a", data=data, auto.key=T)
xyplot(Red.Ants ~ as.numeric(Distance.to.Ant.Mound) + as.factor(Aphid.Treatment), group = Genotype, type = "a", data=data)

# Mixed effect model
(RedAnt.full <- glmer(Red.Ants ~ as.numeric(Distance.to.Ant.Mound)*as.factor(Aphid.Treatment)*as.factor(Genotype) + (1|Block), family = poisson, data=data)) # Full model. Warning message that glm fitted rates numerically 0 occurred

(RedAnt.null <- glmer(Red.Ants ~ (1|Block), family = poisson, data=data))

## Effects on leaf tip folder abundance
# visualize leaf tip folder abundance
bwplot(Leaf.Tip.Folds..pre.July.14.category.~Genotype, data=MaxLTF)
xyplot(Leaf.Tip.Folds..pre.July.14.category.~as.factor(Aphid.Treatment), group = Genotype, type = "a", data=MaxLTF)
xyplot(Leaf.Tip.Folds..pre.July.14.category.~as.numeric(Distance.to.Ant.Mound), group = Genotype, type = "a", data=MaxLTF)

(LTF.full <- glmer(Leaf.Tip.Folds..pre.July.14.category.~as.numeric(Distance.to.Ant.Mound)*as.factor(Aphid.Treatment)*as.factor(Genotype) + (1|Block), family = poisson, data=MaxLTF))

(LTF.geno <- glmer(Leaf.Tip.Folds..pre.July.14.category.~as.factor(Genotype) + (1|Block), family = poisson, data=MaxLTF))

(LTF.null <- glmer(Leaf.Tip.Folds..pre.July.14.category.~ (1|Block), family = poisson, data=MaxLTF))


## Visualize the Spider data
bwplot(Spiders ~ as.factor(Genotype), data=data) # effect of distance to ant mound and genotype on herbivore abundance
xyplot(Spiders ~ as.factor(Aphid.Treatment) + as.numeric(Distance.to.Ant.Mound), group = Genotype, type = "a", data=data) # effect of genotype's interaction with aphid treatment and distance to ant mound on spider abundance

## Effect of treatments on Spider abundance. Genotype alone as well as Genotype's separate interactions with aphid treatment and distance to ant mound appear to be the best models
(spid.full <- glmer(Spiders ~ as.numeric(Distance.to.Ant.Mound)*as.factor(Aphid.Treatment)*as.factor(Genotype) + (1|Block), family = poisson, data=data)) # Full model

(spid.aphidDist <- glmer(Spiders ~ as.factor(Aphid.Treatment)*as.numeric(Distance.to.Ant.Mound) + (1|Block), family = poisson, data=data)) # aphid by distance to ant mound interaction

(spid.genAphid <- glmer(Spiders ~ as.factor(Aphid.Treatment)*as.factor(Genotype) + (1|Block), family = poisson, data=data)) # aphid by genotype interaction

(spid.genDist <- glmer(Spiders ~ as.numeric(Distance.to.Ant.Mound)*as.factor(Genotype) + (1|Block), family = poisson, data=data)) # distance by genotype interaction

(spid.geno <- glmer(Spiders ~ as.factor(Genotype) + (1|Block), family = poisson, data=data)) # only genotype

(spid.dist <- glmer(Spiders ~ as.numeric(Distance.to.Ant.Mound) + (1|Block), family = poisson, data=data)) # only distance to ant mound

(spid.aphid <- glmer(Spiders ~ as.factor(Aphid.Treatment) + (1|Block), family = poisson, data=data)) # only aphid treatment

(spid.null <-  glmer(Spiders ~ (1|Block), family = poisson, data=data)) # null model





## Effect of date and genotype on herbivores
xyplot(herbivores ~ as.numeric(relative.Date), group = Genotype, type = "a", data=raw)

## Effect of data and ant mound distance on herbivores
xyplot(herbivores ~ as.numeric(relative.Date), group = Distance.to.Ant.Mound, type = "a", data=raw)

# testing out a mixed effects models
lme1 <- lmer(herbivores ~ as.factor(Genotype) + (as.factor(Genotype) | Block), data = data)
lme2 <- lmer(herbivores ~ as.factor(Genotype) + (1 | Block), data = data)
anova(lme1, lme2)

################### OLD ################################
## added data
data$PredAbund <- with(data, X._Spiders+X._Fly_Predators+X._Chalcid_Parasitoids + X._Other_Predators)
data$OtherHerbAbund <- with(data, X._Beetles+X._Spittlebugs+X._Leaf_Tip_Folds+X._Leaf_Bundlers+X._Tent_Mines + X._Leaf_Edge_Roller+X._Leaf_Hopper+X._Other_Herbivores)
data$AntAbund <- with(data, X._Red_Ants+X._Black_Ants)

## aphid data only
aphid <- subset(data, Aphid_Treatment == "aphid") 
aphid.NoG <- subset(aphid, Genotype != "G")

## cuttings with at least one aphid
OnePlusAphid <- subset(aphid, X._Green_Aphids > 0)

## Exploratory data analysis
hist(OnePlusAphid$X._Green_Aphids, right =T)
hist(aphid$X._Red_Ants, right =T)
hist(aphid$X._Black_Ants, right = T)

## Effect of Genotype on aphids. GENOTYPE G APPEARS TO BE DRIVING THE GENOTYPE EFFECT ON APHID ABUNDANCE BECAUSE THERE WERE NEVER ANY APHIDS ON G.
plot(X._Green_Aphids~as.factor(Genotype), data=OnePlusAphid)
glm.Aphids <- glm(X._Green_Aphids~as.factor(Genotype), data=OnePlusAphid, family=quasipoisson)
sum.Aphids <- summary(glm.Aphids)
anova(glm.Aphids, test = "F")

## Effect of ant mound distance and genotype on aphids
xyplot(X._Green_Aphids~as.numeric(Distance_to_Ant_Mound), group = Genotype, type = "a", data=aphid)
glm.Inter <- glm(X._Green_Aphids ~ as.factor(Genotype) * as.numeric(Distance_to_Ant_Mound), data = aphid, family = quasipoisson)
sum.Inter <- summary(glm.Inter)
anova(glm.Inter, test = "F")

## Effect of ant mound distance on red ants
xyplot(X._Red_Ants ~ Distance_to_Ant_Mound, groups = Aphid_Treatment, data = data)
glm.Ant <- glm(X._Red_Ants ~ Distance_to_Ant_Mound * Aphid_Treatment, data = data, family = quasipoisson)
sum.Ant <- summary(glm.Ant)

## Effect of red ants on black ant abundance
xyplot(X._Black_Ants ~ X._Red_Ants, data=data)

## Effect of red ants on abundance of all predators
xyplot(PredAbund ~ X._Red_Ants, data=data)

## Effect of black ants on abundance of all predators
xyplot(PredAbund ~ X._Black_Ants, data=data)

## Effect of red ants on herbivore abundance
xyplot(OtherHerbAbund~X._Red_Ants, data=data)

## Effect of black ants on herbivore abundance
xyplot(OtherHerbAbund~X._Black_Ants, data=data)

## Effect of red ants on aphid abundance
xyplot(X._Green_Aphids ~ X._Red_Ants, data=data)

## Effect of black ants on aphid abundance
xyplot(X._Green_Aphids ~ X._Black_Ants, data=data)

## Effect of all ants on aphid abundance
xyplot(X._Green_Aphids ~ AntAbund, data=data)
with(data, cor.test(X._Green_Aphids, AntAbund))

## Effect of data on aphid abundance
xyplot(X._Green_Aphids ~ as.factor(Date), data=data)
