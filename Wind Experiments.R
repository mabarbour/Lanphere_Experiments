############## Wind-by-Genotype Experiment####################


######### Need to update data sheet with twisty gall, stem gall, sawfly, etc. abundances
######### Interesting pattern: wind exposed plants attract light green aphids [which is a pattern I thoght I was observing in the field as well].  Need to think about the story behind this...[would need to remove aphids to test for interactive effects on senescence]
######## Calculate degree of stress? Difference between average number of mature leaves per shoot on unexposed genotype and observed number of mature leaves per shoot.  I could also use the percent browned leaf as well as whether or not a leaf was even collected from the plant. Or difference between average number of mature leaves on unexposed genotype and observed number of mature leaves.

# set working directory
setwd("~/Desktop/Lanphere Experiments/Wind Experiment")

# load libraries
library("lattice")
library("lme4")
library("nlme")
library("reshape2")

# load data
WindData <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Visual Arthropod Surveys & Plant Architecture Data/Wind Plant Trait and Arthropod data.csv")

Wind.leaf.trait.data <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Plant Traits/Wind leaf trait data.csv")
trichomes <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Plant Traits/Trichome Density - Wind Experiment 2012.csv")
DryLeafWt <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Plant Traits/Leaf Weights - Dry - Wind Experiment 2012.csv")

WindLeafMine <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Arthropod Dissections/Leaf Miner Dissections - Wind Experiment 2012.csv", stringsAsFactors=F)

WindLeafTiesOthers <- read.csv("~/Desktop/Lanphere Experiments/Wind Experiment/Arthropod Dissections/Leaf Ties - Wind Experiment 2012.csv")

### manage data
# manage wind data
WindData$Wind.Exposure <- factor(x = WindData$Wind.Exposure, levels = c("Exposed","Unexposed"), labels = c("e","u"))
WindData$Collection.Number <- with(WindData, paste(Wind.Exposure,Block,sep=""))
WindData$Collection.Number <- with(WindData, paste(Collection.Number,Plant.Position,sep="."))
WindLeafMine$Collection.Number <- with(WindLeafMine, paste(Treatment,Block,sep=""))
WindLeafMine$Collection.Number <- with(WindLeafMine, paste(Collection.Number,Plant, sep="."))
WindLeafTiesOthers$Collection.Number <- with(WindLeafTiesOthers, paste(Treatment,Block,sep=""))
WindLeafTiesOthers$Collection.Number <- with(WindLeafTiesOthers, paste(Collection.Number,Plant, sep="."))

# merge and manage data sets
wind <- merge(WindData,Wind.leaf.trait.data,by.x="Collection.Number",by.y="Collection.Number")
wind <- merge(wind,trichomes,by.x="Collection.Number", by.y="Collection.No.")
wind <- merge(wind,DryLeafWt, by.x="Collection.Number", by.y="Collection.No.")
#wind <- merge(wind,WindLeafMine, by.x="Collection.Number", by.y="Collection.Number") # noticed that I lost data Leaf miner data when I merged these data set
wind$LDMC <- with(wind, Dry.Leaf.Wt/Wet.Leaf.Weight)
wind$WC <- with(wind, (Wet.Leaf.Weight-Dry.Leaf.Wt)/Dry.Leaf.Wt)

# manage leaf mine parasitism data
windLM <- merge(WindLeafMine, WindData, by="Collection.Number")
empty <- grep("empty", windLM$Leaf.Tip.Fold)
windLM$Leaf.Tip.Fold <- replace(windLM$Leaf.Tip.Fold, empty,"LTFlarv") # changed all entries with "empty" to "LTFlarv" to indicate they survived.
sp.329 <- grep("329", windLM$Leaf.Tip.Fold)
windLM$Leaf.Tip.Fold <- replace(windLM$Leaf.Tip.Fold, sp.329, "sp.329") # changed all entries with "329" to "sp.329" to indicate they were parasitized by sp. 329.
windLM$LTF.survive <- with(windLM, ifelse(Leaf.Tip.Fold == "LTFlarv", 1,0))
windLM$sp329.attack <- with(windLM, ifelse(Leaf.Tip.Fold == "sp.329",1,0))

# manage leaf tie and other parasitism data
windOther <- merge(WindLeafTiesOthers,WindData,by="Collection.Number")
windOther$live319 <- with(windOther, ifelse(Contents == "living 319", 1,0))
windOther$live320 <- with(windOther, ifelse(Contents == "living 320", 1,0))
windOther$live321 <- with(windOther, ifelse(Contents == "living 321", 1,0))
windOther$live322 <- with(windOther, ifelse(Contents == "living 322", 1,0))


WindData$ShootTotal <- with(WindData, X1.shoot.M+X2.shoot.M+X3.shoot.M+X4.shoot.M+X5.shoot.M+X6.shoot.M+X7.shoot.M+X8.shoot.M+X9.shoot.M+X10.shoot.M)

WindData$LeavesPerSL <- with(WindData,Mature.Green.Leaves/ShootTotal)

WindData$SpidTot <- with(WindData, Spiders.black.with.yellow.legs+Spiders.Other)

WindData$lacewingEggs <- WindData$Parisatoid.cocoons.on.Aphids

WindData$LH <- with(WindData, Leaf.Hoppers.Blue.Brown+Leaf.Hoppers.Camo+LH.nymph.blk.ylw+LH.nymph.other)

WindData$AphidLoad <- with(WindData,Light.Green.Aphids/Mature.Green.Leaves)



################### Exploratory Data Analysis #############################

# GxE effect on Percent browning.  Probably need to redo the analysis by treating them as random effects
with(wind,interaction.plot(Wind.Exposure,Genotype,Percent.Browned, col=1:10)) # clear interactive effect of wind exposure on percent of leaf that is browned.
brown.glm <- glm(Percent.Browned~Wind.Exposure*Genotype,wind,family="gaussian")
summary(brown.glm)
anova(brown.glm,test="Chisq") #strong effect of wind, and a slight interactive effect

# GxE effect on Trichome density. Probably need to redo the analysis by treating them as random effects
with(wind,interaction.plot(Wind.Exposure,Genotype,Trichome.Density)) # clear genotype effect on trichome density
tri.glm <- glm(Trichome.Density~Wind.Exposure*Genotype,wind,family="quasipoisson")
summary(tri.glm)
anova(tri.glm,test="Chisq") # large effect of plant genotype

# GxE effect on Leaf Dry Matter Content (LDMC)
with(wind,interaction.plot(Wind.Exposure,Genotype,LDMC))
LDMC.glm <- glm(LDMC~Wind.Exposure*Genotype,wind,family="gaussian")
summary(LDMC.glm)
anova(LDMC.glm,test="F") # appears to be a marginal main effect of genotype and wind exposure

# GxE effect on water content (WC) of leaves
with(wind, interaction.plot(Wind.Exposure,Genotype,WC))
WC.glm <- glm(WC~Wind.Exposure*Genotype,wind,family="gaussian")
summary(WC.glm)
anova(WC.glm, test="F") # appears to be a marginal main effect of wind exposure, but higher water content at exposed sites???

# GxE effect on leaf tip folder survival. Still not accounting for other stages or parasitism. Concerned about analyzing "glmer" models.  Not accounting for density (per leaf), simply abundance.
aggLTF <- aggregate(LTF.survive~Collection.Number+Genotype+Treatment,FUN=sum, windLM)
aggLTF <- merge(WindData,aggLTF, all.x=T)
aggLTF$LTF.survive[is.na(aggLTF$LTF.survive)] <- 0 # change all of the NA's in the column to zero
with(aggLTF, interaction.plot(Treatment,Genotype,LTF.survive, col=1:10)) # abundance of LTF appears to increase at unexposed sites.
LTF.1 <- glmer(LTF.survive~Wind.Exposure + (1|Genotype) + (1|Block),family="poisson",aggLTF)
summary(LTF.1) 

LTF.2 <- glmer(LTF.survive~Wind.Exposure + (Block|Genotype), family="poisson", aggLTF)
anova(LTF.1,LTF.2) # suggests that there is a block by genotype interaction?

LTF.null <- glmer(LTF.survive~(Block|Genotype), family="poisson",aggLTF)
anova(LTF.null,LTF.1) # suggest that wind exposure has a vary strong effect

summary(LTF.2)

LTF.random <- glmer(LTF.survive~(1|Wind.Exposure)+(1|Genotype)+(1|Block), family="poisson", aggLTF)
summary(LTF.random) # block has no variance explained???

LTF.fixed2 <- glmer(LTF.survive~Wind.Exposure*Genotype + (1|Block), family="poisson", aggLTF)
summary(LTF.fixed)

LTF.fixed1 <- glmer(LTF.survive~Wind.Exposure+Genotype+(1|Block), family="poisson", aggLTF)
summary(LTF.fixed1)

anova(LTF.fixed,LTF.fixed1) # again, marginally significant effect of including the interaction term

glm.LTF <- glm(LTF.survive~Wind.Exposure*Genotype+Block+Mature.Green.Leaves,family="poisson",aggLTF)
summary(glm.LTF) # probelm that most of my p-value estimates are close to 1???  I feel like I have a zero-inflated dataset, but there doesn't appear to be overdispersion...Actually, the residuals don't look good at all.
anova(glm.LTF,test="Chisq") # marginally significant wind by genotype interaction.  Stronger main effects of genotype relative to wind exposure.  Density-mediated effect by some plant trait (leaf abundance by itself doesn't explain much).

# look for correlations between plant traits and LTF abundance
traitLTF <- merge(wind,aggLTF)
subTraitLTF <- subset(traitLTF, select = c(LTF.survive,ShootTotal,Max.Height, Trichome.Density,LDMC,WC,Percent.Browned,Mature.Green.Leaves)) #,Genotype,Wind.Exposure,Block))
splom(subTraitLTF) # may be a bit of relationship with the percent of leaf that is browned (indicator of water stress...)
cor(subTraitLTF) # why is nothing showing up with Mature.Green.Leaves???  Note that LDMC and WC are highly correlated (as expected)
with(subTraitLTF, cor.test(LTF.survive,Percent.Browned)) # significant negative correlation between these two parameters
with(subTraitLTF, cor.test(LTF.survive, ShootTotal)) # marginally significant positive correlation between these two parameters
with(subTraitLTF, cor.test(LTF.survive,Mature.Green.Leaves)) # significant positive correlation between these two parameters
     
# GxE effect on parasitoid sp. 329 abundance.  Right now, I think my abundances are too low for at least the LTFs.
agg329 <- aggregate(sp329.attack~Collection.Number+Genotype+Treatment,FUN=sum,windLM)
with(agg329, interaction.plot(Treatment,Genotype,sp329.attack)) # no clear effect...
with(agg329, table(Treatment,Genotype,sp329.attack))
agg329 <- merge(WindData,agg329, all.x=T)
agg329$sp329.attack[is.na(agg329$sp329.attack)] <- 0 # change all of the NA's in the column to zero

glm.329 <- glm(sp329.attack~Genotype*as.factor(Treatment)+Block,family="poisson",agg329)
summary(glm.329) # problem that estimates are close to 1
anova(glm.329,test="Chisq") # I suspect that my abundances are too low to detect any effects

# explore relationship between abundance of LTF and sp. 329
LTF.329 <- merge(aggLTF,agg329)
aggLTFcorr <- aggregate(LTF.survive ~ Genotype + Wind.Exposure, FUN=mean, LTF.329)
agg329corr <- aggregate(sp329.attack~Genotype+Wind.Exposure,FUN=mean,LTF.329)
LTF.329corr <- merge(aggLTFcorr,agg329corr)

plot(sp329.attack~LTF.survive,LTF.329corr)
with(LTF.329corr, cor.test(sp329.attack,LTF.survive)) # okay, highly significant correlation between parasitoid abundance and abundances of LTFs when they are averaged across genotypes and wind exposure (i.e. 20 points)

# GxE effect on morphospecies 319, 320, 321 abundance. Don't know how the stage structure will work for this population.  But can I argue that it is stage structured?  I do have what appear to be "in-between" stages of larva.
agg319 <- aggregate(live319~Collection.Number+Genotype+Treatment,FUN=sum,windOther)
with(agg319, interaction.plot(Treatment,Genotype,live319))
glm.319 <- glm(live319~Genotype*as.factor(Treatment),family="poisson",agg319)
summary(glm.319)
anova(glm.319, test="Chisq") # I wonder why there is a significant interactive effect...

agg320 <- aggregate(live320~Collection.Number+Genotype+Treatment,FUN=sum,windOther)
with(agg320, interaction.plot(Treatment,Genotype,live320))
glm.320 <- glm(live320~Genotype*as.factor(Treatment),family="poisson",agg320)
summary(glm.320)
anova(glm.320, test="Chisq") # no effect for sp. 320

agg321 <- aggregate(live321~Collection.Number+Genotype+Treatment,FUN=sum,windOther)
with(agg321, interaction.plot(Treatment,Genotype,live321))
glm.321 <- glm(live321~Genotype*as.factor(Treatment),family="poisson",agg321)
summary(glm.321)
anova(glm.321, test="Chisq") # no effect for sp. 321

agg322 <- aggregate(live322~Collection.Number+Genotype+Treatment,FUN=sum,windOther)
with(agg322, interaction.plot(Treatment,Genotype,live322)) 
glm.322 <- glm(live322~Genotype*as.factor(Treatment),family="poisson",agg322)
summary(glm.322)
anova(glm.322, test="Chisq") # no effect for sp. 322

# data sets
data <- WindData
omitDead <- subset(WindData,Dead.!=1)
ltrait <- Wind.leaf.trait.data
WindDissect

# Parasitism.  At present the rate of parasitism didn't vary between treatments, but is simply a function of LTF/LEF density.
with(WindDissect, table(LTF.Parasitized,Treatment))
with(WindDissect, table(LEF.Parasitized,Treatment))


# Percent of collected leaf browned
plot(Percent.Browned~Treatment,ltrait, ylab = "Percent of Leaf Browned") # "w" category was for a leaf from an unknown treatment.

# Survival
with(data,table(Dead.,Genotype,Wind.Exposure))
death.glm <- glm(Dead.~Genotype*Wind.Exposure+Block,family=binomial,data)
summary(death.glm)
anova(death.glm,test="Chisq")

# Shoot Total Length
par(mar=c(5,5,4,2))
with(omitDead, interaction.plot(Wind.Exposure,Genotype,ShootTotal,fun=mean,type="l", col=1:10, legend=F, xlab="",ylab="Shoot Growth (cm)", cex.axis=1.5, cex.lab=1.5, lwd=3)) # shoot growth = sum of all shoot distances measured
shoot.lm <- lm(ShootTotal~Genotype*Wind.Exposure+Block,omitDead) # don't know if this is coded correctly
summary(shoot.lm)
anova(shoot.lm)

# Number of Mature Leaves
with(subset(omitDead,Mature.Green.Leaves > 0), interaction.plot(Wind.Exposure,Genotype,Mature.Green.Leaves,fun=mean, col=1:10, ylab="Number of Mature Leaves")) 
MGL.lm <- lm(Mature.Green.Leaves~Wind.Exposure*Genotype+Block,omitDead)
summary(MGL.lm)
anova(MGL.lm) # interactive effect on the number of mature green leaves. Why?

# Number of Mature Leaves per shoot length 
with(subset(omitDead,LeavesPerSL > 0), interaction.plot(Wind.Exposure,Genotype,LeavesPerSL,fun=mean, col=1:10)) #still not working...
LeavesPerSL.lm <- lm(LeavesPerSL~Genotype*Wind.Exposure+Block,omitDead)
summary(LeavesPerSL.lm)
anova(LeavesPerSL.lm) # Stronger effect of plant genotype

# Plant Height
with(omitDead, interaction.plot(Wind.Exposure,Genotype,Max.Height,fun=mean, col=1:10))
height.lm <- lm(Max.Height~Genotype*Wind.Exposure+Block,WindData)
summary(height.lm)
anova(height.lm)  # relatively stronger effect of host genotype

# Spider Abundance
with(omitDead,interaction.plot(Wind.Exposure,Genotype,SpidTot,fun=mean, col=1:10))
spid.glm <- glm(SpidTot~Wind.Exposure*Genotype+Block,omitDead,family=poisson)
summary(spid.glm)
anova(spid.glm,test="Chisq")

# Light Green Aphid Abundance
with(omitDead,interaction.plot(Wind.Exposure,Genotype,Light.Green.Aphids,fun=mean,col=1:10, ylab = "Light Green Aphid Abundance"))
plot(Light.Green.Aphids~Wind.Exposure,omitDead)
plot(Light.Green.Aphids~Mature.Green.Leaves,omitDead)
plot(Light.Green.Aphids~LeavesPerSL,omitDead)
aphid.glm <- glm(Light.Green.Aphids~Wind.Exposure*Genotype+Block,omitDead,family=quasipoisson)
summary(aphid.glm)
anova(aphid.glm,test="Chisq") # wind expsoure and block have a significant effect on light.green.aphid abundance.
plot(aphid.glm)  # quasipoisson is a much better model
leaves.glm <- glm(Light.Green.Aphids~LeavesPerSL+Max.Height,omitDead,family=quasipoisson)
summary(leaves.glm) # both plant height and number of leaves per shoot length influence light green aphid abundance.
aphid.glm2 <- glm(Light.Green.Aphids~Wind.Exposure+Genotype+Block+LeavesPerSL+Max.Height,omitDead,family=quasipoisson)
summary(aphid.glm2)
anova(aphid.glm2,test="Chisq")

# Light green aphid "load" (number per mature green leaves)
with(omitDead,interaction.plot(Wind.Exposure,Genotype,AphidLoad,fun=mean,col=1:10)) # not really working...
aphidload <- lm(AphidLoad~Wind.Exposure*Genotype,omitDead)
summary(aphidload)
anova(aphidload) # still a significant effect of wind, but not much else.  In other words, the number of mature green leaves doesn't appear to have much consequence on the number of aphids [but perhaps stress?]

# Lacewing eggs
with(omitDead,interaction.plot(Wind.Exposure,Genotype,lacewingEggs,fun=mean,col=1:10))
lacewing.glm <- glm(lacewingEggs~Wind.Exposure*Genotype+Block,omitDead,family=poisson)
summary(lacewing.glm)
anova(lacewing.glm,test="Chisq") # genotype effect on lacewing egg abundance.  However, these abundances are quite low...

# Leafhoppers Blue-Brown. No other leaf hopper morphospecies or nymphs occurred during the survey
with(omitDead,interaction.plot(Wind.Exposure,Genotype,Leaf.Hoppers.Blue.Brown,fun=mean,col=1:10)) # very few in abundance

# Ants
with(omitDead,interaction.plot(Wind.Exposure,Genotype,Black.Ants,fun=mean,col=1:10)) # very few in abundance
