# Why are there zeros for mature shoot total length???

## load required libraries ----
library(dplyr)
library(reshape2)
library(ggplot2)
library(visreg)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
library(lmerTest)
#library(RLRsim)
library(car)

## load required data set ----
wind.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_trait_df.csv') %>%
  tbl_df() %>%
  mutate(block = as.factor(block),
         Year = as.factor(Year),
         plant.position = as.character(plant.position),
         X = as.factor(X))
glimpse(wind.df)

aa.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/ant_aphid_trait_df.csv') %>%
  tbl_df() %>%
  mutate(fact.Ant.Mound.Dist = as.factor(Ant.Mound.Dist))
glimpse(aa.df)

## Wind: Plant architecture 2012 ----

# explore correlations among architecture traits
arch.traits <- c("Height","all.shoot.total.length","all.shoot.count")

# 2012 
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,arch.traits])
height.outlier <- which(wind.df$Height > 80)  # not a typo according to data entry, but definitely seems to be an outlier. Maybe it was written down incorrectly in the field, because it is almost twice the size of the other plants. I've decided to remove it from this analysis.
scatterplotMatrix(filter(wind.df, Year == "2012")[-height.outlier,
                                                  arch.traits])
scatterplotMatrix(log(filter(wind.df, Year == "2012")[-height.outlier,
                                                      arch.traits])) # logging doesn't appear to linearize these relationships anymore.
corr.test(filter(wind.df, Year == "2012")[-height.outlier,arch.traits])

# PCA: 2012
arch.pca.2012.df <- filter(wind.df, Year == "2012")[-height.outlier, ] %>%
  #filter(Year == "2012") %>%
  select(block, genotype, treatment, 
         Height, all.shoot.total.length, all.shoot.count) %>%
  na.omit()

arch.PCA.2012 <- princomp(arch.pca.2012.df[ ,arch.traits], cor = TRUE)
summary(arch.PCA.2012) # only first component has an eigenvalue > 1
plot(arch.PCA.2012)
arch.PCA.2012$loadings # positive values of PC1 indicate higher levels of productivity.
arch.PCA.2012$scores[ ,"Comp.1"]

arch.pca.2012.df <- mutate(arch.pca.2012.df, arch.PC1 = arch.PCA.2012$scores[ ,"Comp.1"]) 

hist(arch.pca.2012.df$arch.PC1) # normally distributed

## Analysis

# arch PC1 2012 
pc1.2012.lmer <- lmer(arch.PC1 ~ genotype*treatment +
                        (1|block) + (1|block:treatment),
                      data = arch.pca.2012.df)
plot(pc1.2012.lmer) # looks pretty good
anova(pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(pc1.2012.lmer, xvar = "genotype", by = "treatment")

# plant height
height.2012 <- lmer(Height ~ genotype*treatment + 
                      (1|block) + (1|block:treatment), 
                    data = arch.pca.2012.df) 
plot(height.2012) # not great
anova(height.2012, ddf = "Kenward-Roger")
visreg(height.2012, xvar = "genotype", by = "treatment")

# shoot count
shoot.count.2012 <- lmer(all.shoot.count ~ genotype*treatment +
                           (1|block) + (1|block:treatment), 
                         data = arch.pca.2012.df)
plot(shoot.count.2012) # okay
anova(shoot.count.2012, ddf = "Kenward-Roger")
visreg(shoot.count.2012, xvar = "genotype", by = "treatment")

# shoot length
shoot.length.2012 <- lmer(all.shoot.total.length ~ genotype*treatment +
                            (1|block) + (1|block:treatment), 
                          data = arch.pca.2012.df)
plot(shoot.length.2012) # okay
anova(shoot.length.2012, ddf = "Kenward-Roger") # GxE 
visreg(shoot.length.2012, xvar = "genotype", by = "treatment")

## Wind: Leaf quality 2012 ----

# explore correlations among architecture traits
LQ.traits <- c("leaf_WC","leaf_trichome.density")

# 2012 
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,LQ.traits])
WC.outlier <- which(wind.df$leaf_WC > 5)  # I've decided to remove it from this analysis, because it is unreasonably large and was likely a measurement error
scatterplotMatrix(filter(wind.df, Year == "2012")[-WC.outlier,
                                                  LQ.traits])
scatterplotMatrix(log(filter(wind.df, Year == "2012")[-WC.outlier,
                                                      LQ.traits]+1)) # logging doesn't appear to linearize these relationships anymore.
corr.test(filter(wind.df, Year == "2012")[-WC.outlier,LQ.traits])

# PCA: 2012
LQ.pca.2012.df <- filter(wind.df, Year == "2012")[-WC.outlier, ] %>%
  #filter(Year == "2012") %>%
  select(block, genotype, treatment, 
         leaf_WC, leaf_trichome.density) %>%
  mutate(log_trichome.density = log(leaf_trichome.density+1),
         log_WC = log(leaf_WC)) %>%
  na.omit()

log.LQ.traits <- c("log_WC","log_trichome.density")

LQ.PCA.2012 <- princomp(LQ.pca.2012.df[ ,log.LQ.traits], cor = TRUE)
summary(LQ.PCA.2012) # only first component has an eigenvalue > 1
plot(LQ.PCA.2012)
LQ.PCA.2012$loadings # positive values of PC1 indicate higher levels of productivity.
LQ.PCA.2012$scores[ ,"Comp.1"]

LQ.pca.2012.df <- mutate(LQ.pca.2012.df, LQ.PC1 = LQ.PCA.2012$scores[ ,"Comp.1"]) 

hist(LQ.pca.2012.df$LQ.PC1) # much better distribution with log transformations

## Analysis

# LQ PC1 
LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ genotype*treatment +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2012.df)
plot(LQ.pc1.2012.lmer) # looks pretty good
anova(LQ.pc1.2012.lmer, ddf = "Kenward-Roger") # only genotype effect
visreg(LQ.pc1.2012.lmer, xvar = "genotype", by = "treatment")

# WC
leaf_WC.2012.lmer <-  lmer(log(leaf_WC) ~ treatment*genotype + 
                             (1|block) + (1|block:treatment), 
                           data = LQ.pca.2012.df)
plot(leaf_WC.2012.lmer) # better with log-transformation
anova(leaf_WC.2012.lmer, ddf = "Kenward-Roger") 
visreg(leaf_WC.2012.lmer, xvar = "genotype", by = "treatment")

# trichome density
trichome.density.2012.lmer <-  lmer(log(leaf_trichome.density+1) ~
                                      treatment*genotype +
                                      (1|block) + (1|block:treatment),
                                    data = LQ.pca.2012.df)
plot(trichome.density.2012.lmer) # not great, but better than no transformation
anova(trichome.density.2012.lmer, ddf = "Kenward-Roger") 
visreg(trichome.density.2012.lmer, xvar = "genotype", by = "treatment")

## Wind plant architecture 2013 ----

# 2013 
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,arch.traits]) # no outliers
scatterplotMatrix(log(filter(wind.df, Year == "2013")[ ,arch.traits])) # logging appears to help a big
corr.test(filter(wind.df, Year == "2013")[ ,arch.traits])
corr.test(log(filter(wind.df, Year == "2013")[ ,arch.traits]))

# PCA: 2013
arch.pca.2013.df <- wind.df %>%
  filter(Year == "2013") %>%
  select(block, genotype, treatment, 
         Height, all.shoot.total.length, all.shoot.count) %>%
  na.omit()

arch.PCA.2013 <- princomp(log(arch.pca.2013.df[ ,arch.traits]), cor = TRUE)
summary(arch.PCA.2013) # only first component has an eigenvalue > 1
plot(arch.PCA.2013)
arch.PCA.2013$loadings # positive values of PC1 indicate higher levels of productivity.
arch.PCA.2013$scores[ ,"Comp.1"]

arch.pca.2013.df <- mutate(arch.pca.2013.df, arch.PC1 = -1*arch.PCA.2013$scores[ ,"Comp.1"]) # multiplied by -1 so that positive values of PC1 indicate higher levels of productivity

hist(arch.pca.2013.df$arch.PC1) # normally distributed with log-transformation

## Analysis

# arch PC1 2013 
pc1.2013.lmer <- lmer(arch.PC1 ~ genotype*treatment +
                        (1|block) + (1|block:treatment),
                      data = arch.pca.2013.df)
plot(pc1.2013.lmer) # looks okay
anova(pc1.2013.lmer, ddf = "Kenward-Roger")
visreg(pc1.2013.lmer, xvar = "genotype", by = "treatment")

# plant height
height.2013 <- lmer((Height) ~ genotype*treatment + 
                      (1|block) + (1|block:treatment), 
                    data = arch.pca.2013.df)
plot(height.2013) # okay with or without log
anova(height.2013, ddf = "Kenward-Roger") # if I don't log, there is a marginally significant GxE
visreg(height.2013, xvar = "genotype", by = "treatment")


# shoot count
shoot.count.2013 <- lmer((all.shoot.count) ~ genotype*treatment +
                           (1|block) + (1|block:treatment), 
                         data = arch.pca.2013.df)
plot(shoot.count.2013) # logging doesn't really help
anova(shoot.count.2013, ddf = "Kenward-Roger")
visreg(shoot.count.2013, xvar = "genotype", by = "treatment")

# shoot length
shoot.length.2013 <- lmer(log(all.shoot.total.length) ~ genotype*treatment +
                            (1|block) + (1|block:treatment),
                          data = arch.pca.2013.df)
plot(shoot.length.2013) # much better with log-transformation
anova(shoot.length.2013, ddf = "Kenward-Roger")
visreg(shoot.length.2013, xvar = "genotype", by = "treatment")

## Wind: Leaf quality traits 2013 ----
LQ.traits.2013 <- c("leaf_WC","SLA","leaf_C_N","leaf_C","leaf_N",
                    "larva.wet.wt.exp1", "larva.wet.wt.exp2")

# explore correlations
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,LQ.traits.2013]) # non-linearity between C:N and N.
WC_nonsense <- which(filter(wind.df, Year == "2013")[ ,"leaf_WC"] < 0) # nonsensical values for leaf WC. Essentially, they suggest that the wet leaf mass is less than the dry leaf mass, which is not possible. Therefore, I'm removing them from further analysis
scatterplotMatrix(filter(wind.df, Year == "2013")[-WC_nonsense,LQ.traits.2013])
scatterplotMatrix(log(filter(wind.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]+1))
corr.test(filter(wind.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]) # larva.wet.wt.exp1 has the lowest sample size, but has a strong positive correlation with larva.wet.wt.exp2; therefore, I'm only going to use the latter for the PCA analysis.
corr.test(log(filter(wind.df, Year == "2013")[-WC_nonsense,LQ.traits.2013]+1))

# PCA
LQ.pca.2013.df <- filter(wind.df, Year == "2013")[-WC_nonsense, ] %>%
  #filter(Year == "2013") %>%
  select(block, genotype, treatment, 
         leaf_WC, SLA, leaf_C_N, larva.wet.wt.exp2) %>%
  mutate(log_CN = log(leaf_C_N),
         log_SLA = log(SLA),
         log_WC = log(leaf_WC),
         log_larva.wet = log(larva.wet.wt.exp2 + 1)) %>%
  na.omit()

pca.LQ.traits.2013 <- c("leaf_WC", "SLA", "leaf_C_N", "larva.wet.wt.exp2")
pca.log.LQ.traits.2013 <- c("log_WC", "SLA", "log_CN", "log_larva.wet")

scatterplotMatrix(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013])
corr.test(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013])

LQ.PCA.2013 <- princomp(LQ.pca.2013.df[ ,pca.log.LQ.traits.2013], cor = TRUE)
summary(LQ.PCA.2013) # only first component has an eigenvalue > 1
plot(LQ.PCA.2013)
LQ.PCA.2013$loadings # Positive values of PC1 have higher leaf water content, whereas lower values have higher SLA and C:N. Positive value of PC2 have higher Spodoptera growth rates, whereas lower values have high C:N.
LQ.PCA.2013$scores[ ,c("Comp.1","Comp.2")]

LQ.pca.2013.df <- mutate(LQ.pca.2013.df, 
                         LQ.PC1 = LQ.PCA.2013$scores[ ,"Comp.1"],
                         LQ.PC2 = LQ.PCA.2013$scores[ ,"Comp.2"]) 

hist(LQ.pca.2013.df$LQ.PC1) # normally distributed
hist(LQ.pca.2013.df$LQ.PC2) # roughly normal

## Analysis
# note that I use the filtered wind.df for the individual traits to preserve the original sample sizes.

# LQ PC1
LQ.pc1.2013.lmer <- lmer(LQ.PC1 ~ genotype*treatment +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2013.df)
plot(LQ.pc1.2013.lmer) # looks good
anova(LQ.pc1.2013.lmer, ddf = "Kenward-Roger") # only genotype effect
visreg(LQ.pc1.2013.lmer, xvar = "genotype", by = "treatment")

# LQ PC2
LQ.pc2.2013.lmer <- lmer(LQ.PC2 ~ genotype*treatment +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2013.df)
plot(LQ.pc2.2013.lmer) # looks okay
anova(LQ.pc2.2013.lmer, ddf = "Kenward-Roger") # only genotype effect
visreg(LQ.pc2.2013.lmer, xvar = "genotype", by = "treatment")

# leaf C:N
leaf_CN.2013.lmer <-  lmer(log(leaf_C_N) ~ treatment*genotype + 
                             (1|block) + (1|block:treatment), 
                           data = filter(wind.df, Year == "2013"))
plot(leaf_CN.2013.lmer) # residuals a bit better with log transformation
anova(leaf_CN.2013.lmer, ddf = "Kenward-Roger") # qualitatively same results with or without log transformation

# SLA
SLA.2013.lmer <-  lmer(SLA ~ treatment*genotype + 
                         (1|block) + (1|block:treatment), 
                       data = filter(wind.df, Year == "2013"))
plot(SLA.2013.lmer) # okay
anova(SLA.2013.lmer, ddf = "Kenward-Roger") 

# WC
leaf_WC.2013.lmer <-  lmer(log(leaf_WC) ~ treatment*genotype + 
                             (1|block) + (1|block:treatment), 
                           data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
plot(leaf_WC.2013.lmer) # better with log-transformation
anova(leaf_WC.2013.lmer, ddf = "Kenward-Roger") 
visreg(leaf_WC.2013.lmer, xvar = "genotype", by = "treatment")

# larva wet weight exp. 2
larva.wet.wt.exp2.lmer <-  lmer(log(larva.wet.wt.exp2+1) ~ treatment*genotype +
                                  (1|block) + (1|block:treatment), 
                                data = filter(wind.df, Year == "2013"))
plot(larva.wet.wt.exp2.lmer) # better with log transformation
anova(larva.wet.wt.exp2.lmer, ddf = "Kenward-Roger") 
visreg(larva.wet.wt.exp2.lmer, xvar = "genotype", by = "treatment")

## Ant-aphid 2012: architecture ----
aa.arch.traits <- c("Height","all.shoot.count","mature.shoot.total.length")

# explore correlations
scatterplotMatrix(filter(aa.df, Year == "2012")[ ,aa.arch.traits])

# PCA: 2012
aa.arch.pca.2012.df <- filter(aa.df, Year == "2012") %>%
  #filter(Year == "2012") %>%
  select(Block, Genotype, Aphid.Treatment, Ant.Mound.Dist, fact.Ant.Mound.Dist,
         Height, mature.shoot.total.length, all.shoot.count) %>%
  na.omit()

aa.arch.PCA.2012 <- princomp(aa.arch.pca.2012.df[ ,aa.arch.traits], cor = TRUE)
summary(aa.arch.PCA.2012) # only first component has an eigenvalue > 1
plot(aa.arch.PCA.2012)
aa.arch.PCA.2012$loadings 
aa.arch.PCA.2012$scores[ ,"Comp.1"]

aa.arch.pca.2012.df <- mutate(aa.arch.pca.2012.df, arch.PC1 = -1*aa.arch.PCA.2012$scores[ ,"Comp.1"]) # multiplied by -1 so that positive values of PC1 indicate higher levels of productivity.

hist(aa.arch.pca.2012.df$arch.PC1) # roughly normal distribution

## Analysis

# arch PC1 2012 
aa.pc1.2012.lmer <- lmer(arch.PC1 ~ Genotype*Aphid.Treatment*Ant.Mound.Dist +
                           (1|Block) + (1|Block:fact.Ant.Mound.Dist),
                      data = aa.arch.pca.2012.df)
plot(aa.pc1.2012.lmer) # looks pretty good
anova(aa.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.pc1.2012.lmer, xvar = "Genotype")

# plant height
height.2012 <- lmer(Height ~ Genotype*Aphid.Treatment*Ant.Mound.Dist + 
                      (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                    data = filter(aa.df, Year == "2012")) 
plot(height.2012) # not great
anova(height.2012, ddf = "Kenward-Roger")
visreg(height.2012, xvar = "Genotype")

# shoot count
aa.shoot.count.2012 <- lmer(log(all.shoot.count) ~ Genotype*Aphid.Treatment*Ant.Mound.Dist + 
                      (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                    data = filter(aa.df, Year == "2012")) 
plot(aa.shoot.count.2012) # better with log transformation
anova(aa.shoot.count.2012, ddf = "Kenward-Roger")
visreg(aa.shoot.count.2012, xvar = "Genotype")

# shoot length. 
aa.mature.length.2012 <- lmer(mature.shoot.total.length ~ Genotype*Aphid.Treatment*Ant.Mound.Dist + 
                      (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                    data = filter(aa.df, Year == "2012")) 
plot(aa.mature.length.2012) # not great
anova(aa.mature.length.2012, ddf = "Kenward-Roger")
visreg(aa.mature.length.2012, xvar = "Genotype")

## Ant-aphid 2012: leaf quality

# explore correlations
scatterplotMatrix(filter(aa.df, Year == "2012")[ ,LQ.traits])
corr.test(filter(aa.df, Year == "2012")[ ,LQ.traits])

# PCA
aa.LQ.pca.2012.df <- filter(aa.df, Year == "2012") %>%
  select(Block, Genotype, Aphid.Treatment, 
         Ant.Mound.Dist, fact.Ant.Mound.Dist, 
         leaf_WC, leaf_trichome.density) %>%
  mutate(log_WC = log(leaf_WC),
         log_trichomes = log(leaf_trichome.density + 1)) %>%
  na.omit()

aa.pca.log.LQ.traits <- c("log_WC", "log_trichomes")

scatterplotMatrix(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])
corr.test(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])

aa.LQ.PCA <- princomp(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits], cor = TRUE)
summary(aa.LQ.PCA) # only first component has an eigenvalue > 1
plot(aa.LQ.PCA)
aa.LQ.PCA$loadings # Positive values of PC1 have higher leaf water content, whereas lower values have higher SLA and C:N. Positive value of PC2 have higher Spodoptera growth rates, whereas lower values have high C:N.
aa.LQ.PCA$scores[ ,"Comp.1"]

aa.LQ.PCA.df <- mutate(aa.LQ.pca.2012.df, 
                         LQ.PC1 = -1*aa.LQ.PCA$scores[ ,"Comp.1"]) # multiply by -1 so positive values indicate higher leaf water content and higher trichome density

hist(aa.LQ.PCA.df$LQ.PC1) # normally distributed (skewed if I don't use the logged traits)

## Analysis

# LQ PC1
aa.LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ Genotype + Aphid.Treatment + Ant.Mound.Dist + 
                              (1|Block) + (1|Block:fact.Ant.Mound.Dist),
                         data = aa.LQ.PCA.df)
plot(aa.LQ.pc1.2012.lmer) 
anova(aa.LQ.pc1.2012.lmer, ddf = "Kenward-Roger") 
visreg(aa.LQ.pc1.2012.lmer, xvar = "Genotype")

## Analysis

# WC
aa.WC.2012 <- lmer(log(leaf_WC) ~ Genotype + Aphid.Treatment + Ant.Mound.Dist + 
                      (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                    data = filter(aa.df, Year == "2012")) 
plot(aa.WC.2012) # not great
anova(aa.WC.2012, ddf = "Kenward-Roger")
visreg(aa.WC.2012, xvar = "Genotype")

# trichome density. Used reduced model since we didn't have a completely factorial design.
aa.trichomes.2012 <- lmer(log(leaf_trichome.density+1) ~ Genotype + Aphid.Treatment + Ant.Mound.Dist + 
                              (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                            data = filter(aa.df, Year == "2012")) 
plot(aa.trichomes.2012) # better with log transformation, but still weird
anova(aa.trichomes.2012, ddf = "Kenward-Roger")
visreg(aa.trichomes.2012, xvar = "Genotype")

## old ----
arch.df.2013 <- na.omit(filter(wind.df, Year == "2013")[ ,c("genotype","treatment",arch.traits.red)])
rda.arch.2013 <- rda(arch.df.2013[ ,-c(1,2)] ~ genotype*treatment, data = arch.df.2013, scale = TRUE)
summary(rda.arch.2013)
plot(rda.arch.2013, display = c("sp","cn"))

## trait correlations
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,9:17])
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,c(9,12,14,17)])
scatterplotMatrix(log(filter(wind.df, Year == "2013")[ ,c(9,12,14,17)]))
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,c(18:20)])
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,c(19,21:25,27)])
corr.test(filter(wind.df, Year == "2013")[ ,c(19,21:25,27)])
wind.df[which(wind.df$Height > 80),] # not a typo according to data entry, but definitely seems to be an outlier...Maybe it was written down incorrectly in the field. It is almost twice the size of the other plants...
corr.test(filter(wind.df, Year == "2012")[ ,9:17])
corr.test(filter(wind.df, Year == "2013")[ ,9:17])

# 2013
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,arch.traits])
scatterplotMatrix(log(filter(wind.df, Year == "2013")[ ,arch.traits]))

corr.test(log(filter(wind.df, Year == "2013")[,arch.traits]))


## PCA on plant architecture traits
library(vegan)
arch.pca.df <- na.omit(filter(wind.df, Year == "2012", Height < 80)[ ,c(9,12,14,17)])#[ ,c(9,10,11,15,16)])
arch.pca <- princomp(arch.pca.df, cor = TRUE)
arch.rda <- rda(arch.pca.df ~ genotype, scale = TRUE, filter(wind.df, Year == "2012", Height < 80))
plot(arch.rda, display = c("sp","cn"))
summary(arch.pca)
loadings(arch.pca)
varimax(loadings(arch.pca))

plot(all.shoot.total.length ~ Height, filter(wind.df, Year == "2012", Height < 80)[ ,9:17])

# calculate means for each year, treatment, and genotype combination for graphing.
trait.summary <- wind.df %>%
  group_by(Year, treatment, genotype) %>%
  summarise_each(funs(mean_na.rm = mean(., na.rm = TRUE))) %>%
  select(-(X:plant.id)) %>%
  melt(id.vars = c("Year","treatment","genotype"))

ggplot(filter(trait.summary, variable %in% c("Height","all.shoot.total.length", "all.shoot.count")), 
       aes(x = treatment, y = value, color = genotype, group = genotype)) +
  geom_line() +
  facet_grid(variable ~ Year, scales = "free")

## Notes for these models. 

# plant survival. with a fully specified fixed effect model, it isn't running. I'm also worried about the correlations between random effects if I try to include (treatment|genotype) in the model
with(filter(wind.df, Year == "2012", Dead == 0), table(treatment, genotype))
with(filter(wind.df, Year == "2013", Dead == 0), table(treatment, genotype))

dead.2012.glmer <- lmer(Dead ~ treatment*genotype + (1|block/treatment),
                    data = filter(wind.df, Year == "2012"))#, family = "binomial")
summary(dead.2012.glmer)
anova(dead.2012.glmer, ddf = "Kenward-Roger")

dead.2013.glmer <- glmer(Dead ~ treatment + (1|genotype) + (1|block/treatment),
                         data = filter(wind.df, Year == "2013"), 
                         family = "binomial")
summary(dead.2013.glmer)
confint(profile(dead.2013.glmer))

# plant height
ggplot(wind.df, aes(x = genotype, y = Height, color = treatment)) +
  geom_boxplot() +
  facet_wrap(~Year, ncol = 2)

height.2012.lmer <- lmer(log(Height) ~ treatment*genotype + (1|block/treatment), 
                         data = filter(wind.df, Year == "2012"))
summary(height.2012.lmer)
plot(height.2012.lmer)
anova(height.2012.lmer, ddf = "Kenward-Roger") # get same qualitative answer whether or not I include the outlying Height data point (88 cm)

height.2013.lmer <- lmer((Height) ~ (1|treatment) + (1|genotype) + (1|Year) + (1|block/treatment), 
                         data = filter(wind.df, Height < 80))
summary(height.2013.lmer)
30.345/(30.345+65.572+66.945+16.697+3.494+7.129)
16.697/(30.345+65.572+66.945+16.697+3.494+7.129)
confint(profile(height.2013.lmer))
AIC(height.2013.lmer)
plot(height.2013.lmer)
anova(height.2013.lmer)
anova(height.2013.lmer, ddf = "Kenward-Roger")

# all.shoot.total.length
ggplot(wind.df, aes(x = genotype, y = all.shoot.total.length, color = treatment)) +
  geom_boxplot() +
  facet_wrap(~Year, ncol = 2)

all.shoot.total.length.2012.lmer <- lmer(all.shoot.total.length ~ treatment*genotype + (1|block/treatment), 
                         data = filter(wind.df, Year == "2012"))
summary(all.shoot.total.length.2012.lmer)
plot(all.shoot.total.length.2012.lmer)
anova(all.shoot.total.length.2012.lmer, ddf = "Kenward-Roger")

all.shoot.total.length.2013.lmer <- lmer(log(all.shoot.total.length) ~ treatment*genotype + (1|block/treatment), 
                         data = filter(wind.df, Year == "2013"))
summary(all.shoot.total.length.2013.lmer)
plot(all.shoot.total.length.2013.lmer) # residual suggest log transformation is appropriate
anova(all.shoot.total.length.2013.lmer, ddf = "Kenward-Roger")

# shoot:height ratio. Don't know if this is capturing what I was hoping
ggplot(wind.df, aes(x = genotype, y = all.shoot.total.length/Height, color = treatment)) +
  geom_boxplot() +
  facet_wrap(~Year, ncol = 2)

shoot.height.2012.lmer <- lmer(all.shoot.total.length/Height ~ treatment*genotype + (1|block/treatment), 
                                         data = filter(wind.df, Year == "2012"))
summary(shoot.height.2012.lmer)
plot(shoot.height.2012.lmer)
anova(shoot.height.2012.lmer, ddf = "Kenward-Roger")

shoot.height.2013.lmer <-  lmer(log(all.shoot.total.length/Height) ~ treatment*genotype + (1|block/treatment), 
                               data = filter(wind.df, Year == "2013"))
summary(shoot.height.2013.lmer)
plot(shoot.height.2013.lmer)
anova(shoot.height.2013.lmer, ddf = "Kenward-Roger")

# all.shoot.count
ggplot(wind.df, aes(x = genotype, y = all.shoot.count, color = treatment)) +
  geom_boxplot() +
  facet_wrap(~Year, ncol = 2)

all.shoot.count.2012.lmer <- lmer(all.shoot.count ~ treatment*genotype + (1|block), 
                               data = filter(wind.df, Year == "2012"))
summary(all.shoot.count.2012.lmer)
plot(all.shoot.count.2012.lmer)
anova(all.shoot.count.2012.lmer, ddf = "Kenward-Roger")

all.shoot.count.2013.lmer <-  lmer(all.shoot.count ~ treatment*genotype + (1|block), 
                                data = filter(wind.df, Year == "2013"))
summary(all.shoot.count.2013.lmer)
plot(all.shoot.count.2013.lmer)
anova(all.shoot.count.2013.lmer, ddf = "Kenward-Roger")


