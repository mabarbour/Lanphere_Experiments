## load required libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
library(lmerTest)
library(RLRsim)
library(car)

## load required data set
wind.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_trait_df.csv') %>%
  tbl_df() %>%
  mutate(block = as.factor(block),
         Year = as.factor(Year),
         plant.position = as.character(plant.position))
glimpse(wind.df)

## trait correlations
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,9:17])
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,c(18:20)])
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,c(19,21:25,27)])
corr.test(filter(wind.df, Year == "2013")[ ,c(19,21:25,27)])
wind.df[which(wind.df$Height > 80),] # not a typo according to data entry, but definitely seems to be an outlier...Maybe it was written down incorrectly in the field. It is almost twice the size of the other plants...
corr.test(filter(wind.df, Year == "2012")[ ,9:17])
corr.test(filter(wind.df, Year == "2013")[ ,9:17])

## PCA on plant architecture traits
library(vegan)
arch.pca.df <- na.omit(filter(wind.df, Year == "2012", Height < 80)[ ,c(9,10,11,15,16)])
arch.pca <- princomp(arch.pca.df, cor = TRUE)
arch.rda <- rda(arch.pca.df, scale = TRUE)
plot(arch.rda, display = 'sp')
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

dead.2013.glmer <- glmer(Dead ~ treatment + (1|genotype) + (1|block),
                         data = filter(wind.df, Year == "2013"), family = "binomial")
summary(dead.2013.glmer)

# plant height
height.2012.lmer <- lmer(Height ~ treatment + (1|genotype) + (0 + treatment|genotype) + (1|block/treatment), 
                         data = filter(wind.df, Year == "2012"))
summary(height.2012.lmer)
plot(height.2012.lmer)
anova(height.2012.lmer, ddf = "Kenward-Roger") # get same qualitative answer whether or not I include the outlying Height data point (88 cm)

height.2013.lmer <- lmer(Height ~ treatment*genotype + (1|block/treatment), 
                         data = filter(wind.df, Year == "2013"))
summary(height.2013.lmer)
plot(height.2013.lmer)
anova(height.2013.lmer)
anova(height.2013.lmer, ddf = "Kenward-Roger")

# all.shoot.total.length
all.shoot.total.length.2012.lmer <- lmer(all.shoot.total.length ~ treatment*genotype + (1|block), 
                         data = filter(wind.df, Year == "2012"))
summary(all.shoot.total.length.2012.lmer)
plot(all.shoot.total.length.2012.lmer)
anova(all.shoot.total.length.2012.lmer, ddf = "Kenward-Roger")

all.shoot.total.length.2013.lmer <- lmer(log(all.shoot.total.length) ~ treatment*genotype + (1|block), 
                         data = filter(wind.df, Year == "2013"))
summary(all.shoot.total.length.2013.lmer)
plot(all.shoot.total.length.2013.lmer) # residual suggest log transformation is appropriate
anova(all.shoot.total.length.2013.lmer, ddf = "Kenward-Roger")

# shoot:height ratio
shoot.height.2012.lmer <- lmer(all.shoot.total.length/Height ~ treatment*genotype + (1|block), 
                                         data = filter(wind.df, Year == "2012"))
summary(shoot.height.2012.lmer)
plot(shoot.height.2012.lmer)
anova(shoot.height.2012.lmer, ddf = "Kenward-Roger")

shoot.height.2013.lmer <-  lmer(log(all.shoot.total.length/Height) ~ treatment*genotype + (1|block), 
                               data = filter(wind.df, Year == "2013"))
summary(shoot.height.2013.lmer)
plot(shoot.height.2013.lmer)
anova(shoot.height.2013.lmer, ddf = "Kenward-Roger")

# all.shoot.count
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

# leaf C:N
leaf_CN.2013.lmer <-  lmer(leaf_C_N ~ treatment*genotype + (1|block), 
                                data = filter(wind.df, Year == "2013"))
summary(leaf_CN.2013.lmer)
plot(leaf_CN.2013.lmer)
anova(leaf_CN.2013.lmer, ddf = "Kenward-Roger") # qualitatively same results with or without log transformation

# SLA
SLA.2013.lmer <-  lmer(SLA ~ treatment*genotype + (1|block), 
                           data = filter(wind.df, Year == "2013"))
summary(SLA.2013.lmer)
plot(SLA.2013.lmer)
anova(SLA.2013.lmer, ddf = "Kenward-Roger") 

# leaf_WC
leaf_WC.2012.lmer <-  lmer(log(leaf_WC) ~ treatment*genotype + (1|block), 
                           data = filter(wind.df, Year == "2012"))
summary(leaf_WC.2012.lmer)
plot(leaf_WC.2012.lmer)
anova(leaf_WC.2012.lmer, ddf = "Kenward-Roger") 

leaf_WC.2013.lmer <-  lmer(leaf_WC ~ treatment*genotype + (1|block), 
                       data = filter(wind.df, Year == "2013"))
summary(leaf_WC.2013.lmer)
plot(leaf_WC.2013.lmer)
anova(leaf_WC.2013.lmer, ddf = "Kenward-Roger") 

# leaf percent.browned

# leaf_trichome.density
leaf_trichome.density.2012.lmer <-  lmer(log(leaf_trichome.density+1) ~ treatment*genotype + (1|block), 
                           data = filter(wind.df, Year == "2012"))
summary(leaf_trichome.density.2012.lmer)
plot(leaf_trichome.density.2012.lmer)
anova(leaf_trichome.density.2012.lmer, ddf = "Kenward-Roger") 

# larva wet weight exp. 1
with(filter(wind.df, Year == "2013", larva.survival.exp1 == 1), table(treatment, genotype)) # lower sample sizes in first vs. second experiment.
larva.wet.wt.exp1.lmer <-  lmer(larva.wet.wt.exp1 ~ treatment*genotype + (1|block), 
                                         data = filter(wind.df, Year == "2013"))
summary(larva.wet.wt.exp1.lmer)
plot(larva.wet.wt.exp1.lmer)
anova(larva.wet.wt.exp1.lmer, ddf = "Kenward-Roger") 

# larva wet weight exp. 2
with(filter(wind.df, Year == "2013", larva.survival.exp2 == 1), table(treatment, genotype))
larva.wet.wt.exp2.lmer <-  lmer(log(larva.wet.wt.exp2+1) ~ treatment*genotype + (1|block), 
                                data = filter(wind.df, Year == "2013"))
summary(larva.wet.wt.exp2.lmer)
plot(larva.wet.wt.exp2.lmer)
anova(larva.wet.wt.exp2.lmer, ddf = "Kenward-Roger") 
