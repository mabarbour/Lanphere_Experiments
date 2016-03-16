
## load required libraries ----
library(dplyr)
library(reshape2)
library(ggplot2)
library(visreg)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
library(lmerTest)
library(piecewiseSEM) # for calculating r2 of mixed-effect models.
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
arch.traits <- c("Height","all.shoot.avg.length","all.shoot.count")

# 2012 
scatterplotMatrix(filter(wind.df, Year == "2012")[ ,arch.traits])
height.outlier <- which(wind.df$Height > 80)  # not a typo according to data entry, but definitely seems to be an outlier. Maybe it was written down incorrectly in the field, because it is almost twice the size of the other plants. I've decided to remove it from this analysis.
scatterplotMatrix(filter(wind.df, Year == "2012")[-height.outlier,
                                                  arch.traits])
scatterplotMatrix(log(filter(wind.df, Year == "2012")[-height.outlier,
                                                      arch.traits])) # logging doesn't appear to linearize these relationships anymore.
corr.test(filter(wind.df, Year == "2012")[-height.outlier,arch.traits]) # strong correlation between total shoot length and height in 2012

# PCA: 2012
arch.pca.2012.df <- filter(wind.df, Year == "2012")[-height.outlier, ] %>%
  select(block, genotype, treatment, 
         Height, all.shoot.avg.length, all.shoot.count) %>%
  na.omit()

arch.PCA.2012 <- princomp(arch.pca.2012.df[ ,arch.traits], cor = TRUE)
summary(arch.PCA.2012) # only first component has an eigenvalue > 1
plot(arch.PCA.2012)
arch.PCA.2012$loadings # positive values of PC1 indicate higher levels of plant biomass.
arch.PCA.2012$scores[ ,"Comp.1"]

arch.pca.2012.df <- mutate(arch.pca.2012.df, arch.PC1 = arch.PCA.2012$scores[ ,"Comp.1"], arch.PC2=arch.PCA.2012$scores[ ,"Comp.2"]) 

hist(arch.pca.2012.df$arch.PC1) # normally distributed
hist(arch.pca.2012.df$arch.PC2)

## Analysis

## arch PC1 2012 
with(arch.pca.2012.df, interaction.plot(treatment, genotype, arch.PC1)) # looks like main effect of wind and genotype

pc1.2012.lmer <- lmer(arch.PC1 ~ treatment*genotype +
                        (1|block) + (1|block:treatment),
                      data = arch.pca.2012.df,
                      contrasts = list(treatment = "contr.sum",
                                       genotype = "contr.sum"))
summary(pc1.2012.lmer)
plot(pc1.2012.lmer) # looks pretty good

#pc1.2012.mains <- update(pc1.2012.lmer, .~. -treatment:genotype)
#pc1.2012.geno <- 
# test for GxE
anova(pc1.2012.lmer,
      update(pc1.2012.lmer, .~. -treatment*genotype) + genotype),
      update(pc1.2012.lmer, .~. -(treatment|genotype)))

anova(pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(pc1.2012.lmer, xvar = "treatment")

# random model
pc1.2012.rand <- lmer(arch.PC1 ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2012.df)
pc1.2012.vars <- as.data.frame(VarCorr(pc1.2012.rand))

# plant height
with(filter(wind.df, Height > 0), interaction.plot(treatment, genotype, Height))

height.lmer <- lmer(Height ~ treatment*genotype*Year +
                        (1|block) + 
                        (1|block:treatment) +
                        (1|plant.id),
                    data = wind.df,
                    contrasts = list(treatment = "contr.sum",
                                     genotype = "contr.sum",
                                     Year = "contr.sum"))
summary(height.lmer)

#anova(height.lmer, update(height.lmer, .~. -(treatment|Year) + (1|Year)))
summary(update(height.2012.lmer, .~. -1))
31.225/26.237 # willows were 1.2x taller in unexposed plots
1 - (26.237/31.225) # willows were 16% shorter in unexposed plots

plot(height.2012.lmer) # looks pretty good


# test for GxE
anova(height.2012.lmer,
      update(height.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(height.2012.lmer, .~. -(treatment|genotype)))

anova(height.lmer, ddf = "Kenward-Roger")
sem.model.fits(height.lmer)
visreg(height.2012.lmer, xvar = "treatment")

# random model
height.2012.rand <- lmer(Height ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2012.df)
summary(height.2012.rand)
20.505/(2.684+0+20.505+11.515+43.435) # 26% of the variance explained by genotype
11.515/(2.684+0+20.505+11.515+43.435) # 15% of the variance explained by wind exposure
height.2012.vars <- as.data.frame(VarCorr(height.2012.rand))

# shoot count
with(filter(wind.df, all.shoot.count > 0), interaction.plot(treatment, Year, all.shoot.count))
with(filter(wind.df, all.shoot.count > 0), interaction.plot(Year, genotype, all.shoot.count))

shoot.count.glmer <- lmer(all.shoot.count ~ treatment*genotype*Year +
                           (1|block) + 
                           (1|block:treatment) +
                           (1|plant.id),
                         data = wind.df,
                         contrasts = list(treatment = "contr.sum",
                                          genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(shoot.count.glmer)
overdisp_fun(shoot.count.glmer) # no overdispersion
plot(shoot.count.glmer) # looks pretty good

visreg(shoot.count.glmer, 
       xvar = "treatment", by = "genotype",
       scale = "response")
visreg(shoot.count.glmer, 
       xvar = "treatment", by = "Year",
       scale = "response")
visreg(shoot.count.glmer, 
       xvar = "Year", by = "genotype",
       scale = "response")

# Likelihood ratio tests
drop1(shoot.count.glmer, test = "Chisq")
drop1(update(shoot.count.glmer, .~. -treatment:genotype:Year), test = "Chisq")
drop1(update(shoot.count.glmer, .~. -treatment*genotype*Year + treatment + genotype + Year), test = "Chisq")

## Calculate R2 for significant predictors

# interaction effects
shoot.count.2x <- update(shoot.count.glmer, .~. -treatment:genotype:Year)
shoot.count.GenoYear <- update(shoot.count.2x, .~. -genotype:Year)
shoot.count.WindYear <- update(shoot.count.2x, .~. -treatment:Year)

shoot.count.2x.fits <- sem.model.fits(list(shoot.count.2x, shoot.count.GenoYear, shoot.count.WindYear))
shoot.count.2x.fits[1, "Marginal"] - shoot.count.2x.fits[2, "Marginal"] # 2.6% of variance explained by GxYear
shoot.count.2x.fits[1, "Marginal"] - shoot.count.2x.fits[3, "Marginal"] # 1.1% of variance explained by wind x year.

# main effects
shoot.count.main <- update(shoot.count.glmer, .~. -treatment*genotype*Year + treatment + genotype + Year)
shoot.count.geno <- update(shoot.count.main, .~. -genotype)
shoot.count.wind <- update(shoot.count.main, .~. -treatment)
shoot.count.year <- update(shoot.count.main, .~. -Year)

shoot.count.main.fits <- sem.model.fits(list(shoot.count.main, shoot.count.geno, shoot.count.wind, shoot.count.year))
shoot.count.main.fits[1, "Marginal"] - shoot.count.main.fits[2, "Marginal"] # 13.1% of variance explained by genotype
shoot.count.main.fits[1, "Marginal"] - shoot.count.main.fits[3, "Marginal"] # 7.4% of variance explained by wind
shoot.count.main.fits[1, "Marginal"] - shoot.count.main.fits[4, "Marginal"] # <1% of variance explained by year

# average shoot length
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(treatment, genotype, all.shoot.avg.length))
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(treatment, Year, all.shoot.avg.length))
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(Year, genotype, all.shoot.avg.length))
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ Year))
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ genotype))
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ treatment))

shoot.length.lmer <- lmer(log(all.shoot.avg.length) ~ treatment*genotype*Year +
                            (1|block) + 
                            (1|block:treatment) +
                            (1|plant.id),
                          data = wind.df,
                          contrasts = list(treatment = "contr.sum",
                                           genotype = "contr.sum",
                                           Year = "contr.sum"))
summary(shoot.length.lmer)
plot(shoot.length.lmer) # looks okay

visreg(shoot.length.lmer, 
       xvar = "treatment", by = "genotype",
       scale = "response")
visreg(shoot.length.lmer, 
       xvar = "treatment", by = "Year",
       scale = "response")
visreg(shoot.length.lmer, 
       xvar = "Year", by = "genotype",
       scale = "response")

# Kenward-Roger tests
anova(shoot.length.lmer, ddf = "Kenward-Roger")

## Calculate R2 for significant predictors

# no interaction effects were significant so I didn't calculate R2

# main effects
shoot.length.main <- update(shoot.length.lmer, .~. -treatment*genotype*Year + treatment + genotype + Year)
shoot.length.geno <- update(shoot.length.main, .~. -genotype)
shoot.length.wind <- update(shoot.length.main, .~. -treatment)
shoot.length.year <- update(shoot.length.main, .~. -Year)

shoot.length.main.fits <- sem.model.fits(list(shoot.length.main, shoot.length.geno, shoot.length.wind, shoot.length.year))
shoot.length.main.fits[1, "Marginal"] - shoot.length.main.fits[2, "Marginal"] # 7.5% of variance explained by genotype
shoot.length.main.fits[1, "Marginal"] - shoot.length.main.fits[3, "Marginal"] # 3.8% of variance explained by wind
shoot.length.main.fits[1, "Marginal"] - shoot.length.main.fits[4, "Marginal"] # 13% of variance explained by year

# shoot length
shoot.length.2012.lmer <- lmer(all.shoot.avg.length ~ treatment + (treatment|genotype) +
                                (1|block) + (1|block:treatment),
                              data = filter(wind.df, Year == "2012"))
summary(shoot.length.2012.lmer)
plot(shoot.length.2012.lmer) # looks pretty good

# test for GxE
anova(shoot.length.2012.lmer,
      update(shoot.length.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(shoot.length.2012.lmer, .~. -(treatment|genotype)))

anova(shoot.length.2012.lmer, ddf = "Kenward-Roger")
visreg(shoot.length.2012.lmer, xvar = "treatment")

# random model
shoot.length.2012.rand <- lmer(all.shoot.total.length ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2012.df)
summary(shoot.length.2012.rand)
1.811/(2.761+1.811+7.731+9.7+28.827) # 4% of the variance explained by genotype
7.731/(2.761+1.811+7.731+9.7+28.827) # 15% of the variance explained by GxE
9.7/(2.761+1.811+7.731+9.7+28.827) # 19% of the variance explained by wind exposure
shoot.length.2012.vars <- as.data.frame(VarCorr(shoot.length.2012.rand))

mean(c(29,26,15)) # 23% mean effect of G and/or GxE in 2012
mean(c(19,4,15)) # 13% mean effect of wind

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
corr.test(filter(wind.df, Year == "2012")[-WC.outlier,LQ.traits]) # positive correlation between leaf_WC and trichome density

# PCA: 2012
LQ.pca.2012.df <- filter(wind.df, Year == "2012")[-WC.outlier, ] %>%
  select(block, genotype, treatment, 
         leaf_WC, leaf_trichome.density) %>%
  mutate(log_trichome.density = log(leaf_trichome.density+1),
         log_WC = log(leaf_WC)) %>%
  na.omit()

log.LQ.traits <- c("log_WC","log_trichome.density")

LQ.PCA.2012 <- princomp(LQ.pca.2012.df[ ,log.LQ.traits], cor = TRUE)
summary(LQ.PCA.2012) # only first component has an eigenvalue > 1
plot(LQ.PCA.2012)
LQ.PCA.2012$loadings # positive values of PC1 indicate higher levels of leaf water content and trichome densities.
LQ.PCA.2012$scores[ ,"Comp.1"]

LQ.pca.2012.df <- mutate(LQ.pca.2012.df, LQ.PC1 = LQ.PCA.2012$scores[ ,"Comp.1"]) 

hist(LQ.pca.2012.df$LQ.PC1) # much better distribution with log transformations

## Analysis

# LQ PC1 
LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ treatment + (treatment|genotype) +
                        (1|block) + (1|block:treatment),
                      data = LQ.pca.2012.df)
summary(LQ.pc1.2012.lmer)
plot(LQ.pc1.2012.lmer) # looks okay

# test for GxE
anova(LQ.pc1.2012.lmer,
      update(LQ.pc1.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(LQ.pc1.2012.lmer, .~. -(treatment|genotype)))

anova(LQ.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc1.2012.lmer, xvar = "treatment")

# random model
LQ.pc1.2012.rand <- lmer(LQ.PC1 ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = LQ.pca.2012.df)
summary(LQ.pc1.2012.rand )
LQ.pc1.2012.vars <- as.data.frame(VarCorr(LQ.pc1.2012.rand))

# WC
leaf_WC.2012.lmer <- lmer(log(leaf_WC) ~ treatment + (treatment|genotype) +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2012.df)
summary(leaf_WC.2012.lmer)
plot(leaf_WC.2012.lmer) # looks okay

# test for GxE
anova(leaf_WC.2012.lmer,
      update(leaf_WC.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(leaf_WC.2012.lmer, .~. -(treatment|genotype)))

anova(leaf_WC.2012.lmer, ddf = "Kenward-Roger")
visreg(leaf_WC.2012.lmer, xvar = "treatment")

# random model
leaf_WC.2012.rand <- lmer(log(leaf_WC) ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = LQ.pca.2012.df)
summary(leaf_WC.2012.rand )
leaf_WC.2012.vars <- as.data.frame(VarCorr(leaf_WC.2012.rand)) %>%
  select(group = grp, var = vcov, sd = sdcor) %>%
  mutate(response = "log.leaf_WC",
         perc.var = round(var/sum(var),2)) %>%
  tbl_df()

## function to calculate variance components
#library(varComp)
#leaf_WC.2012.rand <- varComp(fixed = log(leaf_WC) ~ 1,
 #                            data = LQ.pca.2012.df,
  #                           random = ~ treatment + genotype + block + block:treatment)
#summary(leaf_WC.2012.rand)
#varComp.test(leaf_WC.2012.rand)

# trichome density
leaf_trichome.density.2012.lmer <- lmer(log(leaf_trichome.density+1) ~ treatment + (treatment|genotype) +
                            (1|block) + (1|block:treatment),
                          data = LQ.pca.2012.df)
summary(leaf_trichome.density.2012.lmer)
plot(leaf_trichome.density.2012.lmer) # not great, could consider glmer poisson

# test for GxE
anova(leaf_trichome.density.2012.lmer,
      update(leaf_trichome.density.2012.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(leaf_trichome.density.2012.lmer, .~. -(treatment|genotype)))

anova(leaf_trichome.density.2012.lmer, ddf = "Kenward-Roger")
visreg(leaf_trichome.density.2012.lmer, xvar = "treatment")

# random model
leaf_trichome.density.2012.rand <- lmer(log(leaf_trichome.density+1) ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = LQ.pca.2012.df)
summary(leaf_trichome.density.2012.rand )
leaf_trichome.density.2012.vars <- as.data.frame(VarCorr(leaf_trichome.density.2012.rand)) %>%
  select(group = grp, var = vcov, sd = sdcor) %>%
  mutate(response = "log+1.trichome_density",
         perc.var = round(var/sum(var),2)) %>%
  tbl_df()

## Wind plant architecture 2013 ----

# 2013 
scatterplotMatrix(filter(wind.df, Year == "2013")[ ,arch.traits]) # no outliers
scatterplotMatrix(log(filter(wind.df, Year == "2013")[ ,arch.traits])) # logging appears to help a bit
corr.test(filter(wind.df, Year == "2013")[ ,arch.traits]) # all traits are now highly correlated
corr.test(log(filter(wind.df, Year == "2013")[ ,arch.traits]))

# PCA: 2013
arch.pca.2013.df <- wind.df %>%
  filter(Year == "2013") %>%
  select(block, genotype, treatment, 
         Height, all.shoot.avg.length, all.shoot.count) %>%
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
with(arch.pca.2013.df, interaction.plot(treatment, genotype, arch.PC1))

pc1.2013.lmer <- lmer(arch.PC1 ~ treatment + (treatment|genotype) +
                        (1|block) + (1|block:treatment),
                      data = arch.pca.2013.df)
summary(pc1.2013.lmer)
plot(pc1.2013.lmer) # looks pretty good

# test for GxE
anova(pc1.2013.lmer,
      update(pc1.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(pc1.2013.lmer, .~. -(treatment|genotype)))

anova(pc1.2013.lmer, ddf = "Kenward-Roger")
visreg(pc1.2013.lmer, xvar = "treatment")

# random model
pc1.2013.rand <- lmer(arch.PC1 ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2013.df)
summary(pc1.2013.rand)
0.5305/(0.3961+0.5186+0.2041+0.5305+1.0331) # 20% of the variance for wind exposure
0.2041/(0.3961+0.5186+0.2041+0.5305+1.0331) # 8% of the variance for genotype.
pc1.2013.vars <- as.data.frame(VarCorr(pc1.2013.rand))

# plant height
with(arch.pca.2013.df, interaction.plot(treatment, genotype, Height))
height.2013.lmer <- lmer(Height ~ treatment + (treatment|genotype) +
                           (1|block) + (1|block:treatment),
                         data = filter(wind.df, Year == "2013"))
summary(height.2013.lmer)
summary(update(height.2013.lmer, .~. -1))
1-(11.648/22.62)
coef(update(height.2013.lmer, .~. -1))

plot(height.2013.lmer) # looks okay, note that logging qualitatively alters the conclusion about GxE.

# test for GxE
anova(height.2013.lmer,
      update(height.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(height.2013.lmer, .~. -(treatment|genotype)))

anova(height.2013.lmer, ddf = "Kenward-Roger")
visreg(height.2013.lmer, xvar = "treatment")

# random model
height.2013.rand <- lmer(Height ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2013.df)
summary(height.2013.rand)
56.295/(19.083+16.579+1.176+23.105+56.295+64.671) # wind exposure explained 31% of the variance in plant height
1.176/(19.083+16.579+1.176+23.105+56.295+64.671) # genotype explained <1% of the variance in plant height
23.105/(19.083+16.579+1.176+23.105+56.295+64.671) # GxE explained 13% of the variance in plant height.
height.2013.vars <- as.data.frame(VarCorr(height.2013.rand))

# shoot count
hist(arch.pca.2013.df$all.shoot.count)
hist(sqrt(arch.pca.2013.df$all.shoot.count))
hist(log(arch.pca.2013.df$all.shoot.count))

shoot.count.2013.lmer <- lmer(all.shoot.count ~ treatment + (treatment|genotype) +
                                (1|block) + (1|block:treatment),
                              data = arch.pca.2013.df)
summary(shoot.count.2013.lmer)
summary(update(shoot.count.2013.lmer, .~. -1))
1-(3.4655/5.7586)
plot(shoot.count.2013.lmer)  

# test for GxE
anova(shoot.count.2013.lmer,
      update(shoot.count.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(shoot.count.2013.lmer, .~. -(treatment|genotype))) # marginally close to GxE (should divide p-value by 2)

anova(shoot.count.2013.lmer, ddf = "Kenward-Roger")
visreg(shoot.count.2013.lmer, xvar = "treatment")

# random model
shoot.count.2013.rand <- lmer(all.shoot.count ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2013.df)
summary(shoot.count.2013.rand)
2.5598/(1.7568+1.3476+0.8464+2.5598+7.9179) # wind exposure explained 18% of the variance
0.8464/(1.7568+1.3476+0.8464+2.5598+7.9179) # genotype explained 6% of the variance
shoot.count.2013.vars <- as.data.frame(VarCorr(shoot.count.2013.rand))

# average shoot length
hist(filter(wind.df, Year == "2013")$all.shoot.avg.length)
avg.shoot.length.2013.lmer <- lmer(log(all.shoot.avg.length) ~
                                     treatment + (treatment|genotype) +
                                 (1|block) + (1|block:treatment),
                               data = filter(wind.df, Year == "2013"))
summary(avg.shoot.length.2013.lmer)
summary(lmer(all.shoot.avg.length ~
               treatment + (treatment|genotype) - 1 +
               (1|block) + (1|block:treatment),
             data = filter(wind.df, Year == "2013")))
1-(0.6018/0.8429) # 29% shorter in wind exposed plots
coef(lmer(all.shoot.avg.length ~
               treatment + (treatment|genotype) +
               (1|block) + (1|block:treatment),
             data = filter(wind.df, Year == "2013")))
0.8221796/0.3622077
plot(avg.shoot.length.2013.lmer) # much better with log-transformation

# test for GxE
anova(avg.shoot.length.2013.lmer,
      update(avg.shoot.length.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(avg.shoot.length.2013.lmer, .~. -(treatment|genotype)))

anova(avg.shoot.length.2013.lmer, ddf = "Kenward-Roger")
visreg(avg.shoot.length.2013.lmer, xvar = "treatment")

# random model
avg.shoot.length.2013.rand <- lmer(log(all.shoot.avg.length) ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = filter(wind.df, Year == "2013"))
summary(avg.shoot.length.2013.rand)
0.06784/(0.09846+0.29431+0.06784+0.04201+0.29892) # genotype explained 8% of the total variance
0.04201/(0.09846+0.29431+0.06784+0.04201+0.29892) # wind exposure explained 5% of the total variance
avg.shoot.length.2013.vars <- as.data.frame(VarCorr(avg.shoot.length.2013.rand))

# shoot length
shoot.length.2013.lmer <- lmer(log(all.shoot.total.length) ~ treatment + (treatment|genotype) +
                                 (1|block) + (1|block:treatment),
                               data = arch.pca.2013.df)
summary(shoot.length.2013.lmer)
plot(shoot.length.2013.lmer) # much better with log-transformation

# test for GxE
anova(shoot.length.2013.lmer,
      update(shoot.length.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(shoot.length.2013.lmer, .~. -(treatment|genotype)))

anova(shoot.length.2013.lmer, ddf = "Kenward-Roger")
visreg(shoot.length.2013.lmer, xvar = "treatment")

# random model
shoot.length.2013.rand <- lmer(all.shoot.total.length ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = arch.pca.2013.df)
summary(shoot.length.2013.rand)
shoot.length.2013.vars <- as.data.frame(VarCorr(shoot.length.2013.rand))

mean(c(31,18,5)) # 18% of variance for wind exposure
mean(c(8,6,13)) # 9% of variance for G and/or GxE

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
summary(LQ.PCA.2013) # first two components have an eigenvalue > 1
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
LQ.pc1.2013.lmer <- lmer(LQ.PC1 ~ treatment + (treatment|genotype) +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2013.df)
summary(LQ.pc1.2013.lmer)
plot(LQ.pc1.2013.lmer) # looks okay

# test for GxE
anova(LQ.pc1.2013.lmer,
      update(LQ.pc1.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(LQ.pc1.2013.lmer, .~. -(treatment|genotype)))

anova(LQ.pc1.2013.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc1.2013.lmer, xvar = "treatment")

# random model
LQ.pc1.2013.rand <- lmer(LQ.PC1 ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = LQ.pca.2013.df)
summary(LQ.pc1.2013.rand )
LQ.pc1.2013.vars <- as.data.frame(VarCorr(LQ.pc1.2013.rand))

# WC
leaf_WC.2013.lmer <- lmer(log(leaf_WC) ~ treatment + (treatment|genotype) +
                            (1|block) + (1|block:treatment),
                          data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_WC.2013.lmer)
plot(leaf_WC.2013.lmer)  # better with log-transformation

# test for GxE
anova(leaf_WC.2013.lmer,
      update(leaf_WC.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(leaf_WC.2013.lmer, .~. -(treatment|genotype)))

anova(leaf_WC.2013.lmer, ddf = "Kenward-Roger")
visreg(leaf_WC.2013.lmer, xvar = "treatment")

# random model
leaf_WC.2013.rand <- lmer(log(leaf_WC) ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_WC.2013.rand )
leaf_WC.2013.vars <- as.data.frame(VarCorr(leaf_WC.2013.rand))

# LQ PC2
LQ.pc2.2013.lmer <- lmer(LQ.PC2 ~ treatment + (treatment|genotype) +
                           (1|block) + (1|block:treatment),
                         data = LQ.pca.2013.df)
summary(LQ.pc2.2013.lmer)
plot(LQ.pc2.2013.lmer) # looks okay

# test for GxE
anova(LQ.pc2.2013.lmer,
      update(LQ.pc2.2013.lmer, .~. -(treatment|genotype) + (1|genotype)), # Genotype effect would be significant when dividing p-value by 2, since this test is on the boundary.
      update(LQ.pc2.2013.lmer, .~. -(treatment|genotype)))

anova(LQ.pc2.2013.lmer, ddf = "Kenward-Roger")
visreg(LQ.pc2.2013.lmer, xvar = "treatment")

# random model
LQ.pc2.2013.rand <- lmer(LQ.PC2 ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = LQ.pca.2013.df)
summary(LQ.pc2.2013.rand )
LQ.pc2.2013.vars <- as.data.frame(VarCorr(LQ.pc2.2013.rand))

# leaf C:N
leaf_CN.2013.lmer <- lmer(log(leaf_C_N) ~ treatment + (treatment|genotype) +
                            (1|block) + (1|block:treatment),
                          data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_CN.2013.lmer)
plot(leaf_CN.2013.lmer)   # residuals a bit better with log transformation

# test for GxE
anova(leaf_CN.2013.lmer,
      update(leaf_CN.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(leaf_CN.2013.lmer, .~. -(treatment|genotype)))

anova(leaf_CN.2013.lmer, ddf = "Kenward-Roger")
visreg(leaf_CN.2013.lmer, xvar = "treatment")

# random model
leaf_CN.2013.rand <- lmer(log(leaf_C_N) ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(leaf_CN.2013.rand )
leaf_CN.2013.vars <- as.data.frame(VarCorr(leaf_CN.2013.rand)) %>%
  select(group = grp, var = vcov, sd = sdcor) %>%
  mutate(response = "log.leaf_CN",
         perc.var = round(var/sum(var),2)) %>%
  tbl_df()

# SLA
SLA.2013.lmer <- lmer(SLA ~ treatment + (treatment|genotype) +
                            (1|block) + (1|block:treatment),
                          data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(SLA.2013.lmer)
plot(SLA.2013.lmer)   # residuals a bit better with log transformation

# test for GxE
anova(SLA.2013.lmer,
      update(SLA.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(SLA.2013.lmer, .~. -(treatment|genotype)))

anova(SLA.2013.lmer, ddf = "Kenward-Roger")
visreg(SLA.2013.lmer, xvar = "treatment")

# random model
SLA.2013.rand <- lmer(SLA ~  (1|treatment) + (1|genotype) + (1|block) + (1|block:treatment), data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(SLA.2013.rand )
SLA.2013.vars <- as.data.frame(VarCorr(SLA.2013.rand)) %>%
  select(group = grp, var = vcov, sd = sdcor) %>%
  mutate(response = "SLA",
         perc.var = round(var/sum(var),2)) %>%
  tbl_df()

# larva wet weight exp. 2
hist(filter(wind.df, Year == "2013", larva.wet.wt.exp2 != "NA")[-WC_nonsense, ]$larva.wet.wt.exp2)
hist(log(filter(wind.df, Year == "2013", larva.wet.wt.exp2 != "NA")[-WC_nonsense, ]$larva.wet.wt.exp2+1))

with(filter(wind.df, Year == "2013", larva.wet.wt.exp2 != "NA")[-WC_nonsense, ],
     interaction.plot(treatment, genotype, larva.wet.wt.exp2))

larva.wet.wt.exp2.2013.lmer <- lmer(log(larva.wet.wt.exp2+1) ~ treatment + (treatment|genotype) +
                        (1|block) + (1|block:treatment),
                      data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(larva.wet.wt.exp2.2013.lmer)
plot(larva.wet.wt.exp2.2013.lmer)   # residuals a bit better with log transformation

larva.wet.wt.exp2.2013.glmer <- glmer(larva.wet.wt.exp2 ~ treatment + (treatment|genotype) +
                                      (1|block) + (1|block:treatment),
                                    data = filter(wind.df, Year == "2013")[-WC_nonsense, ], family = "poisson")
summary(larva.wet.wt.exp2.2013.glmer)
plot(larva.wet.wt.exp2.2013.glmer)

# test for GxE
anova(larva.wet.wt.exp2.2013.lmer,
      update(larva.wet.wt.exp2.2013.lmer, .~. -(treatment|genotype) + (1|genotype)),
      update(larva.wet.wt.exp2.2013.lmer, .~. -(treatment|genotype))) # likely significant after dividing p-value by 2 since it is on the boundary

anova(larva.wet.wt.exp2.2013.lmer, ddf = "Kenward-Roger")
visreg(larva.wet.wt.exp2.2013.lmer, xvar = "treatment")

# random model
larva.wet.wt.exp2.2013.rand <- lmer(log(larva.wet.wt.exp2+1) ~  (1|treatment) + (treatment|genotype) + (1|block) + (1|block:treatment), data = filter(wind.df, Year == "2013")[-WC_nonsense, ])
summary(larva.wet.wt.exp2.2013.rand )
larva.wet.wt.exp2.2013.vars <- as.data.frame(VarCorr(larva.wet.wt.exp2.2013.rand)) %>%
  filter(vcov > 0) %>%
  mutate(grp = paste(grp, var1, sep = "_")) %>%
  select(group = grp, var = vcov, sd = sdcor) %>%
  mutate(response = "SLA",
         perc.var = round(var/sum(var),2)) %>%
  tbl_df()
   
  

## Ant-aphid 2012: architecture ----
aa.arch.traits <- c("Height","all.shoot.count","mature.shoot.total.length")

# explore correlations
scatterplotMatrix(filter(aa.df, Year == "2012")[ ,aa.arch.traits])

# PCA: 2012
aa.arch.pca.2012.df <- filter(aa.df, Year == "2012") %>%
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
aa.pc1.2012.lmer <- lmer(arch.PC1 ~ Aphid.Treatment*Ant.Mound.Dist + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                      data = aa.arch.pca.2012.df)
plot(aa.pc1.2012.lmer) # looks pretty good
summary(aa.pc1.2012.lmer)

# test for GxE
anova(aa.pc1.2012.lmer,
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)),
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)),
      update(aa.pc1.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

anova(aa.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.pc1.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.pc1.2012.rand <- lmer(arch.PC1 ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data = aa.arch.pca.2012.df)
summary(aa.pc1.2012.rand)
aa.pc1.2012.vars <- as.data.frame(VarCorr(aa.pc1.2012.rand))


# plant height
aa.height.2012.lmer <- lmer(Height ~ Aphid.Treatment*Ant.Mound.Dist + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data = aa.arch.pca.2012.df)
plot(aa.height.2012.lmer) # not great...
summary(aa.height.2012.lmer)

# test for GxE
anova(aa.height.2012.lmer,
      update(aa.height.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)),
      update(aa.height.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)),
      update(aa.height.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

anova(aa.height.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.height.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.height.2012.rand <- lmer(Height ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                           (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data = aa.arch.pca.2012.df)
summary(aa.height.2012.rand)
aa.height.2012.vars <- as.data.frame(VarCorr(aa.height.2012.rand))

# shoot count
aa.shoot.count.2012.lmer <- lmer(log(all.shoot.count) ~ Aphid.Treatment*Ant.Mound.Dist + 
                              (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                              (1|Block) + (1|Block:Ant.Mound.Dist),
                            data = aa.arch.pca.2012.df)
plot(aa.shoot.count.2012.lmer) # not great...
summary(aa.shoot.count.2012.lmer)

# test for GxE
anova(aa.shoot.count.2012.lmer,
      update(aa.shoot.count.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)),
      update(aa.shoot.count.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)),
      update(aa.shoot.count.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

anova(aa.shoot.count.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.shoot.count.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.shoot.count.2012.rand <- lmer(log(all.shoot.count) ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                              (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                              (1|Block) + (1|Block:Ant.Mound.Dist),
                            data = aa.arch.pca.2012.df)
summary(aa.shoot.count.2012.rand)
aa.shoot.count.2012.vars <- as.data.frame(VarCorr(aa.shoot.count.2012.rand))

# shoot length. 
aa.mature.shoot.total.length.2012.lmer <- lmer(mature.shoot.total.length ~ Aphid.Treatment*Ant.Mound.Dist + 
                                   (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                                   (1|Block) + (1|Block:Ant.Mound.Dist),
                                 data = aa.arch.pca.2012.df)
plot(aa.mature.shoot.total.length.2012.lmer) # not great...
summary(aa.mature.shoot.total.length.2012.lmer)

# test for GxE
anova(aa.mature.shoot.total.length.2012.lmer,
      update(aa.mature.shoot.total.length.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)),
      update(aa.mature.shoot.total.length.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)),
      update(aa.mature.shoot.total.length.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

anova(aa.mature.shoot.total.length.2012.lmer, ddf = "Kenward-Roger") #convergence problems...
visreg(aa.mature.shoot.total.length.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.mature.shoot.total.length.2012.rand <- lmer(mature.shoot.total.length ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                                   (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                                   (1|Block) + (1|Block:Ant.Mound.Dist),
                                 data = aa.arch.pca.2012.df)
summary(aa.mature.shoot.total.length.2012.rand)
aa.mature.shoot.total.length.2012.vars <- as.data.frame(VarCorr(aa.mature.shoot.total.length.2012.rand))

## Ant-aphid 2012: leaf quality ----

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

# decided from correlation plots that the non-logged transformed relationships appear to be more linear, so I've decided to go with those ones.
aa.pca.log.LQ.traits <- c("log_WC", "log_trichomes")
scatterplotMatrix(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])
corr.test(aa.LQ.pca.2012.df[ ,aa.pca.log.LQ.traits])

aa.LQ.PCA <- princomp(aa.LQ.pca.2012.df[ ,c("leaf_WC","leaf_trichome.density")], cor = TRUE) 
summary(aa.LQ.PCA) # only first component has an eigenvalue > 1
plot(aa.LQ.PCA)
aa.LQ.PCA$loadings # PC1 represents the correlation between leaf water content and trichome density.
aa.LQ.PCA$scores[ ,"Comp.1"]

aa.LQ.PCA.df <- mutate(aa.LQ.pca.2012.df, 
                         LQ.PC1 = -1*aa.LQ.PCA$scores[ ,"Comp.1"]) # multiply by -1 so positive values indicate higher leaf water content and higher trichome density

hist(aa.LQ.PCA.df$LQ.PC1) # a bit skewed, but transformation should help
hist(log(aa.LQ.PCA.df$LQ.PC1+2))

## Analysis

with(aa.LQ.PCA.df, table(Genotype, Ant.Mound.Dist, Aphid.Treatment)) # testing 3-way interaction is unwise
# 2-way interactions mostly unbalanced
with(aa.LQ.PCA.df, table(Genotype, Ant.Mound.Dist)) 
with(aa.LQ.PCA.df, table(Genotype, Aphid.Treatment)) 
with(aa.LQ.PCA.df, table(Aphid.Treatment, Ant.Mound.Dist)) # this one isn't too bad.
# 1-way effects
with(aa.LQ.PCA.df, table(Genotype)) 

contr.sum(aa.LQ.PCA.df$Genotype)

## Analysis

# LQ PC1. Ran without 3-way interaction because there were convergence problems
aa.LQ.pc1.2012.lmer <- lmer(LQ.PC1 ~ Aphid.Treatment*Ant.Mound.Dist + 
                           (Aphid.Treatment + Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data =  aa.LQ.PCA.df)
plot(aa.LQ.pc1.2012.lmer) # looks pretty good
summary(aa.LQ.pc1.2012.lmer)

# test for GxE
anova(aa.LQ.pc1.2012.lmer,
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype) + (Ant.Mound.Dist|Genotype)), # replaced with Ant.Mound.Dist because it had the highest variance of the possible interactions
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype) + (1|Genotype)), # likely significance of Genotype after dividing p-value by 2 to account for boundary effect.
      update(aa.LQ.pc1.2012.lmer, .~. -(Aphid.Treatment + Ant.Mound.Dist|Genotype)))

anova(aa.LQ.pc1.2012.lmer, ddf = "Kenward-Roger")
visreg(aa.LQ.pc1.2012.lmer, xvar = "Ant.Mound.Dist", by = "Aphid.Treatment")

# random model
aa.LQ.pc1.2012.rand <- lmer(LQ.PC1 ~ (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                           (Aphid.Treatment + Ant.Mound.Dist|Genotype) +
                           (1|Block) + (1|Block:Ant.Mound.Dist),
                         data =  aa.LQ.PCA.df)
summary(aa.LQ.pc1.2012.rand)

aa.LQ.pc1.2012.vars <- as.data.frame(VarCorr(aa.LQ.pc1.2012.rand))


# WC
aa.leaf_WC.2012.lmer <- lmer(log(leaf_WC) ~ Aphid.Treatment*Ant.Mound.Dist + 
                            (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                            (1|Block) + (1|Block:Ant.Mound.Dist),
                          data = filter(aa.df, Year == "2012"))
summary(aa.leaf_WC.2012.lmer)
plot(aa.leaf_WC.2012.lmer) # looks okay

# test for GxE
anova(aa.leaf_WC.2012.lmer,
      update(aa.leaf_WC.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)), 
      update(aa.leaf_WC.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)), 
      update(aa.leaf_WC.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))

# random model
aa.leaf_WC.2012.rand <- lmer(log(leaf_WC) ~  (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                               (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                               (1|Block) + (1|Block:Ant.Mound.Dist), data = filter(aa.df, Year == "2012"))
summary(aa.leaf_WC.2012.rand )
aa.leaf_WC.2012.vars <- as.data.frame(VarCorr(aa.leaf_WC.2012.rand))

# trichome density
aa.leaf_trichome.density.2012.lmer <- lmer(log(leaf_trichome.density+1) ~ Aphid.Treatment*Ant.Mound.Dist + 
                                             (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                                             (1|Block) + (1|Block:Ant.Mound.Dist),
                                        data = filter(aa.df, Year == "2012"))
summary(aa.leaf_trichome.density.2012.lmer)
plot(aa.leaf_trichome.density.2012.lmer) # not great, could consider glmer poisson

# test for GxE
anova(aa.leaf_trichome.density.2012.lmer,
      update(aa.leaf_trichome.density.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (Aphid.Treatment|Genotype)), 
      update(aa.leaf_trichome.density.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype) + (1|Genotype)), 
      update(aa.leaf_trichome.density.2012.lmer, .~. -(Aphid.Treatment*Ant.Mound.Dist|Genotype)))


anova(aa.leaf_trichome.density.2012.lmer, ddf = "Kenward-Roger")

# random model
aa.leaf_trichome.density.2012.rand <- lmer(log(leaf_trichome.density+1) ~  (1|Aphid.Treatment) + (1|Ant.Mound.Dist) + (1|Aphid.Treatment:Ant.Mound.Dist) + 
                                             (Aphid.Treatment*Ant.Mound.Dist|Genotype) +
                                             (1|Block) + (1|Block:Ant.Mound.Dist), data = filter(aa.df, Year == "2012"))
summary(aa.leaf_trichome.density.2012.rand )
aa.leaf_trichome.density.2012.vars <- as.data.frame(VarCorr(aa.leaf_trichome.density.2012.rand))

