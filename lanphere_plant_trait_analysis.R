
## load required libraries and fucnctions ----
source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(dplyr)
library(reshape2)
library(ggplot2)
#library(visreg)
#library(psych)
#library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
#library(lmerTest)
library(piecewiseSEM) # for calculating r2 of mixed-effect models.
library(lme4)
library(car) # for Anova function
library(afex) # for calculating p-values of mixed effect models
library(broom) # for tidying up model outputs
#library(blme)

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
  mutate(fact.Ant.Mound.Dist = as.factor(Ant.Mound.Dist),
         X = as.factor(X))
glimpse(aa.df)

## Survival analysis ----

# explore 2-way interactions
with(wind.df, table(Dead, genotype, treatment))
with(wind.df, table(Dead, genotype, Year))
with(wind.df, table(Dead, treatment, Year))

# explore main effects
with(wind.df, table(Dead, treatment))
with(wind.df, table(Dead, genotype))
with(wind.df, table(Dead, Year))

# GLMM. Modelled genotype as random effect because the model was having difficulty converging.
surv.glmer <- glmer(Dead ~ treatment + Year + genotype +
                            (1|block) + 
                            (1|block:treatment) +
                            (1|plant.id),
                          data = wind.df,
                          contrasts = list(treatment = "contr.sum",
                                           genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "binomial",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(surv.glmer)
# overdisp_fun() # overdispersion test unnecessary, since the response is binary.
plot(surv.glmer) 

# Likelihood ratio tests for main effects
drop1(surv.glmer, test = "Chisq") # no ExE
drop1(update(surv.glmer, .~. -treatment:Year), test = "Chisq") # only Year effect

# no G or GxE effect
anova(surv.glmer, 
      update(surv.glmer, .~. -(treatment|genotype) + (1|genotype)),
      update(surv.glmer, .~. -(treatment|genotype)))

## Wind: plant height analysis ----
hist(wind.df$Height)
height.outlier <- which(wind.df$Height > 80)  # not a typo according to data entry, but definitely seems to be an outlier. Maybe it was written down incorrectly in the field, because it is almost twice the size of the other plants. I've decided to remove it from this analysis.

# LMMM. 
height.lmer <- lmer(Height ~ treatment*Year*genotype + 
                        (1|block) + 
                        (1|block:treatment) +
                        (1|plant.id),
                    data = wind.df[-height.outlier, ],
                    contrasts = list(treatment = "contr.sum",
                                     Year = "contr.sum",
                                     genotype = "contr.sum"))
summary(height.lmer)
plot(height.lmer) # diamond shaped...

# Kenward-Roger test
height.anova <- anova.table(height.lmer, experiment = "wind")

# calculate variance components for simplified model
height.up <- update(height.lmer, .~. -treatment*Year*genotype + (1|genotype) + treatment*Year)

height.R2 <- var.table(height.up, experiment = "wind")

## Wind: shoot count analysis ----

# GLMM
shoot.count.glmer <- glmer(all.shoot.count ~ treatment*genotype*Year +
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
plot(shoot.count.glmer)

## Wald-tests
shoot.count.anova <- anova.table(shoot.count.glmer, test = "Chisq", type = 3, experiment = "wind")

# calculate variance components for simplified model
shoot.count.up <- update(shoot.count.glmer, .~. -treatment*Year*genotype + (1|genotype) + treatment*Year)

shoot.count.R2 <- var.table(shoot.count.up, experiment = "wind")

## Wind: average shoot length analysis ----

# LMM 
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

# Kenward-Roger tests 
shoot.length.anova <- anova.table(shoot.length.lmer, experiment = "wind")

## Calculate R2 for significant predictors
# no interaction effects were significant so I didn't calculate R2
shoot.length.up <- update(shoot.length.lmer, .~. -treatment*genotype*Year + treatment + Year + (1|genotype))

shoot.length.R2 <- var.table(shoot.length.up, experiment = "wind")

## Wind: leaf water content ----
hist(wind.df$leaf_WC)
WC_nonsense <- which(wind.df$leaf_WC < 0)

leaf_WC.lmer <- lmer(log(leaf_WC) ~ treatment*genotype*Year +
                            (1|block) + 
                            (1|block:treatment) +
                            (1|plant.id),
                         data = wind.df[-WC_nonsense, ],
                         contrast = list(treatment = "contr.sum",
                                         genotype = "contr.sum",
                                         Year = "contr.sum"))
summary(leaf_WC.lmer)  # note that plant.id explains the least amount of variance
plot(leaf_WC.lmer) # looks good

# Kenward-Roger test
leaf_WC.anova <- anova.table(leaf_WC.lmer, experiment = "wind")

## Calculate R2
leaf_WC.up <- update(leaf_WC.lmer, .~. -treatment*genotype*Year + Year + (Year|genotype) - (1|plant.id)) # removed plant.id as random effect because it explained very little of the variance in the big model and it enabled this model to converge. Note though that the results are qualitatively similar if I remove other random effects to enable the model to converge 
summary(leaf_WC.up)

leaf_WC.R2 <- var.table(leaf_WC.up, experiment = "wind")

## Wind: trichome density analysis ----
leaf_trichome.density.glmer <- glmer(
  leaf_trichome.density ~ treatment*genotype +
    (1|block) + 
    (1|block:treatment) +
    (1|plant.id), # need individual-level random effect to account for overdispersion
  data = wind.df, 
  family = "poisson",
  contrasts = list(treatment = "contr.sum",
                   genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(leaf_trichome.density.glmer)
overdisp_fun(leaf_trichome.density.glmer)
plot(leaf_trichome.density.glmer) # becomes rather skewed when I model the individual-level random effect

# Likelihood-ratio tests
trichome.density.anova <- anova.table(leaf_trichome.density.glmer, test = "Chisq", experiment = "wind")

## Calculate R2
leaf_trichome.density.up <- update(leaf_trichome.density.glmer, .~. -treatment*genotype + (1|genotype) - (1|block)) # removed block as well because the updated model with block didn't explain any of the variance (therefore, over-specified). Note that this doesn't change the variance explained by Genotype either way.

leaf_trichome.density.R2 <- var.table(leaf_trichome.density.up, experiment = "wind")

## Wind: leaf C:N analysis ----

# LMM
leaf_CN.lmer <- lmer(log(leaf_C_N) ~ treatment*genotype +
                            (1|block) + (1|block:treatment),
                          data = wind.df,
                     contrasts = list(treatment = "contr.sum",
                                      genotype = "contr.sum"))
summary(leaf_CN.lmer)
plot(leaf_CN.lmer) # looks pretty good

## Kenwar-Roger test
leaf_CN.anova <- anova.table(leaf_CN.lmer, experiment = "wind")

## Calculate R2
leaf_CN.up <- update(leaf_CN.lmer, .~. -treatment*genotype + (1|genotype))

leaf_CN.R2 <- var.table(leaf_CN.up, experiment = "wind")

## Wind: SLA analysis ----
hist(wind.df$SLA)

# LMM
SLA.lmer <- lmer(SLA ~ treatment*genotype +
                            (1|block) + (1|block:treatment),
                          data = wind.df,
                      contrasts = list(treatment = "contr.sum",
                                       genotype = "contr.sum"))
summary(SLA.lmer)
plot(SLA.lmer) # residuals look okay

## Kenward-Roger test
SLA.anova <- anova.table(SLA.lmer, experiment = "wind")

## Calculate R2
SLA.up <- update(SLA.lmer, .~. -treatment*genotype + (1|genotype))

SLA.R2 <- var.table(SLA.up, experiment = "wind")

## Wind: larva wet weight exp. 1 analysis ----
hist(wind.df$larva.wet.wt.exp1)

with(filter(wind.df, larva.wet.wt.exp1 != "NA"),
     interaction.plot(treatment, genotype, larva.wet.wt.exp1))

# convergence problem if I included GxE, so I only modelled the main effects
larva.wet.wt.exp1.glmer <- glmer(
  larva.wet.wt.exp1 ~ treatment + genotype +
    (1|block) + (1|block:treatment),
  data = wind.df, family = "poisson",
  contrasts = list(treatment = "contr.sum",
                   genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(larva.wet.wt.exp1.glmer)
overdisp_fun(larva.wet.wt.exp1.glmer) # no overdispersion
plot(larva.wet.wt.exp1.glmer) # weird residuals...

## Likelihood ratio tests
drop1(larva.wet.wt.exp1.glmer, test = "Chisq") # no G or E

## Wind: larva wet weight exp. 2 analysis ----
hist(wind.df$larva.wet.wt.exp2)

with(filter(wind.df, larva.wet.wt.exp2 != "NA"),
     interaction.plot(treatment, genotype, larva.wet.wt.exp2))

larva.wet.wt.exp2.glmer <- glmer(
  larva.wet.wt.exp2 ~ treatment*genotype +
    (1|block) + (1|block:treatment),
  data = wind.df, family = "poisson",
  contrasts = list(treatment = "contr.sum",
                   genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(larva.wet.wt.exp2.glmer)
overdisp_fun(larva.wet.wt.exp2.glmer) # no overdispersion
plot(larva.wet.wt.exp2.glmer) # weird residuals...

## Likelihood ratio tests
drop1(larva.wet.wt.exp2.glmer, test = "Chisq") # no GxE
drop1(update(larva.wet.wt.exp2.glmer, .~. -treatment:genotype), 
             test = "Chisq") # no genotype or treatment effect

## Ant-aphid: plant height analysis ----

# LMM
aa.height.lmer <- lmer(Height ~ Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype +
                           (1|Block) + (1|Block:fact.Ant.Mound.Dist),
                         data = filter(aa.df, Year == "2012"),
                         contrasts = list(Aphid.Treatment = "contr.sum",
                                          Genotype = "contr.sum"))
summary(aa.height.lmer)
plot(aa.height.lmer) # not great...

## Kenward-Roger test
aa.height.anova <- anova.table(aa.height.lmer, experiment = "ant-aphid") 

## Calculate R2
aa.height.up <- update(aa.height.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Genotype))

aa.height.R2 <- var.table(aa.height.up, experiment = "ant-aphid")

## Ant-aphid: Shoot count analysis ----

# GLMM
aa.shoot.count.glmer <- glmer(
  all.shoot.count ~ Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + 
    (1|Block) + (1|Block:fact.Ant.Mound.Dist),
                             data = filter(aa.df, Year == "2012"),
                             family = "poisson",
                             contrasts = list(Aphid.Treatment = "contr.sum",
                                             Genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
plot(aa.shoot.count.glmer) # not great
summary(aa.shoot.count.glmer)
overdisp_fun(aa.shoot.count.glmer) # no overdispersion

## Wald Chi-square tests
aa.shoot.count.anova <- anova.table(aa.shoot.count.glmer, test = "Chisq", experiment = "ant-aphid")

## Calculate R2 for significant predictors
aa.shoot.count.up <- update(aa.shoot.count.glmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Genotype))

aa.shoot.count.R2 <- var.table(aa.shoot.count.up, experiment = "ant-aphid")

## Ant-Aphid: Average shoot length analysis ----
hist(aa.df$mature.shoot.avg.length)

aa.mature.shoot.avg.length.lmer <- lmer(mature.shoot.avg.length ~ Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Block) + (1|Block:Ant.Mound.Dist),
                                 data = filter(aa.df, Year == "2012"),
                                 contrasts = list(Aphid.Treatment = "contr.sum",
                                                  Genotype = "contr.sum"))
plot(aa.mature.shoot.avg.length.lmer) # not great...
summary(aa.mature.shoot.avg.length.lmer)

## Kenward-Roger
aa.mature.shoot.avg.length.anova <- anova.table(aa.mature.shoot.avg.length.lmer, experiment = "ant-aphid")

## Calculate R2
aa.mature.shoot.avg.length.up <- update(aa.mature.shoot.avg.length.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Genotype))

aa.mature.shoot.avg.length.R2 <- var.table(aa.mature.shoot.avg.length.up, experiment = "ant-aphid")

## Ant-Aphid: leaf water content analysis ----
hist(aa.df$leaf_WC)

## LMM. Had to fix a less complex model because we didn't have enough replication of all treatment combinations
aa.leaf_WC.lmer <- lmer(log(leaf_WC) ~ Aphid.Treatment*scale(Ant.Mound.Dist) + scale(Ant.Mound.Dist)*Genotype + (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                        data = filter(aa.df, Year == "2012"),
                        contrasts = list(Aphid.Treatment = "contr.sum",
                                         Genotype = "contr.sum"))
summary(aa.leaf_WC.lmer)
plot(aa.leaf_WC.lmer) # looks okay

## Kenward-Roger test
aa.leaf_WC.anova <- anova.table(aa.leaf_WC.lmer, experiment = "ant-aphid") # nothing significant, but genotype is the closest.

## Calculate R2
# Genotype was not significant, but I'm still going to calculate its R2.
aa.leaf_WC.up <- update(aa.leaf_WC.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist) - scale(Ant.Mound.Dist)*Genotype + (1|Genotype))

anova(aa.leaf_WC.up, update(aa.leaf_WC.up, .~. -(1|Genotype)))

aa.leaf_WC.R2 <- var.table(aa.leaf_WC.up, experiment = "ant-aphid")

## Ant-aphid: leaf trichome density ----

# GLMM
aa.leaf_trichome.density.glmer <- glmer(leaf_trichome.density ~ Aphid.Treatment*scale(Ant.Mound.Dist) + scale(Ant.Mound.Dist)*Genotype + 
                                          (1|Block) + 
                                          (1|X) +
                                          (1|Block:Ant.Mound.Dist),
                                        data = filter(aa.df, Year == "2012"),
                                        family = "poisson",
                                        contrasts = list(Aphid.Treatment = "contr.sum",
                                                         Genotype = "contr.sum"),
                                        control=glmerControl(optimizer="bobyqa",
                                                             optCtrl=list(maxfun=2e5)))
summary(aa.leaf_trichome.density.glmer)
plot(aa.leaf_trichome.density.glmer) # not great
overdisp_fun(aa.leaf_trichome.density.glmer) # overdispersed so I modelled an individual-level random effect

## Wald Chi-square tests
aa.trichome.density.anova <- anova.table(aa.leaf_trichome.density.glmer, test = "Chisq", experiment = "ant-aphid")

## Calculate R2
aa.leaf_trichome.density.up <- update(aa.leaf_trichome.density.glmer, .~. -Aphid.Treatment*scale(Ant.Mound.Dist) - scale(Ant.Mound.Dist)*Genotype + (1|Genotype))

aa.leaf_trichome.density.R2 <- var.table(aa.leaf_trichome.density.up, experiment = "ant-aphid")

## Print Summary ----
trait.anovas <- bind_rows(height.anova, 
                          shoot.count.anova, 
                          shoot.length.anova,
                          leaf_WC.anova,
                          SLA.anova,
                          trichome.density.anova,
                          leaf_CN.anova,
                          aa.height.anova,
                          aa.shoot.count.anova,
                          aa.mature.shoot.avg.length.anova,
                          aa.leaf_WC.anova,
                          aa.trichome.density.anova) %>%
  dplyr::select(experiment, response, term, test, df, Df.res, statistic, p.value, Pr..Chisq.)

trait.R2s <- bind_rows(height.R2,
                       shoot.count.R2,
                       shoot.length.R2,
                       leaf_WC.R2,
                       SLA.R2,
                       leaf_trichome.density.R2,
                       leaf_CN.R2,
                       aa.height.R2,
                       aa.shoot.count.R2,
                       aa.mature.shoot.avg.length.R2,
                       aa.leaf_WC.R2,
                       aa.leaf_trichome.density.R2)
#bind_cols(trait.anovas, trait.R2s)
write.csv(trait.anovas, "trait.anovas.csv")
write.csv(trait.R2s, "trait.R2s.csv")

ggplot(trait.R2s, aes(x = response, y = var_percent, fill = Factors)) + geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() + 
  facet_grid(~experiment) 

## Old ----
## Wind: Plant architecture 2012

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

## Wind: Leaf quality 2012 

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

## Wind plant architecture 2013

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


## Wind: Leaf quality traits 2013 
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

## Ant-aphid 2012: architecture 
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

