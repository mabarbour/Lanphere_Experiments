
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

anova.table(height.lmer, experiment = "wind")

# calculate variance components for simplified model
height.up <- update(height.lmer, .~. -treatment*Year*genotype + (1|genotype) + treatment*Year)

var.table(height.up, "wind")
plotREsim(REsim(height.up))
sem.model.fits(height.up)

# Kenward-Roger test
height.anova <- anova(height.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "plant height",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "F (Kenward-Roger)",
         Experiment = "wind",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2 for significant predictors
height.2x <- update(height.lmer, .~. -treatment:genotype:Year)
height.main <- update(height.lmer, .~. -treatment*genotype*Year + treatment + genotype + Year)

height.R2 <- 
  data.frame(Factor = c("treatment:Year:genotype","treatment:Year","treatment","genotype","Year"),
             Response = "plant height",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(height.lmer, update(height.2x, .~. -treatment:Year:genotype)),
                                deltaR2(height.2x, update(height.2x, .~. -treatment:Year)),
                                deltaR2(height.main, update(height.main, .~. -treatment)),
                                deltaR2(height.main, update(height.main, .~. -genotype)),
                                deltaR2(height.main, update(height.main, .~. -Year))),2))

## Wind: shoot count analysis ----

# explore 2-way interactions
with(filter(wind.df, all.shoot.count > 0), interaction.plot(treatment, genotype, all.shoot.count))
with(filter(wind.df, all.shoot.count > 0), interaction.plot(treatment, Year, all.shoot.count))
with(filter(wind.df, all.shoot.count > 0), interaction.plot(Year, genotype, all.shoot.count))

# explore main effects
with(filter(wind.df, all.shoot.count > 0), plot(all.shoot.count ~ treatment))
with(filter(wind.df, all.shoot.count > 0), plot(all.shoot.count ~ genotype))
with(filter(wind.df, all.shoot.count > 0), plot(all.shoot.count ~ Year))

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
Anova(shoot.count.glmer, test = "Chisq", type = 3)
overdisp_fun(shoot.count.glmer) # no overdispersion
plot(shoot.count.glmer) # looks okay

# calculate variance components for simplified model
shoot.count.up <- update(shoot.count.glmer, .~. -treatment*Year*genotype + (1|genotype) + treatment*Year)
Anova(shoot.count.up, test = "Chisq", type = 3)

summary(shoot.count.up)
var.calc(shoot.count.up)
plotREsim(REsim(shoot.count.up))
sem.model.fits(shoot.count.up) # not working....



# Likelihood ratio tests
w.sh.3 <- drop1(shoot.count.glmer, test = "Chisq") %>% tidy()
w.sh.2 <- drop1(update(shoot.count.glmer, .~. -treatment:genotype:Year), test = "Chisq") %>% tidy()
w.sh.1 <- drop1(update(shoot.count.glmer, .~. -treatment*genotype*Year + treatment + genotype + Year), test = "Chisq") %>% tidy()

shoot.count.anova <- bind_rows(w.sh.1, w.sh.2, w.sh.3) %>%
  filter(term != "<none>") %>%
  mutate(Response = "shoot count",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "Chi-square",
         Experiment = "wind",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2 for significant predictors
shoot.count.2x <- update(shoot.count.glmer, .~. -treatment:genotype:Year)
shoot.count.main <- update(shoot.count.glmer, .~. -treatment*genotype*Year + treatment + genotype + Year)

shoot.count.R2 <- 
  data.frame(Factor = c("genotype:Year","treatment:Year","treatment","genotype","Year"),
             Response = "shoot count",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(shoot.count.2x, update(shoot.count.2x, .~. -genotype:Year)),
                        deltaR2(shoot.count.2x, update(shoot.count.2x, .~. -treatment:Year)),
                        deltaR2(shoot.count.main, update(shoot.count.main, .~. -treatment)),
                        deltaR2(shoot.count.main, update(shoot.count.main, .~. -genotype)),
                        deltaR2(shoot.count.main, update(shoot.count.main, .~. -Year))),2))

## Wind: average shoot length analysis ----

# plot 2-way interactions
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(treatment, genotype, all.shoot.avg.length))
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(treatment, Year, all.shoot.avg.length))
with(filter(wind.df, all.shoot.avg.length > 0), interaction.plot(Year, genotype, all.shoot.avg.length))

# plot main effects
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ Year))
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ genotype))
with(filter(wind.df, all.shoot.avg.length > 0), plot(all.shoot.avg.length ~ treatment))

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
Anova(shoot.length.lmer, test = "F", type = 3) %>% tidy()

# Kenward-Roger tests 
shoot.length.anova <- anova(shoot.length.lmer, ddf = "Kenward-Roger")   %>%
  tidy() %>%
  mutate(Response = "average shoot length",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "F (Kenward-Roger)",
         Experiment = "wind",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "log") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2 for significant predictors
# no interaction effects were significant so I didn't calculate R2
shoot.length.main <- update(shoot.length.lmer, .~. -treatment*genotype*Year + treatment + genotype + Year)

shoot.length.R2 <- 
  data.frame(Factor = c("treatment","genotype","Year"),
             Response = "average shoot length",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(shoot.length.main, update(shoot.length.main, .~. -treatment)),
                                deltaR2(shoot.length.main, update(shoot.length.main, .~. -genotype)),
                                deltaR2(shoot.length.main, update(shoot.length.main, .~. -Year))),2))

## Wind: leaf water content ----
hist(wind.df$leaf_WC)
WC_nonsense <- which(wind.df$leaf_WC < 0)

plot(leaf_WC ~ genotype, wind.df[-WC_nonsense, ])
with(filter(wind.df[-WC_nonsense, ], leaf_WC > 0), interaction.plot(Year, genotype, leaf_WC))

leaf_WC.lmer <- lmer(log(leaf_WC) ~ treatment*genotype*Year +
                            (1|block) + 
                            (1|block:treatment) +
                            (1|plant.id),
                         data = wind.df[-WC_nonsense, ],
                         contrast = list(treatment = "contr.sum",
                                         genotype = "contr.sum",
                                         Year = "contr.sum"))
summary(leaf_WC.lmer)
plot(leaf_WC.lmer) # looks good

# Kenward-Roger test
leaf_WC.anova <- anova(leaf_WC.lmer, ddf = "Kenward-Roger")  %>%
  tidy() %>%
  mutate(Response = "leaf water content",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "F (Kenward-Roger)",
         Experiment = "wind",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "log") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
leaf_WC.2x <- update(leaf_WC.lmer, .~. -treatment:genotype:Year)
leaf_WC.main <- update(leaf_WC.lmer, .~. -treatment*genotype*Year + treatment + genotype + Year)

leaf_WC.R2 <- 
  data.frame(Factor = c("genotype:Year","genotype"),
             Response = "leaf water content",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(leaf_WC.2x, update(leaf_WC.2x, .~. -genotype:Year)),
                                deltaR2(leaf_WC.main, update(leaf_WC.main, .~. -genotype))),2))

## Wind: trichome density analysis ----
with(filter(wind.df, leaf_trichome.density != "NA"), interaction.plot(treatment, genotype, leaf_trichome.density))

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
w.td.2 <- drop1(leaf_trichome.density.glmer, 
      test = "Chisq") %>% tidy()
w.td.1 <- drop1(update(leaf_trichome.density.glmer, .~. -treatment:genotype), 
      test = "Chisq") %>% tidy()

trichome.density.anova <- bind_rows(w.td.1, w.td.2) %>%
  filter(term != "<none>") %>%
  mutate(Response = "leaf trichome density",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "Chi-square",
         Experiment = "wind",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2
leaf_trichome.density.main <- update(leaf_trichome.density.glmer, .~. -treatment*genotype + treatment + genotype)

leaf_trichome.density.R2 <- 
  data.frame(Factor = "genotype",
             Response = "leaf trichome density",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(leaf_trichome.density.main, 
                                        update(leaf_trichome.density.main, .~. -genotype))),2))

## Wind: leaf C:N analysis ----
hist(wind.df$leaf_C_N)
plot(leaf_C_N ~ genotype, wind.df)

leaf_CN.lmer <- lmer(log(leaf_C_N) ~ treatment*genotype +
                            (1|block) + (1|block:treatment),
                          data = wind.df,
                     contrasts = list(treatment = "contr.sum",
                                      genotype = "contr.sum"))
summary(leaf_CN.lmer)
plot(leaf_CN.lmer) # looks pretty good

## Kenwar-Roger test
leaf_CN.anova <- anova(leaf_CN.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "leaf C:N",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "F (Kenward-Roger)",
         Experiment = "wind",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "log") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
leaf_CN.main <- update(leaf_CN.lmer, .~. -treatment*genotype + treatment + genotype)

leaf_CN.R2 <- 
  data.frame(Factor = c("genotype"),
             Response = "leaf C:N",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(leaf_CN.main, update(leaf_CN.main, .~. -genotype))),2))

## Wind: SLA analysis ----
hist(wind.df$SLA)
SLA.lmer <- lmer(SLA ~ treatment*genotype +
                            (1|block) + (1|block:treatment),
                          data = wind.df,
                      contrasts = list(treatment = "contr.sum",
                                       genotype = "contr.sum"))
summary(SLA.lmer)
plot(SLA.lmer) # residuals look okay

## Kenward-Roger test
SLA.anova <- anova(SLA.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "SLA",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "F (Kenward-Roger)",
         Experiment = "wind",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
SLA.main <- update(SLA.lmer, .~. -treatment*genotype + treatment + genotype)

SLA.R2 <- 
  data.frame(Factor = c("genotype"),
             Response = "SLA",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "wind",
             delta_R2 = round(c(deltaR2(SLA.main, update(SLA.main, .~. -genotype))),2))

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
aa.df <- mutate(aa.df, scale.AntDist = scale(Ant.Mound.Dist))
aa.height.lmer <- mixed(Height ~ Aphid.Treatment*scale.AntDist*Genotype +
                           (1|Block) + (1|Block:fact.Ant.Mound.Dist),
                         data = filter(aa.df, Year == "2012"),
                         contrasts = list(Aphid.Treatment = "contr.sum",
                                          Genotype = "contr.sum"), )
library(afex)
mixed(aa.height.lmer)
plot(aa.height.lmer) # not great...
summary(aa.height.lmer)

Anova(aa.height.lmer, test = "F", type = 3)

aa.height.simp <- update(aa.height.lmer, .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Genotype))
var.calc(aa.height.simp)

## Kenward-Roger test
aa.height.anova <- anova(aa.height.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "plant height",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "F (Kenward-Roger)",
         Experiment = "ant-aphid",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
aa.height.main <- update(aa.height.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + 
                           Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype)

aa.height.R2 <- 
  data.frame(Factor = "Genotype",
             Response = "plant height",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "ant-aphid",
             delta_R2 = round(c(deltaR2(aa.height.main, update(aa.height.main, .~. -Genotype))),2))

## Ant-aphid: Shoot count analysis ----
ggplot(filter(aa.df, Year == "2012"),
       aes(x = Ant.Mound.Dist, y = log(all.shoot.count), color = Aphid.Treatment)) +
  geom_point() + 
  stat_smooth(method = "lm", family = "poisson", se = FALSE)

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

## Likelihood ratio tests
aa.sh.3 <- drop1(aa.shoot.count.glmer, test = "Chisq") %>% tidy()
aa.sh.2 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.Treatment:scale(Ant.Mound.Dist):Genotype), test = "Chisq") %>% tidy()
aa.sh.1 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype), test = "Chisq") %>% tidy()

aa.shoot.count.anova <- bind_rows(aa.sh.1, aa.sh.2, aa.sh.3) %>%
  filter(term != "<none>") %>%
  mutate(Response = "shoot count",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "Chi-square",
         Experiment = "ant-aphid",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2 for significant predictors
aa.shoot.count.2x <- update(aa.shoot.count.glmer, .~. -Aphid.Treatment:scale(Ant.Mound.Dist):Genotype)
aa.shoot.count.main <- update(aa.shoot.count.glmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + 
                           Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype)

aa.shoot.count.R2 <- 
  data.frame(Factor = c("Aphid.Treatment:scale(Ant.Mound.Dist)", "Aphid.Treatment", "Genotype"),
             Response = "shoot count",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "ant-aphid",
             delta_R2 = round(c(deltaR2(aa.shoot.count.2x, update(aa.shoot.count.2x, .~. -Aphid.Treatment:scale(Ant.Mound.Dist))),
                                deltaR2(aa.shoot.count.main, update(aa.shoot.count.main, .~. -Aphid.Treatment)),
                                deltaR2(aa.shoot.count.main, update(aa.shoot.count.main, .~. -Genotype))),2))

## Ant-Aphid: Average shoot length analysis ----
hist(aa.df$mature.shoot.avg.length)

aa.mature.shoot.avg.length.lmer <- lmer(mature.shoot.avg.length ~ Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + (1|Block) + (1|Block:Ant.Mound.Dist),
                                 data = filter(aa.df, Year == "2012"),
                                 contrasts = list(Aphid.Treatment = "contr.sum",
                                                  Genotype = "contr.sum"))
plot(aa.mature.shoot.avg.length.lmer) # not great...
summary(aa.mature.shoot.avg.length.lmer)

## Kenward-Roger
aa.mature.shoot.avg.length.anova <- anova(aa.mature.shoot.avg.length.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "average shoot length",
         Response.type = "phenotype",
         Specific.type = "plant growth",
         Test.type = "F (Kenward-Roger)",
         Experiment = "ant-aphid",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
aa.mature.shoot.avg.length.main <- update(aa.mature.shoot.avg.length.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist)*Genotype + 
                           Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype)

aa.mature.shoot.avg.length.R2 <- 
  data.frame(Factor = "Genotype",
             Response = "average shoot length",
             Response.type = "phenotype",
             Specific.type = "plant growth",
             Experiment = "ant-aphid",
             delta_R2 = round(c(deltaR2(aa.mature.shoot.avg.length.main, update(aa.mature.shoot.avg.length.main, .~. -Genotype))),2))

## Ant-Aphid: leaf water content analysis ----
hist(aa.df$leaf_WC)
aa.leaf_WC.lmer <- lmer(log(leaf_WC) ~ Aphid.Treatment*scale(Ant.Mound.Dist) + scale(Ant.Mound.Dist)*Genotype + (1|Block) + (1|Block:fact.Ant.Mound.Dist), 
                        data = filter(aa.df, Year == "2012"),
                        contrasts = list(Aphid.Treatment = "contr.sum",
                                         Genotype = "contr.sum"))
summary(aa.leaf_WC.lmer)
plot(aa.leaf_WC.lmer) # looks okay

## Kenward-Roger test
aa.leaf_WC.anova <- anova(aa.leaf_WC.lmer, ddf = "Kenward-Roger") %>%
  tidy() %>%
  mutate(Response = "leaf water content",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "F (Kenward-Roger)",
         Experiment = "ant-aphid",
         Model.type = "LMM",
         Error.dist = "gaussian",
         Response.trans = "log") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = NumDF, den_DF = DenDF, 
         Statistic = statistic, P_value = p.value, Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         P_value = round(P_value,3),
         den_DF = round(den_DF,1))

## Calculate R2
# Genotype was not significant, but I'm still going to calculate its R2.
aa.leaf_WC.main <- update(aa.leaf_WC.lmer, 
                         .~. -Aphid.Treatment*scale(Ant.Mound.Dist) - scale(Ant.Mound.Dist)*Genotype + 
                           Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype)

aa.leaf_WC.R2 <- 
  data.frame(Factor = "Genotype",
             Response = "leaf water content",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "ant-aphid",
             delta_R2 = round(c(deltaR2(aa.leaf_WC.main, update(aa.leaf_WC.main, .~. -Genotype))),2))

## Ant-aphid: leaf trichome density ----
ggplot(filter(aa.df, Year == "2012"),
       aes(x = Genotype, y = leaf_trichome.density)) +
  geom_boxplot()

aa.leaf_trichome.density.glmer <- glmer(leaf_trichome.density ~ Aphid.Treatment*scale(Ant.Mound.Dist) + scale(Ant.Mound.Dist)*Genotype + 
                                          (1|Block) + (1|Block:Ant.Mound.Dist) + (1|X),
                                        data = filter(aa.df, Year == "2012"),
                                        family = "poisson",
                                        contrasts = list(Aphid.Treatment = "contr.sum",
                                                         Genotype = "contr.sum"),
                                        control=glmerControl(optimizer="bobyqa",
                                                             optCtrl=list(maxfun=2e5)))
summary(aa.leaf_trichome.density.glmer)
plot(aa.leaf_trichome.density.glmer) # not great
overdisp_fun(aa.leaf_trichome.density.glmer)

## Likelihood ratio tests
aa.td.2 <- drop1(aa.leaf_trichome.density.glmer, test = "Chisq") %>% tidy()

aa.leaf_trichome.density.main <- update(aa.leaf_trichome.density.glmer, 
                                        .~. -Aphid.Treatment*scale(Ant.Mound.Dist) - scale(Ant.Mound.Dist)*Genotype + 
                                          Aphid.Treatment + scale(Ant.Mound.Dist) + Genotype)

aa.td.1 <- drop1(aa.leaf_trichome.density.main, test = "Chisq") %>% tidy()

aa.trichome.density.anova <- bind_rows(aa.td.1, aa.td.2) %>%
  filter(term != "<none>") %>%
  mutate(Response = "leaf trichome density",
         Response.type = "phenotype",
         Specific.type = "leaf quality",
         Test.type = "Chi-square",
         Experiment = "ant-aphid",
         Model.type = "GLMM",
         Error.dist = "poisson",
         Response.trans = "none") %>%
  select(Experiment, Response.type, Specific.type, Response,
         Factor = term, num_DF = df, #den_DF = DenDF, 
         Statistic = LRT, P_value = Pr.Chi., Model.type, 
         Test.type, Error.dist, Response.trans) %>%
  mutate(Statistic = round(Statistic,2),
         #den_DF = round(den_DF,1),
         P_value = round(P_value,3))

## Calculate R2
aa.leaf_trichome.density.R2 <- 
  data.frame(Factor = "Genotype",
             Response = "leaf trichome.density",
             Response.type = "phenotype",
             Specific.type = "leaf quality",
             Experiment = "ant-aphid",
             delta_R2 = round(c(deltaR2(aa.leaf_trichome.density.main, update(aa.leaf_trichome.density.main, .~. -Genotype))),2))

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
                          aa.trichome.density.anova)
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
#write.csv(trait.anovas, "trait.anovas.csv")

ggplot(trait.R2s, aes(x = Response, y = delta_R2, fill = Factor)) + geom_bar(stat = "identity", position = position_dodge()) +
  #coord_flip() + 
  facet_grid(Specific.type~Experiment) 

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

