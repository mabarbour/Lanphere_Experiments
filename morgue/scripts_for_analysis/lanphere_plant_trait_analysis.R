#### ANALYZE PLANT TRAIT DATA IN LANPHERE GXE EXPERIMENT 

## Load required libraries and fucnctions ----
source('scripts_for_analysis/required_libraries.R')

#### Load required data sets ----

## wind plant traits
w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block), 
         Year = as.factor(Year),
         X = as.factor(X)) # observation-level random effect to account for overdispersion
glimpse(w.trait.df)

## examine the relationship between C:N data collected in the common garden and that collected in the wind experiment.
# run 'garden_chemistry.R' first to get 'chem.data'
C_N_test <- w.trait.df %>% group_by(Genotype, Wind.Exposure) %>% summarise(avg.leaf_CN = mean(leaf_C_N, na.rm = TRUE)) %>%
  left_join(., chem.data)

plot(avg.leaf_CN ~ C_N_imputed, filter(C_N_test, Wind.Exposure == "Exposed"))
plot(avg.leaf_CN ~ C_N_imputed, filter(C_N_test, Wind.Exposure == "Unexposed"))
with(filter(C_N_test, Wind.Exposure == "Exposed"), cor.test(avg.leaf_CN, C_N_imputed))
with(filter(C_N_test, Wind.Exposure == "Unexposed"), cor.test(avg.leaf_CN, C_N_imputed))

## ant-aphid above-ground plant traits
aa.trait.df <- read.csv('final_data/ant_aphid_trait_df.csv') %>%
  tbl_df() %>%
  filter(Year == "2012") %>%
  mutate(fact.Ant.mound.dist = as.factor(Ant.mound.dist), 
         Block = as.factor(Block),
         X = as.factor(X), # observation-level random effect to account for overdispersion.
         Plot_code = paste(Block, Ant.mound.dist, sep = "_")) %>%
  select(-Year)
glimpse(aa.trait.df)


#### Examine correlation among plant traits ----
## Wind: phenotypic correlations
# 2012
scatterplotMatrix( ~ Height + all.shoot.avg.length + all.shoot.count + leaf_trichome.density + leaf_WC, data = filter(w.trait.df, Year == "2012"))
corr.test(x = select(filter(w.trait.df, Year == "2012", Wind.Exposure == "Exposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))
corr.test(x = select(filter(w.trait.df, Year == "2012", Wind.Exposure == "Unexposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))

# 2013
scatterplotMatrix( ~ Height + all.shoot.avg.length + all.shoot.count + leaf_WC + leaf_C_N + SLA, data = filter(w.trait.df, Year == "2013"))
corr.test(x = select(filter(w.trait.df, Year == "2013", Wind.Exposure == "Exposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_WC, SLA, leaf_C_N))
corr.test(x = select(filter(w.trait.df, Year == "2013", Wind.Exposure == "Unexposed"), Height, all.shoot.avg.length, all.shoot.count, leaf_WC, SLA, leaf_C_N))

## Ant-aphid: phenotypic correlations
scatterplotMatrix( ~ Height + mature.shoot.avg.length + all.shoot.count + leaf_trichome.density + leaf_WC, data = filter(aa.trait.df, Year == "2012"))
corr.test(x = select(filter(aa.trait.df, Year == "2012"), Height, mature.shoot.avg.length, all.shoot.count, leaf_trichome.density, leaf_WC))

## Wind: plant height analysis ----

# LMMM 
height.lmer <- lmer(Height ~ Wind.Exposure*Year*Genotype +  
                        (1|Block) + 
                        (1|Block:Wind.Exposure) +
                        (1|plant_ID),
                    data = w.trait.df,
                    contrasts = list(Wind.Exposure = "contr.sum",
                                     Year = "contr.sum",
                                     Genotype = "contr.sum"))
summary(height.lmer)

# Anova table
(height.anova <- anova.table(height.lmer, type = 2, experiment = "wind"))

# Effects
Effect(c("Wind.Exposure","Genotype","Year"), height.lmer)
plot(Effect(c("Wind.Exposure","Genotype","Year"), height.lmer))
Effect(c("Wind.Exposure","Year"), height.lmer) # effect of unexposed plots increased from 1.2-fold in 2012 to 2-fold in 2013.
Effect(c("Year"), height.lmer) # plants were 39% shorter in 2013 vs. 2012. 
Effect(c("Wind.Exposure"), height.lmer) # plants were 29% shorter in unexposed plots
Effect(c("Genotype"), height.lmer) # plants varied 2-fold in height among the most disparate Genotypes.

w.height.resp <- as.data.frame(Effect(c("Genotype"), height.lmer)) %>% bind_rows(., as.data.frame(Effect(c("Wind.Exposure","Year"), height.lmer))) %>% mutate(response = "height")

# Calculate R2
height.up <- update(height.lmer, .~. -Wind.Exposure*Year*Genotype + (1|Genotype) + Wind.Exposure*Year) # simplify model

(height.R2 <- var.table(height.up, experiment = "wind"))

## Wind: shoot count analysis ----

# GLMM
shoot.count.glmer <- glmer(all.shoot.count ~ Wind.Exposure*Genotype*Year + 
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           (1|plant_ID),
                         data = w.trait.df,
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(shoot.count.glmer)
overdisp_fun(shoot.count.glmer) # no overdispersion

## Likelihood-ratio tests
(shoot.count.3 <- drop1(shoot.count.glmer, test = "Chisq"))
(shoot.count.2 <- drop1(update(shoot.count.glmer, .~. -Wind.Exposure:Genotype:Year), test = "Chisq"))
(shoot.count.1 <- drop1(update(shoot.count.glmer, .~. -Wind.Exposure:Genotype:Year -Wind.Exposure:Genotype -Wind.Exposure:Year -Genotype:Year), test = "Chisq"))
#(shoot.count.anova <- anova.table(shoot.count.glmer, test = "Chisq", type = 2, experiment = "wind"))

## Calculate effects
Effect(c("Wind.Exposure","Year"), shoot.count.glmer) # effect of unexposed plots increased from 1.1-fold in 2012 to 1.6-fold in 2013.
Effect(c("Year"), shoot.count.glmer)
plot(Effect(c("Year","Genotype"), shoot.count.glmer) )
Effect(c("Wind.Exposure"), shoot.count.glmer) # plants produced 22% fewer shoots.
Effect(c("Genotype"), shoot.count.glmer) # plants varied nearly 2.3-fold in shoot count among the most disparate Genotypes.
shoot.count.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), shoot.count.glmer)) %>% mutate(response = "shoot.count")

shoot.count.GxEy <- as.data.frame(Effect(c("Year","Genotype"), shoot.count.glmer)) %>% mutate(response = "shoot.count")
shoot.count.GxEy %>% ggplot(aes(x = Year, y = fit, color = Genotype, group = Genotype)) + geom_line()

# Calculate variance components for simplified model
shoot.count.up <- update(shoot.count.glmer, .~. -Wind.Exposure*Year*Genotype + (Year|Genotype) + Wind.Exposure*Year)

(shoot.count.R2 <- var.table(shoot.count.up, experiment = "wind"))

## Wind: average shoot length analysis ----

# LMM. 
shoot.length.lmer <- lmer(log(all.shoot.avg.length) ~ Wind.Exposure*Genotype*Year +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|plant_ID),
                          data = w.trait.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"))
summary(shoot.length.lmer)

# Kenward-Roger tests 
(shoot.length.anova <- anova.table(shoot.length.lmer, type = 2, experiment = "wind"))

# effects
Effect(c("Wind.Exposure"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants shoots grew 28% shorter.
Effect(c("Year"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants shoots grew 1.7-fold longer in 2012.
Effect(c("Genotype"), shoot.length.lmer, transformation = list(link = log, inverse = exp)) # plants varied 2.2-fold in shoot length among the most disparate Genotypes.

plot(Effect(c("Wind.Exposure","Genotype"), shoot.length.lmer, transformation = list(link = log, inverse = exp)))

## Calculate R2 for significant predictors
# no interaction effects were significant so I didn't calculate R2
shoot.length.up <- update(shoot.length.lmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

(shoot.length.R2 <- var.table(shoot.length.up, experiment = "wind"))

## Wind: leaf water content ----

leaf_WC.lmer <- lmer(log(leaf_WC) ~ Wind.Exposure*Genotype*Year + #scale(avg.moisture.vwc) +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|plant_ID),
                         data = w.trait.df,
                         contrast = list(Wind.Exposure = "contr.sum",
                                         Genotype = "contr.sum",
                                         Year = "contr.sum"))
summary(leaf_WC.lmer)  # note that plant_ID explains the least amount of variance

# Kenward-Roger test
(leaf_WC.anova <- anova.table(leaf_WC.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype","Year"), leaf_WC.lmer, transformation = list(link = log, inverse = exp)))
plot(Effect(c("Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp)))

WC.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_WC")

## Calculate R2
leaf_WC.up <- update(leaf_WC.lmer, .~. -Wind.Exposure*Genotype*Year + Year + (0+Year|Genotype))

(leaf_WC.R2 <- var.table(leaf_WC.up, experiment = "wind"))

## Wind: trichome density analysis ----
leaf_trichome.density.glmer <- glmer(
  leaf_trichome.density ~ Wind.Exposure*Genotype +
    (1|Block) + 
    (1|Block:Wind.Exposure) +
    (1|plant_ID), # need individual-level random effect to account for overdispersion
  data = w.trait.df, 
  family = "poisson",
  contrasts = list(Wind.Exposure = "contr.sum",
                   Genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(leaf_trichome.density.glmer)
overdisp_fun(leaf_trichome.density.glmer)
plot(leaf_trichome.density.glmer) 

# Likelihood-ratio tests
(td.2 <- drop1(leaf_trichome.density.glmer, test = "Chisq"))
(td.1 <- drop1(update(leaf_trichome.density.glmer, .~. -Wind.Exposure:Genotype), test = "Chisq"))
#trichome.density.anova <- anova.table(leaf_trichome.density.glmer, test = "Chisq", experiment = "wind")

# effects
plot(Effect("Genotype", leaf_trichome.density.glmer)) # willows varied 46-fold in leaf trichome density.

td.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_trichome.density.glmer)) %>% mutate(response = "leaf_trichome.density")

## Calculate R2
leaf_trichome.density.up <- update(leaf_trichome.density.glmer, .~. -Wind.Exposure*Genotype + (1|Genotype)) 

(leaf_trichome.density.R2 <- var.table(leaf_trichome.density.up, experiment = "wind"))

## Wind: leaf C:N analysis ----

# LMM
leaf_CN.lmer <- lmer(log(leaf_C_N) ~ Wind.Exposure*Genotype +
                            (1|Block) + (1|Block:Wind.Exposure),
                          data = w.trait.df,
                     contrasts = list(Wind.Exposure = "contr.sum",
                                      Genotype = "contr.sum"))
print(summary(leaf_CN.lmer), correlation = TRUE)
plot(leaf_CN.lmer) # looks pretty good

## Kenwar-Roger test
(leaf_CN.anova <- anova.table(leaf_CN.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype"), leaf_CN.lmer, transformation = list(link = log, inverse = exp))) # willows varied 1.6-fold in leaf C:N

CN.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), leaf_CN.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_CN")

## Calculate R2
leaf_CN.up <- update(leaf_CN.lmer, .~. -Wind.Exposure*Genotype + (1|Genotype)) 

(leaf_CN.R2 <- var.table(leaf_CN.up, experiment = "wind"))

## Wind: SLA analysis ----
hist(w.trait.df$SLA)

# LMM
SLA.lmer <- lmer(SLA ~ Wind.Exposure*Genotype +
                            (1|Block) + (1|Block:Wind.Exposure),
                          data = w.trait.df,
                      contrasts = list(Wind.Exposure = "contr.sum",
                                       Genotype = "contr.sum"))
summary(SLA.lmer)

## Kenward-Roger test
(SLA.anova <- anova.table(SLA.lmer, type = 2, experiment = "wind"))

# effects
plot(Effect(c("Genotype"), SLA.lmer)) # Genotypes varied 1.5-fold in SLA.

SLA.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), SLA.lmer)) %>% mutate(response = "SLA")

# Calculate R2
SLA.up <- update(SLA.lmer, .~. -Wind.Exposure*Genotype + (1|Genotype))

(SLA.R2 <- var.table(SLA.up, experiment = "wind"))

## Wind: Root C:N analysis ----
rootCN.lmer <- lmer(log(root_CN) ~ Wind.Exposure*Genotype + (1|Block) + (1|Block:Wind.Exposure),
                    data = w.trait.df, # results are qualitatively the same if I filter C:N less than 100
                    contrasts = list(Wind.Exposure = "contr.sum",
                                     Genotype = "contr.sum")) # restricting to biologically reasonable values
summary(rootCN.lmer)
plot(rootCN.lmer)

(rootCN.anova <- anova.table(rootCN.lmer, type = 2, experiment = "wind"))

rootCN.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), rootCN.lmer)) %>% mutate(response = "root_CN")

rootCN.up <- update(rootCN.lmer, .~. -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))

(rootCN.R2 <- var.table(rootCN.up, experiment = "wind"))


## Ant-aphid: plant height analysis ----

# LMM
aa.height.lmer <- lmer(Height ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + 
                           (1|Block) + (1|Block:fact.Ant.mound.dist),
                         data = aa.trait.df, #filter(aa.df, Year == "2012"),
                         contrasts = list(Aphid.treatment = "contr.sum",
                                          Genotype = "contr.sum"))
summary(aa.height.lmer)

## Kenward-Roger test
(aa.height.anova <- anova.table(aa.height.lmer, type = 2, experiment = "ant-aphid"))

# effects
plot(Effect("Genotype", aa.height.lmer)) # willows varied 2-fold in height among the most disparate genotypes.

aa.height.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.height.lmer)) %>% mutate(response = "height")


## Calculate R2
aa.height.up <- update(aa.height.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.height.R2 <- var.table(aa.height.up, experiment = "ant-aphid"))

## Ant-aphid: Shoot count analysis ----

# GLMM
aa.shoot.count.glmer <- glmer(
  all.shoot.count ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + 
    (1|Block) + (1|Block:fact.Ant.mound.dist),
                             data = aa.trait.df, #filter(aa.df, Year == "2012"),
                             family = "poisson",
                             contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
  control=glmerControl(optimizer="bobyqa",
                       optCtrl=list(maxfun=2e5)))
summary(aa.shoot.count.glmer)
overdisp_fun(aa.shoot.count.glmer) # no overdispersion

## Likelihood rato tests
(aa.shoot.count.3 <- drop1(aa.shoot.count.glmer, test = "Chisq"))
(aa.shoot.count.2 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.treatment:scale(Ant.mound.dist):Genotype), test = "Chisq"))
(aa.shoot.count.1 <- drop1(update(aa.shoot.count.glmer, .~. -Aphid.treatment:scale(Ant.mound.dist):Genotype -Aphid.treatment:scale(Ant.mound.dist) -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq"))

# effects
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.shoot.count.glmer)) # don't quite understand these effects
plot(Effect(c("Aphid.treatment"), aa.shoot.count.glmer))
Effect(c("Genotype"), aa.shoot.count.glmer) # 1.8-fold variation

aa.shoot.count.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.shoot.count.glmer)) %>% mutate(response = "shoot.count")
aa.shoot.count.ExE <- as.data.frame(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.shoot.count.glmer)) %>% mutate(response = "shoot.count")

## Calculate R2 for significant predictors
aa.shoot.count.up <- update(aa.shoot.count.glmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + Aphid.treatment*scale(Ant.mound.dist) + (1|Genotype))

(aa.shoot.count.R2 <- var.table(aa.shoot.count.up, experiment = "ant-aphid"))

## Ant-aphid: Average shoot length analysis ----
hist(aa.df$mature.shoot.avg.length)

aa.mature.shoot.avg.length.lmer <- lmer(mature.shoot.avg.length ~ Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Block) + (1|Block:Ant.mound.dist),
                                 data = aa.trait.df, #filter(aa.df, Year == "2012"),
                                 contrasts = list(Aphid.treatment = "contr.sum",
                                                  Genotype = "contr.sum"))
summary(aa.mature.shoot.avg.length.lmer)

## Kenward-Roger
(aa.mature.shoot.avg.length.anova <- anova.table(aa.mature.shoot.avg.length.lmer, type = 2, experiment = "ant-aphid"))

## effects
Effect("Genotype", aa.mature.shoot.avg.length.lmer) # 2.5-fold

aa.shoot.length.GxE <- as.data.frame((Effect(c("Aphid.treatment","Genotype"), aa.mature.shoot.avg.length.lmer))) %>% mutate(response = "shoot.length")

## Calculate R2
aa.mature.shoot.avg.length.up <- update(aa.mature.shoot.avg.length.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.mature.shoot.avg.length.R2 <- var.table(aa.mature.shoot.avg.length.up, experiment = "ant-aphid"))

## Ant-Aphid: leaf water content analysis ----
hist(aa.df$leaf_WC)

## LMM. Had to fix a less complex model because we didn't have enough replication of all treatment combinations
aa.leaf_WC.lmer <- lmer(log(leaf_WC) ~ Aphid.treatment*Genotype*scale(Ant.mound.dist) + #Aphid.treatment*scale(Ant.mound.dist) + 
                          #scale(Ant.mound.dist)*Genotype + 
                          (1|Block) + (1|Block:fact.Ant.mound.dist), 
                        data = aa.trait.df, #filter(aa.df, Year == "2012"),
                        contrasts = list(Aphid.treatment = "contr.sum",
                                         Genotype = "contr.sum"))
summary(aa.leaf_WC.lmer)

## Kenward-Roger test
(aa.leaf_WC.anova <- anova.table(aa.leaf_WC.lmer, type = 2, experiment = "ant-aphid")) # nothing significant, but genotype is the closest.

## effects
plot(Effect("Genotype",aa.leaf_WC.lmer, transformation = list(link = log, inverse = exp))) # no significant effect

aa.WC.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"),aa.leaf_WC.lmer, transformation = list(link = log, inverse = exp))) %>% mutate(response = "leaf_WC")

## Calculate R2
# Genotype was not significant, but I'm still going to calculate its R2.
aa.leaf_WC.up <- update(aa.leaf_WC.lmer, 
                         .~. -Aphid.treatment*scale(Ant.mound.dist) - scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.leaf_WC.R2 <- var.table(aa.leaf_WC.lmer, experiment = "ant-aphid"))

## Ant-aphid: leaf trichome density ----

# GLMM. Again, needed to simplify due to lack of replication
aa.leaf_trichome.density.glmer <- glmer(leaf_trichome.density ~ (Aphid.treatment + Genotype + scale(Ant.mound.dist))^2 +
                                          (1|Block) + 
                                          (1|X) +
                                          (1|Block:Ant.mound.dist),
                                        data = aa.trait.df, 
                                        family = "poisson",
                                        contrasts = list(Aphid.treatment = "contr.sum",
                                                         Genotype = "contr.sum"),
                                        control=glmerControl(optimizer="bobyqa",
                                                             optCtrl=list(maxfun=2e5)))
summary(aa.leaf_trichome.density.glmer)
overdisp_fun(aa.leaf_trichome.density.glmer) # overdispersed so I modelled an individual-level random effect

## Likelihood ratio tests
(aa.td.2 <- drop1(aa.leaf_trichome.density.glmer, test = "Chisq"))
(aa.td.1 <- drop1(update(aa.leaf_trichome.density.glmer, .~. -(Aphid.treatment + Genotype + scale(Ant.mound.dist))^2 + Aphid.treatment + Genotype + scale(Ant.mound.dist)), test = "Chisq"))

with(filter(aa.trait.df, leaf_trichome.density != "NA"), interaction.plot(x.factor = Aphid.treatment, response = leaf_trichome.density, trace.factor = Genotype))
table(filter(aa.trait.df, leaf_trichome.density != "NA")$Genotype, filter(aa.trait.df, leaf_trichome.density != "NA")$Aphid.treatment)

## effects
plot(Effect(c("Genotype","Aphid.treatment"), aa.leaf_trichome.density.glmer))
plot(Effect("Genotype",aa.leaf_trichome.density.glmer)) # 30-fold variation in leaf trichome density

aa.td.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"),aa.leaf_trichome.density.glmer)) %>% mutate(response = "leaf_trichome.density")

## Calculate R2
aa.leaf_trichome.density.up <- update(aa.leaf_trichome.density.glmer, .~. -Aphid.treatment*scale(Ant.mound.dist) - scale(Ant.mound.dist)*Genotype + (1|Genotype)) 

(aa.leaf_trichome.density.R2 <- var.table(aa.leaf_trichome.density.glmer, experiment = "ant-aphid"))

