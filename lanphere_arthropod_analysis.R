## load required libraries ----
library(dplyr)
library(reshape2)
library(ggplot2)
library(psych)
library(pbkrtest) # for some reason, I have to hve pbkrtest loaded with lmerTest for it to run the Kenward-Roger test appropriately.
library(lmerTest)
#library(RLRsim)
library(car)
library(vegan)
library(mvabund)

## upload datasets
aa.aphid.GR <- read.csv("~/Documents/Lanphere_Experiments/final_data/ant_aphid_Aphis_popgrowth_df.csv") %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  select(Date_rel:plant_code, Aphis.growth.rate)

glimpse(aa.aphid.GR)

aa.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/ant_aphid_arthropod_df.csv')

wind.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_arthropod_df.csv')



## exploratory plots ----

## wind arthropod community
wind.arth.nonzero <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Spider))
wind.arth.rich <- rowSums(
  select(wind.arth.df, Gracilliaridae_miner:Spider) > 0)

wind.arth.df <- mutate(wind.arth.df, wind.arth.nonzero = wind.arth.nonzero,
                       wind.arth.rich = wind.arth.rich)

# 2012
wind.arth.2012 <- wind.arth.df %>%
  select(Year:plant_code, Gracilliaridae_miner:Spider, wind.arth.rich, wind.arth.nonzero) %>%
  filter(Year == 2012, wind.arth.nonzero > 0, Genotype != "J")

wind.comm.2012 <- wind.arth.2012 %>%
  select(Gracilliaridae_miner:Spider) 
colSums(wind.comm.2012)/sum(colSums(wind.comm.2012)) # dominant groups: aphids, leaf miners, leaf tiers, and spiders

wind.2012.adonis <- adonis(wind.comm.2012 ~ Genotype*Wind.Exposure, data = wind.arth.2012, strata = wind.arth.2012$Block, method = "bray")
wind.2012.adonis # consider using euclidean distances to preserve all of the samples in the dataset.

wind.12.cap <- capscale(wind.comm.2012 ~ Genotype + Wind.Exposure + Condition(Block), data = wind.arth.2012, distance = "bray")
summary(wind.12.cap)
plot(wind.12.cap, display = c("cn","sp"))

anova(betadisper(vegdist(wind.comm.2012, method = "bray"), 
                 group = wind.arth.2012$Genotype,
                 bias.adjust = TRUE))
boxplot(betadisper(vegdist(wind.comm.2012, method = "bray"), 
                 group = wind.arth.2012$Genotype,
                 bias.adjust = TRUE))

wind.2012.rda <- rda(decostand(wind.comm.2012,
                               method = "hellinger") ~ 
                       Wind.Exposure +
                       Condition(Block),
                     data = wind.arth.2012)
summary(wind.2012.rda)
plot(wind.2012.rda, display = c("cn","sp"))

# note that using wind.arth.2012 greatly reduces the size of the dataset.
wind.arth.2012.lmer <- wind.arth.df %>%
  select(Year:plant_code, Gracilliaridae_miner:Spider, wind.arth.nonzero, wind.arth.rich) %>%
  filter(Year == 2012)

rich.12 <- glmer((wind.arth.rich>0) ~ Wind.Exposure*Genotype + 
                  (1|Block/Wind.Exposure),
                family = "binomial",
                wind.arth.2012.lmer)
summary(rich.12)
anova(rich.12, ddf = "Kenward-Roger")
plot(rich.12)
plot((wind.arth.rich) ~ Wind.Exposure + Genotype, wind.arth.2012.lmer)

arth.abund.12 <- lmer(log(wind.arth.nonzero+1) ~ Wind.Exposure*Genotype + 
                        (1|Block/Wind.Exposure),
                      wind.arth.2012.lmer)
anova(arth.abund.12, ddf = "Kenward-Roger")
plot(arth.abund.12)
plot(wind.arth.nonzero ~ Wind.Exposure, wind.arth.2012.lmer)

hist(log(wind.arth.2012$Gracilliaridae_miner+1))
w.Grac.12 <- lmer(log(Gracilliaridae_miner+1) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.Grac.12)
plot(w.Grac.12)
anova(w.Grac.12, type = 3, ddf = "Kenward-Roger")
plot(log(Gracilliaridae_miner+1) ~ Genotype, wind.arth.2012)
table(wind.arth.2012$Genotype)
plot(log(Gracilliaridae_miner+1) ~ Wind.Exposure, wind.arth.2012)

w.Grac.12.glmer <- glmer((Gracilliaridae_miner>0) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure), family = "binomial",
                  wind.arth.2012.lmer)
summary(w.Grac.12.glmer)
plot(w.Grac.12.glmer)
anova(w.Grac.12.glmer, update(w.Grac.12.glmer, .~. -Wind.Exposure))
anova(w.Grac.12, type = 3, ddf = "Kenward-Roger")

w.spid.12 <- lmer(log(Spider+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.spid.12)
anova(w.spid.12, type = 3, ddf = "Kenward-Roger")
plot(w.spid.12)

w.Aphid.12 <- lmer(log(Aphididae+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2012)
summary(w.Aphid.12)
anova(w.Aphid.12, type = 3, ddf = "Kenward-Roger")
plot(w.Aphid.12)
plot(log(Aphididae+1) ~ Wind.Exposure, wind.arth.2012)

w.Aphid.12 <- glmer(Aphididae ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                    family = "poisson",
                   wind.arth.2012.lmer)
summary(w.Aphid.12)
anova(w.Aphid.12, type = 3, ddf = "Kenward-Roger")
plot(w.Aphid.12)

w.Tort.12 <- lmer(log(Tortricidiae_leaftier+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                  wind.arth.2012)
summary(w.Tort.12)
anova(w.Tort.12, type = 3, ddf = "Kenward-Roger")
plot(w.Tort.12)

anova(betadisper(vegdist(wind.comm.2012), 
                 group = wind.arth.2012$Wind.Exposure,
                 bias.adjust = TRUE))
anova(betadisper(vegdist(wind.comm.2012), 
                 group = wind.arth.2012$Genotype,
                 bias.adjust = TRUE)) # violates assumption of adonis...

# 2013
wind.arth.2013 <- wind.arth.df %>%
  select(Year:plant_code, Gracilliaridae_miner:Spider) %>%
  filter(Year == 2013, wind.arth.nonzero > 0)

wind.comm.2013 <- wind.arth.2013 %>%
  select(Gracilliaridae_miner:Spider) 
colSums(wind.comm.2013)/sum(colSums(wind.comm.2013)) # dominant groups: leaf miners, leaf tiers, galls, and spiders

wind.2013.adonis <- adonis(wind.comm.2013 ~ Genotype*Wind.Exposure, data = wind.arth.2013, strata = wind.arth.2013$Block)
wind.2013.adonis # interesting... only an effect of wind exposure

w.spid.13 <- lmer(log(Spider+1) ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.spid.13)
anova(w.spid.13, type = 3, ddf = "Kenward-Roger")
plot(w.spid.13)

hist(wind.arth.2013$Cecidomyiidae_gall)
w.Cecid.13 <- lmer(log(Cecidomyiidae_gall+1) ~ Wind.Exposure + Genotype + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.Cecid.13)
anova(w.Cecid.13, type = 3, ddf = "Kenward-Roger")
plot(w.Cecid.13)
plot(log(Cecidomyiidae_gall+1) ~ Wind.Exposure + Genotype, wind.arth.2013)

w.Cecid.13.glmer <- glmer((Cecidomyiidae_gall>0) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure), family = "binomial",
                   wind.arth.2013)
summary(w.Cecid.13.glmer)

hist(log(wind.arth.2013$Tortricidiae_leaftier+1))
w.Tort.13 <- lmer(log(Tortricidiae_leaftier+1) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure),
                   wind.arth.2013)
summary(w.Tort.13)
anova(w.Tort.13, type = 3, ddf = "Kenward-Roger")
plot(w.Tort.13)
plot(log(Tortricidiae_leaftier+1) ~ Wind.Exposure, wind.arth.2013)

w.Tort.13.glmer <- glmer((Tortricidiae_leaftier>0) ~ Wind.Exposure + (1|Genotype) + (1|Block/Wind.Exposure), family = "binomial",
                  wind.arth.2013)
summary(w.Tort.13)
anova(w.Tort.13, type = 3, ddf = "Kenward-Roger")

plot(betadisper(vegdist(wind.comm.2013), 
                 group = wind.arth.2013$Wind.Exposure,
                 bias.adjust = TRUE)) # violates assumption of adonis
plot(betadisper(vegdist(wind.comm.2013), 
                 group = wind.arth.2013$Genotype,
                 bias.adjust = TRUE)) # doesn't violates assumption of adonis...


## ant-aphid growth rates ----
# no relationship with distance to ant mound
ggplot(aa.aphid.GR, aes(x = Genotype, y = Aphis.growth.rate, color = Genotype)) +
  geom_boxplot() +
  geom_point(aes(fill = Genotype), position = position_jitterdodge(jitter.width = 2)) + 
  #stat_smooth(se = FALSE) +
  facet_wrap(~Date_rel, nrow = 2)

# residuals look a bit binomial in their distribution. I think this is because the different durations put different lower limits to the decline in population growth rate. Now, I'm thinking the best thing to do would be to just used Aphid densities over time. It is a lot more of an intuitive response variable, and I could account for zeros with a poisson or maybe zero-inflated poisson.
aphid.GR.lmer <- lmer(Aphis.growth.rate ~ Ant.mound.dist + Genotype*Date_rel + (1|Block/Ant.mound.dist) + (1|plant_code),
                      aa.aphid.GR)
summary(aphid.GR.lmer)
#confint(profile(aphid.GR.lmer))
plot(aphid.GR.lmer) # things that need to be accounted for.
anova(aphid.GR.lmer, ddf = "Kenward-Roger")


#mutate(non.leaftier.herb.abund = LTF_Caloptilia + tentmine_Phyllonorycter + gall_R_rigidae + gall_R_salicisbrassicoides + gall_Pontania + gall_Aculus + leafhopper_C_reductus + leafhopper_green + leafhopper_unk + sawfly_larva + caterpillar_looper + caterpillar_LB + caterpillar_unk + red_scale + psyllid + grasshopper,
#      total.herb.abund = non.leaftier.herb.abund + leaftier_Tortricid,
#     total.omniv.abund = stinkbug + ant_F_obscuripes,
#    total.pred.abund = spider_Theridion + spider_BY + spider_NW + spider_Tetragnathid + spider_CS + spider_Larionoides)

#herbs <- colnames(wind.2013.vis.df.max)[5:21]
#omnivs <- colnames(wind.2013.vis.df.max)[22:23]
#preds <- colnames(wind.2013.vis.df.max)[24:29]

#wind.2013.max.gg <- wind.2013.vis.df.max %>%
# filter(Dead < 1) %>%
#gather(Species, Abundance, LTF_Caloptilia:total.pred.abund)

#library(ggplot2)
#ggplot(filter(wind.2013.max.gg, Species %in% c("leaftier_Tortricid", "non.leaftier.herb.abund", "total.pred.abund")),
#      aes(x = Genotype, y = Abundance, color = Wind.Exposure)) +
#geom_boxplot() + 
#facet_wrap(~Species, nrow = 3)

# graph doesn't look too convincing. Aphid population growth is virtually always negative...
ggplot(filter(aa.2012.vis.aphidgrowth, Aphid.Treatment == "aphid"),
       aes(x = Genotype, y = Aphid.growth.rate, color = Genotype)) +
  geom_boxplot() +
  facet_wrap(~relative.Date, ncol = 2)

ggplot(filter(aa.2012.vis.alive, Aphid.Treatment == "aphid"),
       aes(x = total.aphid.abund, y = Red.Ants)) +
  geom_point() +
  stat_smooth(method = "loess")

library(lme4)
aphid.lmer <- glmer(Red.Ants.pres ~ total.aphid.abund + (1|Genotype) + (1|relative.Date) + (1|Block/Distance.to.Ant.Mound), aa.2012.vis.alive, family = "binomial")
summary(aphid.lmer)
confint(profile(aphid.lmer))
anova(aphid.lmer)
plot(aphid.lmer)
plot(total.aphid.abund ~ relative.Date, aa.2012.vis.alive)


