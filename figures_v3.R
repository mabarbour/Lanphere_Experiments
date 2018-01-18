
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)
library(coda)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## GET ANT-APHID DATA ----
ant.aphid.trait.SD <- read.csv('output_brms/ant.aphid_SDs.csv') %>%
  gather(key = term, value = posterior_SD, -sample, -(Experiment:Response)) %>%
  filter(Response %in% c("Trait PC1", "Trait PC2")) %>%
  group_by(term, Response, Year) %>%
  summarise(mode = Mode(round(posterior_SD,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,2]) %>%
  ungroup()

ant.aphid.trait.SD.modes <- ant.aphid.trait.SD %>%
  select(term, Year, Response, mode) %>%
  spread(key=Response, value=mode) %>%
  filter(term != "sd_Intercept") %>% # non-informative
  data.frame()
ant.aphid.trait.SD.modes


ant.aphid.arth.SD <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  select(Experiment, Year, Response, sample, sd_sc.Trait.PC1, sd_sc.Trait.PC2, 
         sd_sc.Aphid.treatment, sd_sc.Ant.mound.dist, sd_sc.Aphid.treatment.sc.Ant.mound.dist, 
         sd_Genotype__Intercept, sd_Genotype__sc.Aphid.treatment, sd_Genotype__sc.Ant.mound.dist, sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist,
         sd_Block__Intercept, sd_Plot_code__Intercept) %>%
  gather(key = term, value = posterior_SD, -sample, -(Experiment:Response)) %>%
  filter(Response == "Arthropod Richness", Experiment == "Ant-Aphid") %>%
  group_by(term, Response, Year) %>%
  summarise(mode = Mode(round(posterior_SD,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,2]) 
ant.aphid.arth.SD.modes <- ant.aphid.arth.SD %>%
  select(term, Year, Response, mode) %>%
  spread(key=Response, value=mode) %>%
  filter(term != "sd_Intercept") %>% # non-informative
  data.frame()
ant.aphid.arth.SD.modes

aa.trait.2012 <- ant.aphid.trait.SD.modes %>%
  filter(Year == "2012", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",7)))

aa.trait.rich.2012 <- ant.aphid.arth.SD.modes %>%
  filter(term %in% c("sd_sc.Trait.PC1","sd_sc.Trait.PC2"), Year == "2012") %>%
  spread(term, Arthropod.Richness)

aa.trait.2012$Trait.PC1 <- aa.trait.2012$Trait.PC1*aa.trait.rich.2012$sd_sc.Trait.PC1
aa.trait.2012$Trait.PC2 <- aa.trait.2012$Trait.PC2*aa.trait.rich.2012$sd_sc.Trait.PC2
aa.trait.2012 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term))  +
  geom_bar(stat="identity", color = "black", position = position_stack()) +
  scale_fill_manual(values=cbPalette, name="Source of\nTrait variation", labels=c("Block", "Genotype (G)", "G x Ant", "G x Aphid", "G x Aphid x Ant", "Plot", "Ant", "Aphid", "Aphid x Ant")) +
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Trait effect size on arthropod richness")

## GET WIND DATA ----
wind.trait.SD <- read.csv('output_brms/wind_SDs.csv') %>%
  gather(key = term, value = posterior_SD, -sample, -(Experiment:Response)) %>%
  filter(Response %in% c("Trait PC1", "Trait PC2", "Root C:N")) %>%
  group_by(term, Response, Year) %>%
  summarise(mode = Mode(round(posterior_SD,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,2]) %>%
  ungroup()

wind.trait.SD.modes <- wind.trait.SD %>%
  select(term, Year, Response, mode) %>%
  spread(key=Response, value=mode) %>%
  filter(term != "sd_Intercept") %>% # non-informative
  data.frame()
wind.trait.SD.modes


wind.arth.SD <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  select(Experiment, Year, Response, sample, sd_sc.Trait.PC1, sd_sc.Trait.PC2, 
         sd_sc.Wind.Exposure, sd_Genotype__Intercept, sd_Genotype__sc.Wind.Exposure, 
         sd_Block__Intercept, sd_Plot_code__Intercept) %>%
  gather(key = term, value = posterior_SD, -sample, -(Experiment:Response)) %>%
  filter(Response == "Arthropod Richness", Experiment == "Wind") %>%
  group_by(term, Response, Year) %>%
  summarise(mode = Mode(round(posterior_SD,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,2]) 
wind.arth.SD.modes <- wind.arth.SD %>%
  select(term, Year, Response, mode) %>%
  spread(key=Response, value=mode) %>%
  filter(term != "sd_Intercept") %>% # non-informative
  data.frame()
wind.arth.SD.modes

wind.below.SD <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  select(Experiment, Year, Response, sample, sd_sc.log.Root.CN, sd_sc.Wind.Exposure, 
         sd_Genotype__Intercept, sd_Genotype__sc.Wind.Exposure, 
         sd_Block__Intercept, sd_Plot_code__Intercept) %>%
  gather(key = term, value = posterior_SD, -sample, -(Experiment:Response)) %>%
  filter(Response %in% c("Fungi Rarefied Richness","Bacteria Rarefied Richness")) %>%
  group_by(term, Response) %>%
  summarise(mode = Mode(round(posterior_SD,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(posterior_SD), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(posterior_SD), prob=0.95)[ ,2]) 
wind.below.SD.modes <- wind.below.SD %>%
  select(term, Response, mode) %>%
  spread(key=Response, value=mode) %>%
  data.frame()
wind.below.SD.modes


## Plot sources of Trait PC1 and PC2 effects on Arthropod Richness in Wind 2012
w.trait.2012 <- wind.trait.SD.modes %>%
  select(-Root.C.N) %>%
  filter(Year == "2012", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.rich.2012 <- wind.arth.SD.modes %>%
  filter(term %in% c("sd_sc.Trait.PC1","sd_sc.Trait.PC2"), Year == "2012") %>%
  spread(term, Arthropod.Richness)

w.trait.2012$Trait.PC1 <- w.trait.2012$Trait.PC1*w.trait.rich.2012$sd_sc.Trait.PC1
w.trait.2012$Trait.PC2 <- w.trait.2012$Trait.PC2*w.trait.rich.2012$sd_sc.Trait.PC2
w.trait.2012 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term))  +
  geom_bar(stat="identity", color = "black", position = position_stack()) +
  scale_fill_manual(values=cbPalette[c(1,5,6,2,4)], name="Source of\nTrait variation", labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Trait effect size on arthropod richness")

## Plot sources of Trait PC1 and PC2 effects on Arthropod Richness in Wind 2013
w.trait.2013 <- wind.trait.SD.modes %>%
  select(-Root.C.N) %>%
  filter(Year == "2013", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.rich.2013 <- wind.arth.SD.modes %>%
  filter(term %in% c("sd_sc.Trait.PC1","sd_sc.Trait.PC2"), Year == "2013") %>%
  spread(term, Arthropod.Richness)

w.trait.2013$Trait.PC1 <- w.trait.2013$Trait.PC1*w.trait.rich.2013$sd_sc.Trait.PC1
w.trait.2013$Trait.PC2 <- w.trait.2013$Trait.PC2*w.trait.rich.2013$sd_sc.Trait.PC2
w.trait.2013 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term))  +
  geom_bar(stat="identity", color = "black", position = position_stack()) +
  scale_fill_manual(values=cbPalette[c(1,5,6,2,4)], name="Source of\nTrait variation", labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Trait effect size on arthropod richness")

## Plot sources of Root C:N effects on Fungi and Bacteria Rarefied Richness in Wind 2013
w.root.2013 <- wind.trait.SD.modes %>%
  select(-Trait.PC1, -Trait.PC2) %>%
  filter(Year == "2013", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.microbe.2013 <- filter(wind.below.SD.modes, term == "sd_sc.log.Root.CN")  

w.root.Fungi.2013 <- w.root.2013  
w.root.Fungi.2013$Root.C.N <- w.root.2013$Root.C.N*w.trait.microbe.2013$Fungi.Rarefied.Richness

w.root.Bacteria.2013 <- w.root.2013
w.root.Bacteria.2013$Root.C.N <- w.root.2013$Root.C.N*w.trait.microbe.2013$Bacteria.Rarefied.Richness

w.root.microbes.2013 <- bind_rows(mutate(w.root.Fungi.2013, Community = "Fungi rarefied richness"),
                                  mutate(w.root.Bacteria.2013, Community = "Bacteria rarefied richness"))

w.root.microbes.2013 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -Community) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term))  +
  geom_bar(stat="identity", color = "black", position = position_stack()) +
  scale_fill_manual(values=cbPalette[c(1,5,6,2,4)], name="Source of\nRoot C:N variation", labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_wrap(~Community) +
  ylab("Root C:N effect size on Root-Microbes") +
  xlab("Intraspecific variation")


## LIKELY MORGUE ----
heritable.wind.arth.2012 <- wind.trait.SD.modes %>%
  filter(Year=="2012", term=="sd_Genotype__Intercept") %>%
  mutate(Variation.Source = "Heritable traits")
heritable.wind.arth.2012$Trait.PC1 <- heritable.wind.arth.2012$Trait.PC1 * filter(wind.arth.SD.modes, Year=="2012", term=="sd_sc.Trait.PC1")$Arthropod.Richness
heritable.wind.arth.2012$Trait.PC2 <- heritable.wind.arth.2012$Trait.PC2 * filter(wind.arth.SD.modes, Year=="2012", term=="sd_sc.Trait.PC2")$Arthropod.Richness
heritable.wind.arth.2012 <- heritable.wind.arth.2012 %>%
  gather(traits, effect_size, -term, -Year, -Variation.Source)
heritable.wind.arth.2012.unk <- wind.arth.SD.modes %>%
  filter(Year=="2012", term=="sd_Genotype__Intercept") %>%
  mutate(traits = NA, Variation.Source = "Heritable traits") %>%
  rename(effect_size = Arthropod.Richness)
heritable.wind.arth.2012 <- bind_rows(heritable.wind.arth.2012, heritable.wind.arth.2012.unk)

plasticity.wind.arth.2012 <- wind.trait.SD.modes %>%
  filter(Year=="2012", term%in%c("sd_sc.Wind.Exposure","sd_Genotype__sc.Wind.Exposure","sd_Block__Intercept","sd_Plot_code__Intercept")) %>%
  mutate(Variation.Source = "Plastic traits")
plasticity.wind.arth.2012$Trait.PC1 <- plasticity.wind.arth.2012$Trait.PC1 * filter(wind.arth.SD.modes, Year=="2012", term=="sd_sc.Trait.PC1")$Arthropod.Richness
plasticity.wind.arth.2012$Trait.PC2 <- plasticity.wind.arth.2012$Trait.PC2 * filter(wind.arth.SD.modes, Year=="2012", term=="sd_sc.Trait.PC2")$Arthropod.Richness
plasticity.wind.arth.2012 <- plasticity.wind.arth.2012 %>%
  gather(traits, effect_size, -term, -Year, -Variation.Source)

environment.wind.arth.2012 <- wind.arth.SD.modes %>%
  filter(Year=="2012", term%in%c("sd_sc.Wind.Exposure","sd_Genotype__sc.Wind.Exposure","sd_Block__Intercept","sd_Plot_code__Intercept")) %>%
  mutate(Variation.Source = "Direct/Indirect Environment", traits = NA) %>%
  rename(effect_size = Arthropod.Richness)
environment.wind.arth.2012

bind_rows(heritable.wind.arth.2012, plasticity.wind.arth.2012, environment.wind.arth.2012) %>%
  ggplot(., aes(x = Variation.Source, y=effect_size, fill=traits)) +
  geom_bar(stat="identity", position = position_stack())

bind_rows(heritable.wind.arth.2012, plasticity.wind.arth.2012, environment.wind.arth.2012) %>%
  ggplot(., aes(x = traits, y=effect_size, fill=Variation.Source)) +
  geom_bar(stat="identity", position = position_stack())

bind_rows(heritable.wind.arth.2012, plasticity.wind.arth.2012, environment.wind.arth.2012) %>%
  ggplot(., aes(x = traits, y=effect_size, fill=term)) +
  geom_bar(stat="identity", position = position_stack()) +
  coord_flip() + 
  facet_wrap(~Variation.Source, ncol=1)

cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bind_rows(heritable.wind.arth.2012, plasticity.wind.arth.2012, environment.wind.arth.2012) %>%
  filter(traits %in% c("Trait.PC1","Trait.PC2")) %>%
  ggplot(., aes(x = Variation.Source, y=effect_size, fill=term)) +
  geom_bar(stat="identity", color = "black", position = position_stack()) +
  scale_fill_manual(values=cbPalette[c(1,5,6,2,4)]) +
  #coord_flip() +
  facet_wrap(~traits, ncol=2)
