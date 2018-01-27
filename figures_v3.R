
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)
library(coda)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# color-blind friendly palette
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

aa.partition.2012 <- aa.trait.2012
aa.partition.2012$Trait.PC1 <- aa.partition.2012$Trait.PC1*aa.trait.rich.2012$sd_sc.Trait.PC1
aa.partition.2012$Trait.PC2 <- aa.partition.2012$Trait.PC2*aa.trait.rich.2012$sd_sc.Trait.PC2
aa.partition.2012$term_alt <- factor(c("Block","Genotype (G)","G x E","G x E", "G x E", "Plot", "Ant", "Aphid", "Aphid x Ant"),
                                     levels=c("Genotype (G)","Aphid","Ant","Aphid x Ant","G x E","Block","Plot"))

group_by(aa.partition.2012, Variation) %>% summarise(total.indirect.effect = sum(Trait.PC1) + sum(Trait.PC2))
0.296/0.146 # 2x
aa.pal <- c("#009E73","#F0E442","#D55E00","#CC79A7","#E69F00","#000000","#999999")
plot_aa.partition.2012 <- aa.partition.2012 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -term_alt) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt))  +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse=TRUE)) +
  scale_fill_manual(values=aa.pal, name="Source of\nTrait variation") + #  labels=c("Block", "Genotype (G)", "G x Ant", "G x Aphid", "G x Aphid x Ant", "Plot", "Ant", "Aphid", "Aphid x Ant")
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Indirect effect on arthropod richness")

save_plot(filename = "fig_1_ant-aphid_arthropods.png", plot = plot_aa.partition.2012, base_height = 6, base_width = 8.5)

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
wind.pal <- c("#009E73","#56B4E9","#E69F00","#000000","#999999")

w.trait.2012 <- wind.trait.SD.modes %>%
  select(-Root.C.N) %>%
  filter(Year == "2012", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.rich.2012 <- wind.arth.SD.modes %>%
  filter(term %in% c("sd_sc.Trait.PC1","sd_sc.Trait.PC2"), Year == "2012") %>%
  spread(term, Arthropod.Richness)

w.partition.2012 <- w.trait.2012
w.partition.2012$Trait.PC1 <- w.trait.2012$Trait.PC1*w.trait.rich.2012$sd_sc.Trait.PC1
w.partition.2012$Trait.PC2 <- w.trait.2012$Trait.PC2*w.trait.rich.2012$sd_sc.Trait.PC2
w.partition.2012$term_alt <- factor(c("Block","Genotype (G)","G x E","Plot", "Wind"),
                                     levels=c("Genotype (G)","Wind","G x E","Block","Plot"))

group_by(w.partition.2012, Variation) %>% summarise(total.indirect.effect = sum(Trait.PC1) + sum(Trait.PC2))
0.1305/0.0865 # 1.5x
w.partition.2012 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -term_alt) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt))  +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse=TRUE)) +
  scale_fill_manual(values=wind.pal, name="Source of\nTrait variation") + 
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Indirect effect on arthropod richness")

## Plot sources of Trait PC1 and PC2 effects on Arthropod Richness in Wind 2013
w.trait.2013 <- wind.trait.SD.modes %>%
  select(-Root.C.N) %>%
  filter(Year == "2013", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.rich.2013 <- wind.arth.SD.modes %>%
  filter(term %in% c("sd_sc.Trait.PC1","sd_sc.Trait.PC2"), Year == "2013") %>%
  spread(term, Arthropod.Richness)

w.partition.2013 <- w.trait.2013
w.partition.2013$Trait.PC1 <- w.trait.2013$Trait.PC1*w.trait.rich.2013$sd_sc.Trait.PC1
w.partition.2013$Trait.PC2 <- w.trait.2013$Trait.PC2*w.trait.rich.2013$sd_sc.Trait.PC2
w.partition.2013$term_alt <- factor(c("Block","Genotype (G)","G x E","Plot", "Wind"),
                                    levels=c("Genotype (G)","Wind","G x E","Block","Plot"))
group_by(w.partition.2013, Variation) %>% summarise(total.indirect.effect = sum(Trait.PC1) + sum(Trait.PC2))
0.4643/0.1290 # 3.6x
w.partition.2013 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -term_alt) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt))  +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=wind.pal, name="Source of\nTrait variation") + #, labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_wrap(~Traits, ncol=2) +
  xlab("Intraspecific variation") + ylab("Indirect effect on arthropod richness")

plot_w.partitions <- bind_rows(w.partition.2012, w.partition.2013) %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -term_alt) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt))  +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=wind.pal, name="Source of\nTrait variation") + #, labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_grid(Year~Traits) +
  xlab("Intraspecific variation") + ylab("Indirect effect on arthropod richness")

save_plot(filename = "fig_2_wind_arthropods.png", plot = plot_w.partitions, base_height = 8.5, base_width = 8.5)

## Plot sources of Root C:N effects on Fungi and Bacteria Rarefied Richness in Wind 2013
w.root.2013 <- wind.trait.SD.modes %>%
  select(-Trait.PC1, -Trait.PC2) %>%
  filter(Year == "2013", term != "sigma") %>%
  mutate(Variation = c("Plastic", "Heritable", rep("Plastic",3)))

w.trait.microbe.2013 <- filter(wind.below.SD.modes, term == "sd_sc.log.Root.CN")  

w.root.Fungi.2013 <- w.root.2013  
w.root.Fungi.2013$Root.C.N <- w.root.2013$Root.C.N*w.trait.microbe.2013$Fungi.Rarefied.Richness
w.root.Fungi.2013$term_alt <- factor(c("Block","Genotype (G)","G x E","Plot", "Wind"),
                                    levels=c("Genotype (G)","Wind","G x E","Block","Plot"))
group_by(w.root.Fungi.2013, Variation) %>% summarise(total.indirect.effect = sum(Root.C.N))
0.161/0.023 # 7x

w.root.Bacteria.2013 <- w.root.2013
w.root.Bacteria.2013$Root.C.N <- w.root.2013$Root.C.N*w.trait.microbe.2013$Bacteria.Rarefied.Richness
w.root.Bacteria.2013$term_alt <- factor(c("Block","Genotype (G)","G x E","Plot", "Wind"),
                                     levels=c("Genotype (G)","Wind","G x E","Block","Plot"))
group_by(w.root.Bacteria.2013, Variation) %>% summarise(total.indirect.effect = sum(Root.C.N))
0.035/0.005 # 7x

w.root.microbes.2013 <- bind_rows(mutate(w.root.Fungi.2013, Community = "Root-fungi"),
                                  mutate(w.root.Bacteria.2013, Community = "Root-bacteria"))

plot_w.root.microbes.2013 <- w.root.microbes.2013 %>%
  gather(Traits, Effect_Size, -term, -Year, -Variation, -Community, -term_alt) %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt))  +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse=TRUE)) +
  scale_fill_manual(values=wind.pal, name="Source of\nRoot C:N variation") + #, labels=c("Block", "Genotype (G)", "G x Wind", "Plot", "Wind")) +
  facet_wrap(~Community) +
  ylab("Indirect effect on microbial rarefied richness") +
  xlab("Intraspecific variation")
save_plot(filename = "fig_3_microbes.png", plot = plot_w.root.microbes.2013, base_height = 6, base_width = 8.5)


var_data <- bind_rows(gather(mutate(aa.partition.2012, Community = "Arthropods", Experiment_Year = "Ant-Aphid 2012"), Traits, Effect_Size, -term, -Year, -Variation, -Community, -term_alt, -Experiment_Year),
          gather(mutate(w.partition.2012, Community = "Arthropods", Experiment_Year = "Wind 2012"), Traits, Effect_Size, -term, -Year, -Variation, -Community, -term_alt, -Experiment_Year),
          gather(mutate(w.partition.2013, Community = "Arthropods", Experiment_Year = "Wind 2013"), Traits, Effect_Size, -term, -Year, -Variation, -Community, -term_alt, -Experiment_Year),
          gather(mutate(w.root.microbes.2013, Experiment_Year = "Wind 2013"), Traits, Effect_Size, -term, -Year, -Variation, -Community, -term_alt, -Experiment_Year)) %>%
  filter(Community %in% c("Arthropods","Root-fungi")) # only keeping Root-Fungi because it Root C:N had the largest effect on it, plus the relative importance of heritable vs. plastic trait variation is the same as for root-microbes
unique(var_data$term_alt)
var_data$term_alt <- factor(var_data$term_alt)
levels(var_data$term_alt)
var_data$Traits <- factor(var_data$Traits)
levels(var_data$Traits) 
var_data$Traits <- factor(var_data$Traits, levels=c("Trait.PC1","Trait.PC2","Root.C.N"), labels=c("Trait PC1","Trait PC2", "Root C:N"))
var_data$term_alt_ord <- factor(var_data$term_alt, levels = c("Genotype (G)", "Aphid", "Ant", "Aphid x Ant", "Wind", "G x E", "Block", "Plot"))

lanphere_pal <- c("#009E73","#F0E442","#D55E00","#CC79A7","#56B4E9","#E69F00","#000000","#999999")

plot_var <- var_data %>%
  ggplot(., aes(x=Variation, y=Effect_Size, fill=term_alt_ord)) +
  geom_bar(stat="identity", color = "black", position = position_stack(reverse=T)) +
  scale_fill_manual(values = lanphere_pal, name = "Source of\nTrait Variation") +
  facet_grid(Experiment_Year~Traits) +
  ylab("Indirect effect on species richness") +
  xlab("")
save_plot(filename = "fig_heritable_vs_plastic.png", plot = plot_var, base_height = 6, base_width = 8.5)


## PCA
pc.df <- data.frame(Trait.PC1 = c(0.51, 0.43, 0.64, -0.12, -0.36, 0.55, 0.47, 0.68, -0.08, 0.09, 0.46, 0.49, 0.45, -0.36, -0.46, -0.04),
                    Trait.PC2 = c(-0.49, 0.25, -0.14, -0.72, -0.39, -0.18, 0.29, -0.07, 0.64, 0.69, -0.29, -0.27, -0.20, -0.52, -0.42, 0.59),
                    Traits = c(rep(c("Plant height", "Shoot count", "Shoot length", "Leaf trichome density", "Leaf water content"),2),c("Plant height","Shoot count","Shoot length","Leaf water content","Specific leaf area","Leaf C:N")),
                    Experiment = c(rep("Ant-Aphid",5),rep("Wind",5),rep("Wind",6)),
                    Year = c(rep("2012",5),rep("2012",5),rep("2013",6)))
gather.pc.df <- gather(pc.df, Principal_Component, Value, -Traits, -Experiment, -Year) %>%
  unite(Experiment_Year, Experiment, Year, sep = " ")
gather.pc.df$Traits_ord <- factor(gather.pc.df$Traits, levels=rev(c("Plant height","Shoot count","Shoot length","Leaf water content","Leaf trichome density","Specific leaf area","Leaf C:N")))
plot_pca <- ggplot(gather.pc.df, aes(x=Traits_ord, y=Value)) +
  geom_hline(yintercept = 0, linetype="solid") +
  geom_bar(stat="identity") + 
  coord_flip() + 
  facet_grid(Experiment_Year ~ Principal_Component, drop=T, scales="free_y") +
  xlab("") + ylab("Loadings")
save_plot(filename = "fig_pca.png", plot = plot_pca, base_height=6, base_width=8.5)
