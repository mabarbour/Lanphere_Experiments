
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)
library(coda)

## GET DATA ----
wind.var.df <- read.csv('output_brms/wind_SDs.csv') %>% 
  transmute("Genotype (G)" = sd_Genotype__Intercept^2 /        (sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Genotype__sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Wind =       sd_sc.Wind.Exposure^2 /           (sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Genotype__sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "G x Wind" =     sd_Genotype__sc.Wind.Exposure^2 / (sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Genotype__sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =    sd_Block__Intercept^2 /           (sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Genotype__sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =     sd_Plot_code__Intercept^2 /       (sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Genotype__sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "experimental_design")

wind.var.soil <- read.csv('output_brms/wind_SDs.csv') %>% 
  filter(Response%in%c("scale(log(Soil PC1 + min(Soil PC1) + 1))", "scale(Soil PC2)")) %>%
  transmute(Wind =       sd_sc.Wind.Exposure^2 /           (sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sigma),
            Block =      sd_Block__Intercept^2 /           (sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "experimental_design")

aa.var.df <- read.csv('output_brms/ant.aphid_SDs.csv') %>% 
  transmute("Genotype (G)" =    sd_Genotype__Intercept^2 /                            (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Aphid =             sd_sc.Aphid.treatment^2 /                             (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Ant =               sd_sc.Ant.mound.dist^2 /                              (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Aphid x Ant" =     sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 /           (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "G x Aphid" =       sd_Genotype__sc.Aphid.treatment^2 /                   (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "G x Ant" =         sd_Genotype__sc.Ant.mound.dist^2 /                    (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "G x Aphid x Ant" = sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 / (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =             sd_Block__Intercept^2 /                               (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =              sd_Plot_code__Intercept^2 /                           (sd_Genotype__Intercept^2 + sd_sc.Aphid.treatment^2 + sd_sc.Ant.mound.dist^2 + sd_sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment^2 + sd_Genotype__sc.Ant.mound.dist^2 + sd_Genotype__sc.Aphid.treatment.sc.Ant.mound.dist^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "experimental_design")

aa.arth.trait.var <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  filter(Experiment == "Ant-Aphid") %>%
  transmute("Trait PC1" =       sd_sc.Trait.PC1^2 /                                   (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Trait PC2" =       sd_sc.Trait.PC2^2 /                                   (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Genotype (G)" =    sd_Genotype__Intercept^2 /                            (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =             sd_Block__Intercept^2 /                               (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =              sd_Plot_code__Intercept^2 /                           (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "trait_reg")

w.arth.trait.var.2012 <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  filter(Experiment == "Wind", Year == "2012") %>%
  transmute("Trait PC1" =       sd_sc.Trait.PC1^2 /                                   (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Trait PC2" =       sd_sc.log.trans.Trait.PC2^2 /                         (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Genotype (G)" =    sd_Genotype__Intercept^2 /                            (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Wind =              sd_sc.Wind.Exposure^2/                                (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =             sd_Block__Intercept^2 /                               (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =              sd_Plot_code__Intercept^2 /                           (sd_sc.Trait.PC1^2 + sd_sc.log.trans.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "trait_reg") 

w.arth.trait.var.2013 <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  filter(Experiment == "Wind", Year == "2013", Response == "scale(log(Arthropod Richness + 1))") %>%
  transmute("Trait PC1" =       sd_sc.Trait.PC1^2 /                                   (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Trait PC2" =       sd_sc.Trait.PC2^2 /                                   (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Genotype (G)" =    sd_Genotype__Intercept^2 /                            (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Wind =              sd_sc.Wind.Exposure^2/                                (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =             sd_Block__Intercept^2 /                               (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =              sd_Plot_code__Intercept^2 /                           (sd_sc.Trait.PC1^2 + sd_sc.Trait.PC2^2 + sd_Genotype__Intercept^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "trait_reg") 

w.microbe.trait.soil.var <- read.csv("output_brms/lanphere_trait_regs.csv") %>%
  filter(Experiment == "Wind", Year == "2013", Response != "scale(log(Arthropod Richness + 1))") %>%
  transmute("Soil PC1" =        sd_sc.log.trans.Soil.PC1^2 /                          (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Soil PC2" =        sd_sc.Soil.PC2^2 /                                    (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            "Root C:N" =        sd_sc.log.Root.CN^2/                                  (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Wind =              sd_sc.Wind.Exposure^2/                                (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Block =             sd_Block__Intercept^2 /                               (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            Plot =              sd_Plot_code__Intercept^2 /                           (sd_sc.log.trans.Soil.PC1^2 + sd_sc.Soil.PC2^2 + sd_sc.log.Root.CN^2 + sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sd_Plot_code__Intercept^2 + sigma),
            sample = sample, Experiment = Experiment, Year = Year, Response = Response) %>%
  gather(key = term, value = percent.variance, -(sample:Response)) %>%
  mutate(Model = "trait_reg")

var.df <- bind_rows(wind.var.df, aa.var.df) %>% #, aa.arth.trait.var, w.arth.trait.var.2012, w.arth.trait.var.2013, w.microbe.trait.soil.var) %>%
  unite(Experiment_Year, Experiment, Year, sep = " ", remove = F)
var.df$term <- factor(var.df$term)
levels(var.df$term)
var.df$term_ord <- factor(var.df$term, levels=c("Plot","Block","G x Aphid x Ant","G x Ant","G x Aphid","G x Wind","Aphid x Ant","Ant","Aphid","Wind","Genotype (G)")) # "Soil PC2","Soil PC1", "Root C:N","Trait PC2","Trait PC1"
levels(var.df$term_ord)

## FIGURE: VARIANCE EXPLAINED IN PLANT TRAITS ----
plot_traits <- filter(var.df, Model%in%"experimental_design", Response%in%c("Trait PC1","Trait PC2","Root C:N"))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

summary_traits <- plot_traits %>%
  group_by(Experiment_Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])
summary_traits$Response_ord <- factor(summary_traits$Response, levels = c("Trait PC1","Trait PC2", "Root C:N"))

theme_set(theme_grey())
traits_gg <- ggplot(summary_traits, aes(x=term_ord, y=mode)) + # , fill=Response_ord
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=3, shape=21, fill="grey", color="black", position=position_dodge(width=0.75)) + #
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) +
  #scale_shape_manual(values = c(21,22,23), name="Plant trait") +
  facet_grid(Experiment_Year ~ Response_ord, scales="free_y") + 
  theme_cowplot()# + theme(panel.border = element_rect(colour="black")) # not working for some reason
traits_gg 
save_plot(filename = "fig_trait_variance.png", plot = traits_gg, base_height = 6, base_width = 8)

# something weird is going on when I try to use wind.var.soil in above var.df. It keeps Genotype and G x Wind and reports NA
var.soil <- unite(wind.var.soil, Experiment_Year, Experiment, Year, sep = " ", remove = F)
var.soil$term <- factor(var.soil$term)
levels(var.soil$term)
var.soil$term_ord <- factor(var.soil$term, levels=c("Block","Wind"))
levels(var.soil$term_ord)

## FIGURE 1: COMMUNITY RICHNESS ----
plot_richness <- filter(var.df, Model%in%"experimental_design", Response%in%c("scale(log(Arthropod Richness + 1))", "scale(Fungi Rarefied Richness)", "scale(Bacteria Rarefied Richness)"))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

summary_richness <- plot_richness %>%
  group_by(Experiment_Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])
summary_richness$Response_ord <- factor(summary_richness$Response, levels = c("scale(log(Arthropod Richness + 1))","scale(Fungi Rarefied Richness)", "scale(Bacteria Rarefied Richness)"))

richness_gg <- ggplot(summary_richness, aes(x=term_ord, y=mode, shape=Response_ord)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=3, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), labels=c(0, 10, 20, 30)) +
  scale_shape_manual(values = c(21,22,23), name="Community response", labels=c("Arthropod richness","Fungi rarefied richness","Bacteria rarefied richness")) +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y") 
richness_gg

## SUPP MAT 1 ----

# Color-blind friendly palette:
#cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_traits <- filter(var.df, Model%in%"experimental_design", Response%in%c("scale(Trait PC1)", "scale(log(Trait PC2 + min(Trait PC2) + 1))", "scale(Trait PC2)", "scale(log(Root C:N))"))
plot_traits$Response[which(plot_traits$Response == "scale(log(Trait PC2 + min(Trait PC2) + 1))")] <- "scale(Trait PC2)" # for labeling purposes
unique(plot_traits$Response)

summary_traits <- plot_traits %>%
  group_by(Experiment_Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])
summary_traits$Response_ord <- factor(summary_traits$Response, levels = c("scale(Trait PC1)", "scale(Trait PC2)", "scale(log(Root C:N))"))

traits_gg <- ggplot(summary_traits, aes(x=term_ord, y=mode, shape=Response_ord)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) +
  scale_shape_manual(values = c(21,22,23,24), name="Trait response", labels=c("Trait PC1","Trait PC2","Root C:N")) +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y")  
traits_gg

ggplot(summary_traits, aes(x=term_ord, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) +
  #scale_shape_manual(values = c(21,22,23,24), name="Trait response", labels=c("Trait PC1","Trait PC2","Root C:N")) +
  facet_grid(Experiment_Year~Response_ord, scales="free_y")

summary_soil <- var.soil %>%
  group_by(Experiment_Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])
summary_soil$Response_ord <- factor(summary_soil$Response, levels = c("scale(log(Soil PC1 + min(Soil PC1) + 1))", "scale(Soil PC2)"))

soil_gg <- ggplot(summary_soil, aes(x=term_ord, y=mode, shape=Response_ord)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), labels=c(0, 20, 40, 60)) +
  scale_shape_manual(values = c(21,22), name="Soil response", labels=c("Soil PC1","Soil PC2")) +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y")
soil_gg

plot_grid(traits_gg, soil_gg, ncol=1, align = 'v', labels = "AUTO", rel_heights=c(1,0.25))

## TRAIT - COMMUNITY RELATIONSHIPS ----

summary_aa_trait.reg <- filter(var.df, Model%in%"trait_reg", Experiment == "Ant-Aphid") %>% #Experiment == "Ant-Aphid", Year=="2012") %>%
  group_by(Experiment_Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])
  
aa.trait.reg_gg <- ggplot(summary_aa_trait.reg, aes(x=term_ord, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) #+
  #scale_shape_manual(values = c(21,22)) + #, name="Trait response", labels=c("Trait PC1","Trait PC2","Root C:N")) #+
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y")  
aa.trait.reg_gg

summary_w_trait.reg <- filter(var.df, Model=="trait_reg", Experiment == "Wind", Response == "scale(log(Arthropod Richness + 1))") %>% #Experiment == "Ant-Aphid", Year=="2012") %>%
  group_by(Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])

w.trait.reg_gg <- ggplot(summary_w_trait.reg, aes(x=term_ord, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) #+
  #scale_shape_manual(values = c(21,22)) + #, name="Trait response", labels=c("Trait PC1","Trait PC2","Root C:N")) #+
  facet_wrap(~Year, ncol=1, scales="free_y")  
w.trait.reg_gg

summary_w.below_trait.reg <- filter(var.df, Model%in%"trait_reg", Response != "scale(log(Arthropod Richness + 1))") %>% #Experiment == "Ant-Aphid", Year=="2012") %>%
  group_by(Experiment, Year, Response, term, term_ord) %>%
  summarise(mode = Mode(round(percent.variance,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(percent.variance), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(percent.variance), prob=0.95)[ ,2])

w.below.trait.reg_gg <- ggplot(summary_w.below_trait.reg, aes(x=term_ord, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Variance Explained (%)") +
  xlab("") +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 10, 20, 30, 40, 50)) #+
  #scale_shape_manual(values = c(21,22)) + #, name="Trait response", labels=c("Trait PC1","Trait PC2","Root C:N")) #+
  facet_wrap(~Response, ncol=1, scales="free_y")  
w.below.trait.reg_gg

plot_grid(aa.trait.reg_gg, w.trait.reg_gg, w.below.trait.reg_gg, ncol=3, align="hv")

## NEED RELEVANT TRAIT/PC DATA FOR VARIANCE EXPLAINED AND LOADINGS
## NEED REGRESSIONS FROM BRM ANALYSIS FOR REGRESSION PLOTS

plot_VarExp <- function(Resp, Exper, Yr){
  response_VarExp <- filter(varcomp.df, Response%in%Resp, Experiment%in%Exper, Year%in%Yr)
  response_VarExp$term_group_ord <- factor(response_VarExp$term_group_ord, levels=rev(levels(response_VarExp$term_group_ord)))
  
  ggplot(response_VarExp, aes(x=Response, y=VarComp_mean*100, fill=term_group_ord)) +
    geom_bar(stat = "identity", position=position_stack(reverse=T), color=NA, show.legend = F) + 
    coord_flip() +
    geom_text(aes(label=term_group_ord), stat = "identity", 
              position=position_stack(reverse=T, vjust=0.5),
              angle=45, show.legend=F, size=3, check_overlap=T) + 
    scale_fill_grey(start=0.5, end=0.9)# +
  #xlab("") + 
  #scale_x_continuous(labels=NULL, breaks=NULL)
}

## ANT-APHID 2012 -- ARTHROPOD RICHNESS & TRAIT PC1 ----
VarExp_aa.trait.PC1.12 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Ant-Aphid", Yr = "2012") + ylab("Trait PC1 explained variance (%)")+
  theme_cowplot(font_size=10)

Arrow_aa.trait.PC1.12 <- aa.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC1, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(-0.1,-0.1,-0.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC1_aa.12 <- plot(marginal_effects(trait.rich.aa.2012.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                                points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_aa.trait.PC1.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_aa.trait.PC1.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

aa.traitreg.pc1 <- plot_grid(ArthRich.TraitPC1_aa.12, VarExp_aa.trait.PC1.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 

## ANT-APHID 2012 -- ARTHROPOD RICHNESS & TRAIT PC2 ----
VarExp_aa.trait.PC2.12 <- plot_VarExp(Resp = "Trait.PC2", Exper = "Ant-Aphid", Yr = "2012") + ylab("Trait PC2 explained variance (%)")+
  theme_cowplot(font_size=10)

Arrow_aa.trait.PC2.12 <- aa.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC2, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(1.1,-0.1,1.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC2_aa.12 <- plot(marginal_effects(trait.rich.aa.2012.brm, effects="trait.PC2", probs=c(0.025,0.975)), 
                                points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC2") +
  geom_segment(data=Arrow_aa.trait.PC2.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_aa.trait.PC2.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

aa.traitreg.pc2 <- plot_grid(ArthRich.TraitPC2_aa.12, VarExp_aa.trait.PC2.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 


## WIND 2012 -- ARTRHOPOD RICHNESS & TRAIT PC1 ----
VarExp_w.trait.PC1.12 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Wind", Yr = "2012") + ylab("Trait PC1 explained variance (%)") +
  theme_cowplot(font_size=10)

Arrow_w.trait.PC1.12 <- wind.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC1, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(-0.1,-0.1,-0.1,1.1,-0.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC1_w.12 <- plot(marginal_effects(trait.rich.wind.2012.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                               points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_w.trait.PC1.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.trait.PC1.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

w.traitreg.pc1 <- plot_grid(ArthRich.TraitPC1_w.12, VarExp_w.trait.PC1.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 


## WIND 2013 -- ARTRHOPOD RICHNESS & TRAIT PC1 ----
VarExp_w.trait.PC1.13 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Wind", Yr = "2013")+ ylab("Trait PC1 explained variance (%)")

Arrow_w.trait.PC1.13 <- wind.trait.PCA.2013 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=6), 
         xend=PC1, yend=seq(-1,-0.25,length.out=6), 
         hjust=c(-0.1,-0.1,-0.1,1.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "SLA", "Water content", "C:N"))

ArthRich.TraitPC1_w.13 <- plot(marginal_effects(trait.rich.wind.2013.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                               points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_w.trait.PC1.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.trait.PC1.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

w.traitreg.pc1.2013 <- plot_grid(ArthRich.TraitPC1_w.13, VarExp_w.trait.PC1.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- BACTERIA RAREFIED RICHNESS & SOIL PC2 ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.soil.PC2.13 <- plot_VarExp(Resp = "Soil.PC2", Exper = "Wind", Yr = "2013") + ylab("Soil PC2 explained variance (%)")

Arrow_w.soil.PC2.13 <- wind.soil.PCA %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=17), 
         xend=PC2, yend=seq(-1,-0.25,length.out=17), 
         hjust=c(1.1,1.1,-0.1,-0.1,1.1,-0.1,1.1,-0.1,1.1,1.1,1.1,-0.1,1.1,1.1,-0.1,-0.1,-0.1),
         label=c("NO_3","NH_4", "Ca","Mg","K","P","Fe","Mn","Cu","Zn","B","S","Pb","Al","Cd","% TOM","Moisture"))

BactRareRich.SoilPC2_w.13 <- plot(marginal_effects(bact.trait.rarerich.wind.2013.brm, effects="soil.PC2", probs=c(0.025,0.975)), 
                                  points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Bacteria rarefied richness") + xlab("Soil PC2") +
  geom_segment(data=Arrow_w.soil.PC2.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.soil.PC2.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

bact.reg.pc2 <- plot_grid(BactRareRich.SoilPC2_w.13, VarExp_w.soil.PC2.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- FUNGI RAREFIED RICHNESS & SOIL PC1 ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.soil.PC1.13 <- plot_VarExp(Resp = "Soil.PC1", Exper = "Wind", Yr = "2013")+ ylab("Soil PC1 explained variance (%)")

Arrow_w.soil.PC1.13 <- wind.soil.PCA %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=17), 
         xend=PC1, yend=seq(-1,-0.25,length.out=17), 
         hjust=c(1.1,1.1,rep(-0.1,15)),
         label=c("NO_3","NH_4", "Ca","Mg","K","P","Fe","Mn","Cu","Zn","B","S","Pb","Al","Cd","% TOM","Moisture"))

FungRareRich.SoilPC2_w.13 <- plot(marginal_effects(trait.rarerich.wind.2013.brm, effects="soil.PC1", probs=c(0.025,0.975)), 
                                  points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Fungi rarefied richness") + xlab("Soil PC1") +
  geom_segment(data=Arrow_w.soil.PC1.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.soil.PC1.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

fung.reg.pc1 <- plot_grid(FungRareRich.SoilPC2_w.13, VarExp_w.soil.PC1.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- FUNGI RAREFIED RICHNESS & ROOT C:N ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.root.CN.13 <- plot_VarExp(Resp = "Root C:N", Exper = "Wind", Yr = "2013") + ylab("Root C:N explained variance (%)")

FungRareRich.RootCN_w.13 <- plot(marginal_effects(trait.rarerich.wind.2013.brm, effects="log_root_CN", probs=c(0.025,0.975)), 
                                 points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Fungi rarefied richness") + xlab("log(Root C:N)") 

fung.reg.rootCN <- plot_grid(FungRareRich.RootCN_w.13, VarExp_w.root.CN.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))



## MULTIPANEL PLOTS ----

aa.traits.var.2012 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2"), Exper = c("Ant-Aphid"), Yr = "2012")+xlab("") + ylab("Variance explained (%)")
w.traits.var.2012 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2"), Exper = c("Wind"), Yr = "2012")+xlab("") + ylab("Variance explained (%)")

## NEED TO ADD WIND TRAITS INTO PLOT, PROBABLY CAN'T DO WITH CURRENT PLOT_VAREXP FUNCTION. MAYBE JUST SELECT THE DATA FRAME AND THEN FEED INTO THE PLOTTING FUNCTION...OR GET RID OF THE FUNCTION BECAUSE I WILL ONLY NEED IT 2X AND NOT 5X LIKE I ORIGINALLY THOUGHT...
test <- plot_grid(ArthRich.TraitPC1_aa.12, ArthRich.TraitPC2_aa.12, 
                  ArthRich.TraitPC1_w.12, aa.traits.var.2012,
                  #VarExp_aa.trait.PC1.12, VarExp_aa.trait.PC2.12, VarExp_w.trait.PC1.12, 
                  align='hv')
save_plot("test.pdf", test, base_height = 8.5, base_width = 8.5)

traits.var.2013 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2", "Root C:N"), Exper = "Wind", Yr = "2013")+xlab("") + ylab("Variance explained (%)") #+ xlab("Plant traits") 
soil.var <- plot_VarExp(Resp=c("Soil.PC1", "Soil.PC2"), Exper = "Wind", Yr = "2013")+xlab("")  + ylab("Variance explained (%)") #+ xlab("Soil properties")


test.wind.2013 <- plot_grid(ArthRich.TraitPC1_w.13, BactRareRich.SoilPC2_w.13, 
                            FungRareRich.RootCN_w.13, FungRareRich.SoilPC2_w.13, 
                            traits.var.2013, soil.var,
                            ncol=2, align = 'hv', labels = "AUTO")
save_plot("test.wind.2013.pdf", test.wind.2013, base_height = 11, base_width = 8.5)
#plot_grid(w.traitreg.pc1.2013, fung.reg.rootCN, bact.reg.pc2, fung.reg.pc1, ncol=2, align = 'hv', labels = "AUTO")

