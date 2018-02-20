wind.var.soil <- read.csv('output_brms/wind_SDs.csv') %>% 
  filter(Response%in%c("scale(log(Soil PC1 + min(Soil PC1) + 1))", "scale(Soil PC2)")) %>%
  transmute(Wind =       sd_sc.Wind.Exposure^2 /           (sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sigma),
            Block =      sd_Block__Intercept^2 /           (sd_sc.Wind.Exposure^2 + sd_Block__Intercept^2 + sigma),
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


plot(marginal_effects(bact.trait.rarerich.wind.2013.brm, effects = "sc.Soil.PC2"), points=T)
plot(marginal_effects(trait.rarerich.wind.2013.brm, effects = "sc.log.trans.Soil.PC1"), points=T)


## WIND SOIL PC1 ANALYSIS ----
hist(w.soil$sc.log.trans.Soil.PC1)

soil.PC1.wind.brm <- general_brm(sc.log.trans.Soil.PC1~sc.Wind.Exposure+(1|Block), data=w.soil)
summary(soil.PC1.wind.brm) 

y_w.soil.PC1 <- w.soil$sc.log.trans.Soil.PC1
yrep_w.soil.PC1 <- posterior_predict(soil.PC1.wind.brm, nsamples=100)
#launch_shinystan(soil.PC1.wind.brm)


## WIND SOIL PC2 ANALYSIS ----
hist(w.soil$sc.Soil.PC2)

soil.PC2.wind.brm <- general_brm(sc.Soil.PC2~sc.Wind.Exposure+(1|Block), data=w.soil)
summary(soil.PC2.wind.brm)

y_w.soil.PC2 <- w.soil$sc.Soil.PC2
yrep_w.soil.PC2 <- posterior_predict(soil.PC2.wind.brm, nsamples=100)
#launch_shinystan(soil.PC2.wind.brm)

## LOOK AT NUTRIENT RELATIONSHIPS
belowground <- left_join(select(fungal.df, Block:fungal.rarerich), select(bacteria.df, plant_ID, bacteria.rarerich), by="plant_ID") %>%
  left_join(., transmute(w.soil, Plot_code, sc.log.trans.Soil.PC1=as.numeric(sc.log.trans.Soil.PC1), sc.Soil.PC2=as.numeric(sc.Soil.PC2), Total.N, NO3.N, NH4.N)) %>%
  left_join(., transmute(w.trait.2013, plant_ID, sc.log.Root.CN = as.numeric(sc.log.Root.CN), root_C.perc, root_N.perc))

ggplot(belowground, aes(x=log(root_N.perc), y=sc.log.Root.CN)) + geom_point()
ggplot(belowground, aes(x=root_C.perc, y=sc.log.Root.CN)) + geom_point()

plot_bg <- group_by(belowground, Plot_code) %>%
  summarise_at(vars(fungal.rarerich, bacteria.rarerich, sc.log.Root.CN, root_N.perc, Total.N, NO3.N, NH4.N, sc.log.trans.Soil.PC1, sc.Soil.PC2), mean, na.rm=T) %>%
  ungroup()

ggplot(w.trait.2013, aes(x=sc.log.Root.CN, y=log(leaf_C_N))) + geom_point() + geom_smooth()

ggplot(plot_bg, aes(x=NH4.N, y=sc.log.Root.CN)) + geom_point() + geom_smooth(method="lm")
ggplot(plot_bg, aes(x=fungal.rarerich, y=sc.log.Root.CN)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.log.trans.Soil.PC1, y=sc.log.Root.CN)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.Soil.PC2, y=sc.log.Root.CN)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.log.trans.Soil.PC1, y=fungal.rarerich)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.log.trans.Soil.PC1, y=bacteria.rarerich)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.Soil.PC2, y=fungal.rarerich)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.Soil.PC2, y=bacteria.rarerich)) + geom_point() + geom_smooth()

ggplot(plot_bg, aes(x=sc.log.Root.CN, y=fungal.rarerich)) + geom_point() + geom_smooth()
ggplot(plot_bg, aes(x=sc.log.Root.CN, y=bacteria.rarerich)) + geom_point() + geom_smooth()

ggplot(belowground, aes(x=sc.log.Root.CN, y=fungal.rarerich)) + geom_point() + geom_smooth()
ggplot(belowground, aes(x=sc.log.Root.CN, y=bacteria.rarerich)) + geom_point() + geom_smooth()


#get_BLUPs <- function(brm_model){
#  BLUP_50_interval <- tidy(brm_model, par_type="varying", prob=0.5) %>%
#    rename(BLUP=estimate, SD=std.error, lower_50=lower, upper_50=upper)
#  BLUP_95_interval <- tidy(brm_model, par_type="varying", prob=0.95) %>%
#    rename(BLUP=estimate, SD=std.error, lower_95=lower, upper_95=upper)
#  BLUP_df <- left_join(BLUP_50_interval, BLUP_95_interval)
#  return(BLUP_df)
#}

#get_VarComps <- function(brm_model, Distrib_Var){
#  VarComp_50_interval <- tidy(brm_model, par_type="hierarchical", prob=0.5, robust=F) %>%
#    rename(SD_mean=estimate, SD_SD=std.error, lower_50=lower, upper_50=upper)
#  VarComp_95_interval <- tidy(brm_model, par_type="hierarchical", prob=0.95, robust=F) %>%
#    rename(SD_mean=estimate, SD_SD=std.error, lower_95=lower, upper_95=upper)
#  VarComp_df <- left_join(VarComp_50_interval, VarComp_95_interval) %>%
#    mutate(VarComp_mean=SD_mean^2/(sum(SD_mean^2)+Distrib_Var))
#  return(VarComp_df)
#}

#get_FixedEffects <- function(brm_model){
#FE_50_interval <- tidy(brm_model, par_type = "non-varying", prob=0.5, robust=F) %>% rename(lower_50=lower, upper_50=upper)
#FE_95_interval <- tidy(brm_model, par_type = "non-varying", prob=0.95, robust=F) %>% rename(lower_95=lower, upper_95=upper)
#FE_df <- left_join(FE_50_interval, select(FE_95_interval, term, lower_95, upper_95))
#return(FE_df)
#}

composition_plot <- function(composition_data, term){
  require(cowplot)
  
  get.comp <- composition_data %>%
    group_by_(term, "Species") %>%
    summarise_at(vars(Abundance), sum) 
  
  get.mean.abund <- get.comp %>%
    group_by(Species) %>%
    summarise_at(vars(Abundance), mean) %>%
    arrange(desc(Abundance))
  get.comp$Species <- factor(get.comp$Species, levels=unique(get.mean.abund$Species))
  
  ggplot(get.comp, aes(x=Species, y=Abundance)) + geom_bar(stat = "identity") + facet_wrap(as.formula(paste("~", term))) + coord_flip()
}

## ANT-APHID ARTHROPOD COMPOSITION 2012 ANALYSIS ----

aa.comp <- data.frame(select(aa.arth.df, sc.Aphid.treatment, sc.Ant.mound.dist, Genotype, Block, Plot_code), scale(log(select(aa.arth.df, Gracilliaridae_miner:Spider)+1)))
colnames(aa.comp)[c(5,6,7,8,9,16)] <- c("Plot","Gracilliaridae","Tortricidiae", "Cecidomyiidae", "Tenthredinidae", "Formica")

arth.comp.aa.2012.brm <- brm(cbind(Gracilliaridae, Tortricidiae, Cecidomyiidae, Tenthredinidae,
                                   Cicadellidae, Cercopidae, Psyllidae, Orthoptera, Aphididae, Syrphidae, Formica, Spider)~
                               sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot), 
                             data=aa.comp, control=list(adapt_delta=0.99))
summary(arth.comp.aa.2012.brm) # still need to re-run with higher adapt_delta

aa.comp.betas <- posterior_samples(arth.comp.aa.2012.brm, pars = "^b") %>%
  gather(term, betas) %>%
  separate(term, into = c("type", "species", "treatment"), sep = "_", remove=F) %>%
  group_by(term, type, species, treatment) %>%
  summarise(mode = Mode(round(betas,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(betas), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(betas), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(betas), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(betas), prob=0.95)[ ,2]) %>%
  arrange(mode, species, treatment) %>% filter(treatment != "Intercept")
aa.comp.betas$term_ord <- factor(aa.comp.betas$term, levels = aa.comp.betas$term)

ggplot(aa.comp.betas, aes(x=species, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Beta Coefficient") +
  xlab("") +
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~treatment, ncol=1, scales = "free_y")

aa.comp.vars <- posterior_samples(arth.comp.aa.2012.brm, pars = "^sd") %>%
  gather(term, sd) %>%
  separate(term, into = c("type", "treatment", "blank","species", "slope.intercept"), sep = "_", remove=F) %>%
  group_by(term, type, species, treatment, slope.intercept) %>%
  summarise(mode = Mode(round(sd,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(sd), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(sd), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(sd), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(sd), prob=0.95)[ ,2]) %>%
  arrange(mode) # %>% filter(treatment != "Intercept")
aa.comp.vars$term_ord <- factor(aa.comp.vars$term, levels = aa.comp.vars$term)

ggplot(filter(aa.comp.vars, treatment == "Genotype"), aes(x=species, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Standard Deviations") +
  xlab("") +
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_grid(treatment~slope.intercept, scales = "free_y")

ggplot(filter(aa.comp.vars, treatment != "Genotype"), aes(x=species, y=mode)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", shape = 21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Standard Deviations") +
  xlab("") +
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_grid(slope.intercept~treatment, scales = "free_y")

#summary(arth.comp.aa.2012.brm)
#y_aa.arth.comp.2012 <- aa.arth.2012.comp$Occurrence
#yrep_aa.arth.comp.2012 <- posterior_predict(arth.comp.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.comp.aa.2012.brm)
# Aphis farinosa ## CONSIDER A HURDLE POISSON MODEL, THAT WAY I CAN CHECK ESTIMATE PROBABILITY OF THERE BEING AN APHID AND THEN THE EFFECT ON ABUNDANCE
# RERUN 
#hist(filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid")$Abundance)
#arth.Aphis.aa.2012.brm <- general_brm(Abundance~c.Ant.mound.dist+(1+c.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid"), family=poisson(link="log"))
#y_aa.arth.Aphis.2012 <- filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid")$Abundance
#yrep_aa.arth.Aphis.2012 <- posterior_predict(arth.Aphis.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.Aphis.aa.2012.brm)
# Formica obscuripes
# RERUN 
#hist(filter(aa.arth.2012.comp, Species=="ant_F_obscuripes")$Abundance)
#arth.Fobscuripes.aa.2012.brm <- general_brm(Abundance~Aphid.treatment*c.Ant.mound.dist+(1+Aphid.treatment*c.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=filter(aa.arth.2012.comp, Species=="ant_F_obscuripes"), family=poisson(link="log"))
#y_aa.arth.Fobscuripes.2012 <- filter(aa.arth.2012.comp, Species=="ant_F_obscuripes")$Abundance
#yrep_aa.arth.Fobscuripes.2012 <- posterior_predict(arth.Fobscuripes.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.Fobscuripes.aa.2012.brm)


## MAY BE USEFUL... ----
lanphere_BLUPs <- bind_rows(
  mutate(get_BLUPs(trait.PC1.wind.2012.brm), Experiment="Wind", Year="2012", Trait="Trait.PC1"),
  mutate(get_BLUPs(trait.PC2.wind.2012.brm), Experiment="Wind", Year="2012", Trait="Trait.PC2"),
  mutate(get_BLUPs(trait.PC1.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC1"),
  mutate(get_BLUPs(trait.PC2.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC2")
)
write_csv(lanphere_BLUPs, path="output_brms/lanphere_BLUPs.csv")
lanphere_BLUPs <- bind_rows(
  mutate(get_BLUPs(trait.PC1.wind.2012.brm), Experiment="Wind", Year="2012", Response="Trait PC1"),
  mutate(get_BLUPs(trait.PC2.wind.2012.brm), Experiment="Wind", Year="2012", Response="log(Trait PC2 + min(Trait PC2) + 1)"),
  mutate(get_BLUPs(trait.PC1.wind.2013.brm), Experiment="Wind", Year="2013", Response="Trait PC1"),
  mutate(get_BLUPs(trait.PC2.wind.2013.brm), Experiment="Wind", Year="2013", Response="Trait PC2"),
  mutate(get_BLUPs(soil.PC1.wind.brm), Experiment="Wind", Year="2013", Response="log(Soil PC1 + min(Soil PC1) + 1)"),
  mutate(get_BLUPs(soil.PC2.wind.brm), Experiment="Wind", Year="2013", Response="Soil PC2"),
  mutate(get_BLUPs(root_CN.wind.2013.brm), Experiment="Wind", Year="2013", Response="log(Root C:N)"),
  mutate(get_BLUPs(trait.PC1.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Response="Trait PC1"),
  mutate(get_BLUPs(trait.PC2.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Response="Trait PC2"),
  mutate(get_BLUPs(arth.rich.wind.2012.brm), Experiment="Wind", Year="2012", Response="Arthropod Richness"),
  mutate(get_BLUPs(arth.comp.wind.2012.brm), Experiment="Wind", Year="2012", Response="Arthropod Composition"),
  mutate(get_BLUPs(arth.rich.wind.2013.brm), Experiment="Wind", Year="2013", Response="Arthropod Richness"),
  mutate(get_BLUPs(arth.comp.wind.2013.brm), Experiment="Wind", Year="2013", Response="Arthropod Composition"),
  mutate(get_BLUPs(fungal.rarerich.wind.2013.brm), Experiment="Wind", Year="2013", Response="scale(Fungi Rarefied Richness)"),
  mutate(get_BLUPs(bacteria.rarerich.wind.2013.brm), Experiment="Wind", Year="2013", Response="scale(Bacteria Rarefied Richness)"),
  mutate(get_BLUPs(arth.rich.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness"),
  mutate(get_BLUPs(arth.comp.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Composition")
)
write_csv(lanphere_BLUPs, path="output_brms/lanphere_BLUPs.csv")

#lanphere_VarComps <- bind_rows(
#  mutate(get_VarComps(trait.PC1.wind.2012.brm, Distrib_Var=tidy(trait.PC1.wind.2012.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2012", Response="Trait PC1"),
#  mutate(get_VarComps(trait.PC2.wind.2012.brm, Distrib_Var=tidy(trait.PC2.wind.2012.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2012", Response="log(Trait PC2 + min(Trait PC2) + 1)"),
#  mutate(get_VarComps(trait.PC1.wind.2013.brm, Distrib_Var=tidy(trait.PC1.wind.2013.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="Trait PC1"),
#  mutate(get_VarComps(trait.PC2.wind.2013.brm, Distrib_Var=tidy(trait.PC2.wind.2013.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="Trait PC2"),
#  mutate(get_VarComps(soil.PC1.wind.brm, Distrib_Var=tidy(soil.PC1.wind.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="log(Soil PC1 + min(Soil PC1) + 1)"),
#  mutate(get_VarComps(soil.PC2.wind.brm, Distrib_Var=tidy(soil.PC2.wind.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="Soil PC2"),
#  mutate(get_VarComps(root_CN.wind.2013.brm, Distrib_Var=tidy(root_CN.wind.2013.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="log(Root C:N)"),
#  mutate(get_VarComps(trait.PC1.aa.2012.brm, Distrib_Var=tidy(trait.PC1.aa.2012.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Ant-Aphid", Year="2012", Response="Trait PC1"),
#  mutate(get_VarComps(trait.PC2.aa.2012.brm, Distrib_Var=tidy(trait.PC2.aa.2012.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Ant-Aphid", Year="2012", Response="Trait PC2"),
#  mutate(get_VarComps(arth.rich.wind.2012.brm, Distrib_Var=log(1/exp(tidy(arth.rich.wind.2012.brm, parameters = "b_Intercept")$estimate)+1)), 
#         Experiment="Wind", Year="2012", Response="Arthropod Richness"),
#  mutate(get_VarComps(arth.comp.wind.2012.brm, Distrib_Var=pi^2/3), 
#         Experiment="Wind", Year="2012", Response="Arthropod Composition"),
# mutate(get_VarComps(arth.rich.wind.2013.brm, Distrib_Var=log(1/exp(tidy(arth.rich.wind.2013.brm, parameters = "b_Intercept")$estimate)+1)), 
#         Experiment="Wind", Year="2013", Response="Arthropod Richness"),
#  mutate(get_VarComps(arth.comp.wind.2013.brm, Distrib_Var=pi^2/3), 
#         Experiment="Wind", Year="2013", Response="Arthropod Composition"),
#  mutate(get_VarComps(fungal.rarerich.wind.2013.brm, Distrib_Var=tidy(fungal.rarerich.wind.2013.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="scale(Fungi Rarefied Richness)"),
#  mutate(get_VarComps(bacteria.rarerich.wind.2013.brm, Distrib_Var=tidy(bacteria.rarerich.wind.2013.brm, parameters = "sigma")$estimate^2), 
#         Experiment="Wind", Year="2013", Response="scale(Bacteria Rarefied Richness)"),
#  mutate(get_VarComps(arth.rich.aa.2012.brm, Distrib_Var=log(1/exp(tidy(arth.rich.aa.2012.brm, parameters = "b_Intercept")$estimate)+1)), 
#         Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness"),
#  mutate(get_VarComps(arth.comp.aa.2012.brm, Distrib_Var=pi^2/3), 
#         Experiment="Ant-Aphid", Year="2012", Response="Arthropod Composition")
#)
#write_csv(lanphere_VarComps, path="output_brms/lanphere_VarComps.csv")

#lanphere_trait_regs <- bind_rows(
#  mutate(posterior_samples(trait.rich.aa.2012.brm, pars = "^b"), Experiment="Ant-Aphid", Year="2012", Response="scale(log(Arthropod Richness + 1))"),
#  mutate(posterior_samples(trait.rich.wind.2012.brm, pars = "^b"), Experiment="Wind", Year="2012", Response="scale(log(Arthropod Richness + 1))"),
#  mutate(posterior_samples(trait.rich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="scale(log(Arthropod Richness + 1))"),
#  mutate(posterior_samples(trait.rarerich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="scale(Fungi Rarefied Richness)"),
#  mutate(posterior_samples(bact.trait.rarerich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="scale(Bacteria Rarefied Richness)")
#)
#write_csv(lanphere_trait_regs, path="output_brms/lanphere_trait_regs.csv")
