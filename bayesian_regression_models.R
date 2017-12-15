library(tidyverse)
library(brms)
library(broom)

## LOAD PLANT TRAIT DATA ----

w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code=paste(Block, Wind.Exposure, sep="_"))
w.trait.2012 <- filter(w.trait.df, Year=="2012")
w.trait.2013 <- filter(w.trait.df, Year=="2013")

aa.trait.df <- read.csv('final_data/ant_aphid_trait_df.csv') %>% tbl_df() %>% mutate(Block=as.factor(Block), Ant.mound.dist=as.factor(Ant.mound.dist), Plot_code=paste(Block, Ant.mound.dist, sep="_"))
aa.trait.2012 <- filter(aa.trait.df, Year=="2012")

## upload datasets ----

## ant-aphid: aphid growth rates
aa.aphid.GR <- read.csv('final_data/ant_aphid_Aphis_popgrowth_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) %>%
  select(Date_rel:plant_ID, Aphis.growth.rate)
glimpse(aa.aphid.GR)

## ant-aphid: arthropod community
aa.arth.df <- read.csv('final_data/ant_aphid_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(X = as.factor(X),
         Block = as.factor(Block),
         fact.Ant.mound.dist = as.factor(Ant.mound.dist),
         ord.Ant.mound.dist = ordered(Ant.mound.dist),
         GxE = C(interaction(Genotype,Aphid.treatment), "contr.sum")) #%>% #,
#Plot_code = paste(Block, fact.Ant.mound.dist, sep = "_")) %>%

aa.arth.names <- colnames(select(aa.arth.df, Gracilliaridae_miner:Spider))

aa.arth.2012 <- aa.arth.df

# subset of data where plants had at least one arthropod individual
aa.arth.12.pos <- aa.arth.df %>%
  filter(total.abund > 1) 


## wind: arthropod community
# dead plants have already been removed
wind.arth.df <- read.csv('final_data/wind_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(X = as.factor(X),
         Block = as.factor(Block),
         Year = as.factor(Year),
         Plot_code = interaction(Block, Wind.Exposure),
         GxE = C(interaction(Genotype, Wind.Exposure), contr = "contr.sum", how.many = 9))  # create a new variable for the interaction to permit testing of main effects with type 3 sum of squares. 
#contrasts(wind.arth.df$Wind.Exposure) <- "contr.sum" # important for calculating variance explained
wind.arth.names <- colnames(select(wind.arth.df, Gracilliaridae_miner:Spider)) # for subsetting community data
length(wind.arth.names)

w.arth.2012 <- wind.arth.df %>% filter(Year == "2012")
w.arth.2013 <- wind.arth.df %>% filter(Year == "2013")

# 2012 dataset
# focus dataset on aggregated Family/Order arthropod groupings
w.arth.12.full <- wind.arth.df %>%
  filter(Year == "2012") %>%
  select(Block:plant_ID, Plot_code, GxE, total.abund,
         Gracilliaridae_miner:Spider) 

# subset of data where plants had at least one arthropod individual
w.arth.12.pos2 <- w.arth.12.full %>%
  filter(total.abund > 1) 


# 2013 dataset
# same structure as 2012 dataset
w.arth.13.full <- wind.arth.df %>%
  filter(Year == "2013") %>%
  select(Block:plant_ID, Plot_code, GxE, total.abund,
         Gracilliaridae_miner:Spider) 

w.arth.13.pos2 <- w.arth.13.full %>%
  filter(total.abund > 1)

## Fungal community

f.taxa <- read.csv("final_data/fungi_taxa_table.csv") %>% tbl_df() %>% rename(OTU_ID = X)

fungal.df <- read.csv("final_data/fungal.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich)))

length(f.OTUs)
length(which(colSums(fungal.df[ ,f.OTUs]>0)>10))

## Bacterial community

bacteria.df <- read.csv("final_data/bacteria.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich)))

length(b.OTUs)
length(which(colSums(bacteria.df[ ,b.OTUs]>0)>10))

## wind soil
w.soil <- read.csv('final_data/wind_soil_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) # necessary for modelling as a random effect model


## FUNCTIONS FOR ANALYSIS ----

general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=family, 
      prior=prior(normal(0,1), class=sd),
      control=list(adapt_delta=0.99),
      chains=1)
  # all other brm parameters correspond to the defaults
}

get_BLUPs <- function(brm_model){
  BLUP_50_interval <- tidy(brm_model, par_type="varying", prob=0.5) %>%
    rename(BLUP=estimate, SD=std.error, lower_50=lower, upper_50=upper)
  BLUP_95_interval <- tidy(brm_model, par_type="varying", prob=0.95) %>%
    rename(BLUP=estimate, SD=std.error, lower_95=lower, upper_95=upper)
  BLUP_df <- left_join(BLUP_50_interval, BLUP_95_interval)
  return(BLUP_df)
}

get_VarComps <- function(brm_model, Distrib_Var){
  VarComp_50_interval <- tidy(brm_model, par_type="hierarchical", prob=0.5, robust=F) %>%
    rename(SD_mean=estimate, SD_SD=std.error, lower_50=lower, upper_50=upper)
  VarComp_95_interval <- tidy(brm_model, par_type="hierarchical", prob=0.95, robust=F) %>%
    rename(SD_mean=estimate, SD_SD=std.error, lower_95=lower, upper_95=upper)
  VarComp_df <- left_join(VarComp_50_interval, VarComp_95_interval) %>%
    mutate(VarComp_mean=SD_mean^2/(sum(SD_mean^2)+Distrib_Var))
  return(VarComp_df)
}

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


## WIND TRAIT PC1 2012 ANALYSIS ----
hist(w.trait.2012$trait.PC1)
trait.PC1.wind.2012.brm <- general_brm(trait.PC1~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.trait.2012, family=gaussian(link="identity"))

y_w.trait.PC1.2012 <- w.trait.2012$trait.PC1
yrep_w.trait.PC1.2012 <- posterior_predict(trait.PC1.wind.2012.brm, nsamples=100)
#launch_shinystan(trait.PC1.wind.2012.brm)


## WIND TRAIT PC2 2012 ANALYSIS ----
hist(w.trait.2012$trait.PC2)
w.trait.2012$trait.PC2.trans <- w.trait.2012$trait.PC2-min(w.trait.2012$trait.PC2)+1 # make so minimum value is 1
hist(log(w.trait.2012$trait.PC2.trans))
trait.PC2.wind.2012.brm <- general_brm(log(trait.PC2.trans)~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.trait.2012, family=gaussian(link="identity"))

y_w.trait.PC2.2012 <- log(w.trait.2012$trait.PC2.trans)
yrep_w.trait.PC2.2012 <- posterior_predict(trait.PC2.wind.2012.brm, nsamples=100)
#launch_shinystan(trait.PC2.wind.2012.brm)


## ANT-APHID TRAIT PC1 2012 ANALYSIS ----
hist(aa.trait.2012$trait.PC1)
trait.PC1.aa.2012.brm <- general_brm(trait.PC1~(1|Genotype*Aphid.treatment*Ant.mound.dist)+(1|Block)+(1|Plot_code), data=aa.trait.2012, family=gaussian(link="identity"))

y_aa.trait.PC1.2012 <- aa.trait.2012$trait.PC1
yrep_aa.trait.PC1.2012 <- posterior_predict(trait.PC1.aa.2012.brm, nsamples=100)
#launch_shinystan(trait.PC1.aa.2012.brm)


## ANT-APHID TRAIT PC2 2012 ANALYSIS ----
hist(aa.trait.2012$trait.PC2)
trait.PC2.aa.2012.brm <- general_brm(trait.PC2~(1|Genotype*Aphid.treatment*Ant.mound.dist)+(1|Block)+(1|Plot_code), data=aa.trait.2012, family=gaussian(link="identity"))

y_aa.trait.PC2.2012 <- aa.trait.2012$trait.PC2
yrep_aa.trait.PC2.2012 <- posterior_predict(trait.PC2.aa.2012.brm, nsamples=100)
#launch_shinystan(trait.PC2.aa.2012.brm)


## WIND TRAIT PC1 2013 ANALYSIS ----
hist(w.trait.2013$trait.PC1)
trait.PC1.wind.2013.brm <- general_brm(trait.PC1~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.trait.2013, family=gaussian(link="identity"))

y_w.trait.PC1.2013 <- w.trait.2013$trait.PC1
yrep_w.trait.PC1.2013 <- posterior_predict(trait.PC1.wind.2013.brm, nsamples=100)
#launch_shinystan(trait.PC1.wind.2013.brm)


## WIND TRAIT PC2 2013 ANALYSIS ----
hist(w.trait.2013$trait.PC2)
trait.PC2.wind.2013.brm <- general_brm(trait.PC2~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.trait.2013, family=gaussian(link="identity"))

y_w.trait.PC2.2013 <- w.trait.2013$trait.PC2
yrep_w.trait.PC2.2013 <- posterior_predict(trait.PC2.wind.2013.brm, nsamples=100)
#launch_shinystan(trait.PC2.wind.2013.brm)

## WIND ROOT C:N 2013 ANALYSIS ----

## NOTE THAT SOME MINIMUMS ARE BELOW ZERO, WHICH I SHOULD FIX

plot(w.trait.2013$root_N.perc ~ w.trait.2013$root_CN)
hist(log(w.trait.2013$root_N.perc))
hist(w.trait.2013$root_C.perc)

root_Nperc.wind.2013.brm <- general_brm(root_N.perc~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=filter(w.trait.2013, root_N.perc>0), family=gaussian(link="identity"))

y_root_Nperc.2013 <- filter(w.trait.2013, root_N.perc>0)$root_N.perc
yrep_root_Nperc.2013 <- posterior_predict(root_Nperc.wind.2013.brm, nsamples=100)
launch_shinystan(root_Nperc.wind.2013.brm)

## WIND ARTHROPOD RICHNESS 2012 ANALYSIS ----
hist(w.arth.2012$total.rich)
arth.rich.wind.2012.brm <- general_brm(total.rich~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.arth.2012, family=poisson(link="log"))

y_w.arth.rich.2012 <- w.arth.2012$total.rich
yrep_w.arth.rich.2012 <- posterior_predict(arth.rich.wind.2012.brm, nsamples=100)
#launch_shinystan(arth.rich.wind.2012.brm)


## WIND ARTHROPOD COMPOSITION 2012 ANALYSIS ----
w.arth.2012.comp <- select(w.arth.2012, X:spider_Larionoides, Plot_code, GxE) %>%
  gather(key=Species, value=Abundance, ant_F_obscuripes:spider_Larionoides) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

composition_plot(w.arth.2012.comp, term="Wind.Exposure")
composition_plot(w.arth.2012.comp, term="Genotype")
composition_plot(w.arth.2012.comp, term="GxE")

hist(w.arth.2012.comp$Occurrence)
arth.comp.wind.2012.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=w.arth.2012.comp, family=bernoulli(link="logit"))

y_w.arth.comp.2012 <- w.arth.2012.comp$Occurrence
yrep_w.arth.comp.2012 <- posterior_predict(arth.comp.wind.2012.brm, nsamples=100)
#launch_shinystan(arth.comp.wind.2012.brm)


## ANT-APHID ARTHROPOD RICHNESS 2012 ANALYSIS ----
hist(aa.arth.2012$total.rich)
arth.rich.aa.2012.brm <- general_brm(total.rich~(1|Genotype*Aphid.treatment*Ant.mound.dist)+(1|Block)+(1|Plot_code), data=aa.arth.2012, family=poisson(link="log"))

y_aa.arth.rich.2012 <- aa.arth.2012$total.rich
yrep_aa.arth.rich.2012 <- posterior_predict(arth.rich.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.rich.aa.2012.brm)


## ANT-APHID ARTHROPOD COMPOSITION 2012 ANALYSIS ----
aa.arth.2012.comp <- select(aa.arth.2012, X:LTF_Caloptilia) %>%
  gather(key=Species, value=Abundance, aphid_Aphis:LTF_Caloptilia) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

hist(aa.arth.2012.comp$Occurrence)
arth.comp.aa.2012.brm <- general_brm(Occurrence~(1|Genotype*Aphid.treatment*Ant.mound.dist*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=aa.arth.2012.comp, family=bernoulli(link="logit"))

y_aa.arth.comp.2012 <- aa.arth.2012.comp$Occurrence
yrep_aa.arth.comp.2012 <- posterior_predict(arth.comp.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.comp.aa.2012.brm)

# Aphis farinosa ## CONSIDER A HURDLE POISSON MODEL, THAT WAY I CAN CHECK ESTIMATE PROBABILITY OF THERE BEING AN APHID AND THEN THE EFFECT ON ABUNDANCE
hist(filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid")$Abundance)
arth.Aphis.aa.2012.brm <- general_brm(Abundance~(1|Genotype*Ant.mound.dist)+(1|Block)+(1|Plot_code), data=filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid"), family=poisson(link="log"))

y_aa.arth.Aphis.2012 <- filter(aa.arth.2012.comp, Species=="aphid_Aphis", Aphid.treatment=="aphid")$Abundance
yrep_aa.arth.Aphis.2012 <- posterior_predict(arth.Aphis.aa.2012.brm, nsamples=100)
launch_shinystan(arth.Aphis.aa.2012.brm)

# Formica obscuripes
hist(filter(aa.arth.2012.comp, Species=="ant_F_obscuripes")$Abundance)
arth.Fobscuripes.aa.2012.brm <- general_brm(Abundance~(1|Genotype*Aphid.treatment*Ant.mound.dist)+(1|Block)+(1|Plot_code), data=filter(aa.arth.2012.comp, Species=="ant_F_obscuripes"), family=poisson(link="log"))

y_aa.arth.Fobscuripes.2012 <- filter(aa.arth.2012.comp, Species=="ant_F_obscuripes")$Abundance
yrep_aa.arth.Fobscuripes.2012 <- posterior_predict(arth.Fobscuripes.aa.2012.brm, nsamples=100)
launch_shinystan(arth.Fobscuripes.aa.2012.brm)

## WIND ARTHROPOD RICHNESS 2013 ANALYSIS ----
hist(w.arth.2013$total.rich)
arth.rich.wind.2013.brm <- general_brm(total.rich~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=w.arth.2013, family=poisson(link="log"))

y_w.arth.rich.2013 <- w.arth.2013$total.rich
yrep_w.arth.rich.2013 <- posterior_predict(arth.rich.wind.2013.brm, nsamples=100)
#launch_shinystan(arth.rich.wind.2013.brm)


## WIND FUNGAL RAREFIED-RICHNESS 2013 ANALYSIS ----
hist(scale(fungal.df$fungal.rarerich))
fungal.rarerich.wind.2013.brm <- general_brm(scale(fungal.rarerich)~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=fungal.df, family=gaussian(link="identity"))

y_w.fungal.rarerich.2013 <- as.numeric(scale(fungal.df$fungal.rarerich))
yrep_w.fungal.rarerich.2013 <- posterior_predict(fungal.rarerich.wind.2013.brm, nsamples=100)
#launch_shinystan(fungal.rarerich.wind.2013.brm)

## WIND FUNGAL COMPOSITION 2013 ANALYSIS ----
fungal.comp <- select(fungal.df, X:plant_ID, OTU_1347:OTU_713) %>%
  gather(key=Species, value=Abundance, OTU_1347:OTU_713) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

composition_plot(fungal.comp, term="Wind.Exposure")
composition_plot(fungal.comp, term="Genotype")
#composition_plot(fungal.comp, term="GxE")

hist(fungal.comp$Occurrence)
fungal.comp.wind.2013.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=fungal.comp, family=bernoulli(link="logit"))

y_w.fungal.comp.2012 <- fungal.comp$Occurrence
yrep_w.fungal.comp.2012 <- posterior_predict(fungal.comp.wind.2013.brm, nsamples=100)
#launch_shinystan(fungal.comp.wind.2013.brm)

## WIND BACTERIA COMPOSITION 2013 ANALYSIS ----
bacteria.comp <- select(bacteria.df, X:plant_ID, OTU_1347:OTU_713) %>%
  gather(key=Species, value=Abundance, OTU_1347:OTU_713) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

#composition_plot(bacteria.comp, term="Wind.Exposure")
#composition_plot(bacteria.comp, term="Genotype")
#composition_plot(bacteria.comp, term="GxE")

hist(bacteria.comp$Occurrence)
bacteria.comp.wind.2013.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=bacteria.comp, family=bernoulli(link="logit"))

y_w.bacteria.comp.2012 <- bacteria.comp$Occurrence
yrep_w.bacteria.comp.2012 <- posterior_predict(bacteria.comp.wind.2013.brm, nsamples=100)
#launch_shinystan(bacteria.comp.wind.2013.brm)

## WIND BACTERIAL RAREFIED-RICHNESS 2013 ANALYSIS ----
hist(scale(bacteria.df$bacteria.rarerich))
bacteria.rarerich.wind.2013.brm <- general_brm(scale(bacteria.rarerich)~(1|Genotype*Wind.Exposure)+(1|Block)+(1|Plot_code), data=bacteria.df, family=gaussian(link="identity"))

y_w.bacteria.rarerich.2013 <- scale(bacteria.df$bacteria.rarerich)
yrep_w.bacteria.rarerich.2013 <- posterior_predict(bacteria.rarerich.wind.2013.brm, nsamples=100)
#launch_shinystan(bacteria.rarerich.wind.2013.brm)

## WIND ARTHROPOD COMPOSITION 2013 ANALYSIS ----
w.arth.2013.comp <- select(w.arth.2013, X:spider_Larionoides, Plot_code, GxE) %>%
  gather(key=Species, value=Abundance, ant_F_obscuripes:spider_Larionoides) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

composition_plot(w.arth.2013.comp, term="Wind.Exposure")
composition_plot(w.arth.2013.comp, term="Genotype")
composition_plot(w.arth.2013.comp, term="GxE")

hist(w.arth.2013.comp$Occurrence)
arth.comp.wind.2013.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=w.arth.2013.comp, family=bernoulli(link="logit"))

y_w.arth.comp.2013 <- w.arth.2013.comp$Occurrence
yrep_w.arth.comp.2013 <- posterior_predict(arth.comp.wind.2013.brm, nsamples=100)
#launch_shinystan(arth.comp.wind.2013.brm)

## WIND SOIL PC1 ANALYSIS ----
hist(w.soil$soil.PC1)
hist(log(w.soil$soil.PC1-min(w.soil$soil.PC1)+1))
w.soil$soil.PC1.trans <- w.soil$soil.PC1-min(w.soil$soil.PC1)+1
soil.PC1.wind.brm <- general_brm(log(soil.PC1.trans)~(1|Wind.Exposure)+(1|Block), data=w.soil, family=gaussian(link="identity"))

y_w.soil.PC1 <- log(w.soil$soil.PC1.trans)
yrep_w.soil.PC1 <- posterior_predict(soil.PC1.wind.brm, nsamples=100)
#launch_shinystan(soil.PC1.wind.brm)

## WIND SOIL PC2 ANALYSIS ----
hist(w.soil$soil.PC2)
soil.PC2.wind.brm <- general_brm(soil.PC2~(1|Wind.Exposure)+(1|Block), data=w.soil, family=gaussian(link="identity"))

y_w.soil.PC2 <- w.soil$soil.PC2
yrep_w.soil.PC2 <- posterior_predict(soil.PC2.wind.brm, nsamples=100)
#launch_shinystan(soil.PC2.wind.brm)

## TIDY AND SAVE OUTPUT ----
lanphere_trait_BLUPs_2012 <- bind_rows(
  mutate(get_BLUPs(trait.PC1.wind.2012.brm), Experiment="Wind", Year="2012", Trait="Trait.PC1"),
  mutate(get_BLUPs(trait.PC2.wind.2012.brm), Experiment="Wind", Year="2012", Trait="Trait.PC2"),
  mutate(get_BLUPs(trait.PC1.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC1"),
  mutate(get_BLUPs(trait.PC2.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC2")
)
write_csv(lanphere_trait_BLUPs_2012, path="output_brms/lanphere_trait_BLUPs_2012.csv")

lanphere_trait_VarComps_2012 <- bind_rows(
  mutate(get_VarComps(trait.PC1.wind.2012.brm, Distrib_Var=tidy(trait.PC1.wind.2012.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2012", Trait="Trait.PC1"),
  mutate(get_VarComps(trait.PC2.wind.2012.brm, Distrib_Var=tidy(trait.PC2.wind.2012.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2012", Trait="Trait.PC2"),
  mutate(get_VarComps(trait.PC1.wind.2013.brm, Distrib_Var=tidy(trait.PC1.wind.2013.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Trait.PC1"),
  mutate(get_VarComps(trait.PC2.wind.2013.brm, Distrib_Var=tidy(trait.PC2.wind.2013.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Trait.PC2"),
  mutate(get_VarComps(soil.PC1.wind.brm, Distrib_Var=tidy(soil.PC1.wind.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Soil.PC1"),
  mutate(get_VarComps(soil.PC2.wind.brm, Distrib_Var=tidy(soil.PC2.wind.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Soil.PC2"),
  mutate(get_VarComps(root_Nperc.wind.2013.brm, Distrib_Var=tidy(root_Nperc.wind.2013.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Root %N"),
  mutate(get_VarComps(trait.PC1.aa.2012.brm, Distrib_Var=tidy(trait.PC1.aa.2012.brm, parameters = "sigma")$estimate^2), 
         Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC1"),
  mutate(get_VarComps(trait.PC2.aa.2012.brm, Distrib_Var=tidy(trait.PC2.aa.2012.brm, parameters = "sigma")$estimate^2), 
         Experiment="Ant-Aphid", Year="2012", Trait="Trait.PC2")
)
write_csv(lanphere_trait_VarComps_2012, path="output_brms/lanphere_trait_VarComps_2012.csv")

lanphere_arthropods_VarComps_2012 <- bind_rows(
  mutate(get_VarComps(arth.rich.wind.2012.brm, Distrib_Var=log(1/exp(tidy(arth.rich.wind.2012.brm, parameters = "b_Intercept")$estimate)+1)), 
         Experiment="Wind", Year="2012", Trait="Arthropod Richness"),
  mutate(get_VarComps(arth.comp.wind.2012.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2012", Trait="Arthropod Composition"),
  mutate(get_VarComps(arth.rich.wind.2013.brm, Distrib_Var=log(1/exp(tidy(arth.rich.wind.2013.brm, parameters = "b_Intercept")$estimate)+1)), 
         Experiment="Wind", Year="2013", Trait="Arthropod Richness"),
  mutate(get_VarComps(arth.comp.wind.2013.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Arthropod Composition"),
  mutate(get_VarComps(fungal.rarerich.wind.2013.brm, Distrib_Var=tidy(fungal.rarerich.wind.2013.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Fungal Rarefied Richness"),
  mutate(get_VarComps(bacteria.rarerich.wind.2013.brm, Distrib_Var=tidy(bacteria.rarerich.wind.2013.brm, parameters = "sigma")$estimate^2), 
         Experiment="Wind", Year="2013", Trait="Bacterial Rarefied Richness"),
  mutate(get_VarComps(arth.rich.aa.2012.brm, Distrib_Var=log(1/exp(tidy(arth.rich.aa.2012.brm, parameters = "b_Intercept")$estimate)+1)), 
         Experiment="Ant-Aphid", Year="2012", Trait="Arthropod Richness"),
  mutate(get_VarComps(arth.comp.aa.2012.brm, Distrib_Var=pi^2/3), 
         Experiment="Ant-Aphid", Year="2012", Trait="Arthropod Composition"),
  mutate(get_VarComps(fungal.comp.wind.2013.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Fungal Composition"),
  mutate(get_VarComps(bacteria.comp.wind.2013.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Bacterial Composition")
)
write_csv(lanphere_arthropods_VarComps_2012, path="output_brms/lanphere_arthropods_VarComps_2012.csv")