
## LOAD LIBRARIES ----
library(tidyverse)
library(brms)
library(rstan)
library(broom)
library(parallel)


## SET OPTIONS ----
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


## WIND TRAIT DATA ----
w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code = paste(Block, Wind.Exposure, sep="_"))

w.trait.2012 <- filter(w.trait.df, Year=="2012") %>%
  mutate(sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.Trait.PC1 = scale(trait.PC1),
         trans.trait.PC2 = trait.PC2-min(trait.PC2)+1,          # make positive to enable log-transformation
         sc.log.trans.Trait.PC2 = scale(log(trans.trait.PC2)))  # log-transform to normalize prior to scaling

w.trait.2013 <- filter(w.trait.df, Year=="2013") %>%
  mutate(sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.Trait.PC1 = scale(trait.PC1),
         sc.Trait.PC2 = scale(trait.PC2),
         sc.log.Root.CN = scale(log(root_CN)))                  # log-transform to normalize prior to scaling


## ANT-APHID TRAIT DATA ----
aa.trait.df <- read.csv('final_data/ant_aphid_trait_df.csv') %>% tbl_df() %>% mutate(Block=as.factor(Block), Plot_code=paste(Block, Ant.mound.dist, sep="_"))

aa.trait.2012 <- filter(aa.trait.df, Year=="2012") %>%
  mutate(sc.Aphid.treatment = scale(as.numeric(Aphid.treatment)),
         sc.Ant.mound.dist=scale(Ant.mound.dist),
         sc.Trait.PC1 = scale(trait.PC1),
         sc.Trait.PC2 = scale(trait.PC2))


## ANT-APHID ARTHROPOD COMMUNITY DATA ----
aa.arth.df <- read.csv('final_data/ant_aphid_arthropod_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block), 
         sc.Aphid.treatment=scale(as.numeric(Aphid.treatment)),
         sc.Ant.mound.dist=scale(Ant.mound.dist),
         sc.log1.Arthropod.Rich = scale(log(total.rich+1)))          # log(x+1) to normalize prior to scaling


## WIND ARTHROPOD COMMUNITY DATA ----
# (dead plants have already been removed)
wind.arth.df <- read.csv('final_data/wind_arthropod_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code = interaction(Block, Wind.Exposure))          

w.arth.2012 <- wind.arth.df %>% filter(Year == "2012") %>%
  mutate(sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.log1.Arthropod.Rich = scale(log(total.rich+1)))         # log(x+1) to normalize prior to scaling

w.arth.2013 <- wind.arth.df %>% filter(Year == "2013") %>%
  mutate(sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.log1.Arthropod.Rich = scale(log(total.rich+1)))         # log(x+1) to normalize prior to scaling


## WIND FUNGAL COMMUNITY DATA ----
fungal.df <- read.csv("final_data/fungal.df.csv") %>% tbl_df() %>% 
  mutate(Block = as.factor(Block),
         sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.Fungi.Rarerich = scale(fungal.rarerich))


## WIND BACTERIA COMMUNITY DATA ----
bacteria.df <- read.csv("final_data/bacteria.df.csv") %>% tbl_df() %>% 
  mutate(Block = as.factor(Block), 
         sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         sc.Bacteria.Rarerich = scale(bacteria.rarerich))


## WIND SOIL DATA ----
w.soil <- read.csv('final_data/wind_soil_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block),
         sc.Wind.Exposure = scale(as.numeric(Wind.Exposure)),
         soil.PC1.trans = soil.PC1-min(soil.PC1)+1,             # make positive to enable log-transformation
         sc.log.trans.Soil.PC1 = scale(log(soil.PC1.trans)),
         sc.Soil.PC2 = scale(soil.PC2)) 


## FUNCTIONS FOR ANALYSIS ----
general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=gaussian(link="identity"), 
      prior=c(prior(normal(0,1), class=b),
              prior(normal(0,1), class=sd)),
      control=list(adapt_delta=0.99),
      chains=4)
  # all other brm parameters correspond to the defaults
}

posterior_SDs <- function(brm_model, df, FE_formula){
  Fixed_Effects <- posterior_samples(brm_model, pars = "^b")
  sample_size <- dim(Fixed_Effects)[1]
  
  get_model_matrix <- model.matrix(as.formula(FE_formula), data=df)
  
  model_size <- dim(get_model_matrix)[2]
  
  FE_SD_list <- list()
  for(i in 1:model_size){
    FE_SD_vector <- c()
    for(j in 1:sample_size){
      FE_SD_vector[j] <- sd(as.vector(Fixed_Effects[j,i] * t(get_model_matrix)[i, ]))
    }
    FE_SD_list[[i]] <- FE_SD_vector
  }
  FE_SD_df <- as.data.frame(FE_SD_list)
  colnames(FE_SD_df) <- gsub("b_", "sd_", colnames(Fixed_Effects))
  
  SD_df <- data.frame(sample = 1:4000,
                      FE_SD_df,                                         # fixed effects
                      posterior_samples(brm_model, pars = "^sd"),       # random effects
                      posterior_samples(brm_model, pars = "sigma"))     # residual variance  
  return(SD_df)
}


## WIND TRAIT PC1 2012 ANALYSIS ----

trait.PC1.wind.2012.brm <- general_brm(sc.Trait.PC1~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(trait.PC1.wind.2012.brm)

## WIND TRAIT PC2 2012 ANALYSIS ----

trait.PC2.wind.2012.brm <- brm(sc.log.trans.Trait.PC2~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(trait.PC2.wind.2012.brm) 

## WIND 2012 - SEPARATE TRAIT ANALYSES ----

trichomes.w.2012 <- general_brm(scale(log(leaf_trichome.density+1))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(trichomes.w.2012) 

leafWC.w.2012 <- general_brm(scale(log(leaf_WC))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(leafWC.w.2012) 

height.w.2012 <- general_brm(scale(Height)~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(height.w.2012) 

shoot.count.w.2012 <- general_brm(scale(log(all.shoot.count))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(shoot.count.w.2012) 

shoot.length.w.2012 <- general_brm(scale(all.shoot.avg.length)~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012)
summary(shoot.length.w.2012) 

## ANT-APHID TRAIT PC1 2012 ANALYSIS ----

trait.PC1.aa.2012.brm <- general_brm(sc.Trait.PC1~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(trait.PC1.aa.2012.brm) 

## ANT-APHID TRAIT PC2 2012 ANALYSIS ----

trait.PC2.aa.2012.brm <- general_brm(sc.Trait.PC2~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(trait.PC2.aa.2012.brm) 

## ANT-APHID 2012 - SEPARATE TRAIT ANALYSES ----
hist(log(aa.trait.2012$mature.shoot.avg.length))

trichomes.aa.2012 <- general_brm(scale(log(leaf_trichome.density+1))~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(trichomes.aa.2012) 

leafWC.aa.2012 <- general_brm(scale(log(leaf_WC))~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(leafWC.aa.2012) 

height.aa.2012 <- general_brm(scale(Height)~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(height.aa.2012) 

shoot.count.aa.2012 <- general_brm(scale(log(all.shoot.count))~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(shoot.count.aa.2012) 

shoot.length.aa.2012 <- general_brm(scale(log(mature.shoot.avg.length+1))~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012)
summary(shoot.length.aa.2012) 

## WIND TRAIT PC1 2013 ANALYSIS ----

trait.PC1.wind.2013.brm <- general_brm(sc.Trait.PC1~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(trait.PC1.wind.2013.brm)

## WIND TRAIT PC2 2013 ANALYSIS ----

trait.PC2.wind.2013.brm <- general_brm(sc.Trait.PC2~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(trait.PC2.wind.2013.brm) 

## WIND ROOT C:N 2013 ANALYSIS ----

root_CN.wind.2013.brm <- general_brm(sc.log.Root.CN~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013) 
summary(root_CN.wind.2013.brm) 

belowground <- left_join(transmute(w.trait.2013, Block, sc.Wind.Exposure=scale(as.numeric(Wind.Exposure)), Genotype, plant_ID, Plot_code=paste(Block,ifelse(Wind.Exposure=="Exposed","E","U"),sep="."), sc.log.Root.CN = as.numeric(sc.log.Root.CN), root_C.perc, root_N.perc),
                         transmute(w.soil, Plot_code, sc.log.trans.Soil.PC1=as.numeric(sc.log.trans.Soil.PC1), sc.Soil.PC2=as.numeric(sc.Soil.PC2), Total.N, NO3.N, NH4.N)) # %>%
 # left_join(., transmute(fungal.df, fungal_rarerich))

root_CN.wind.2013 <- general_brm(sc.log.Root.CN~sc.Wind.Exposure+sc.log.trans.Soil.PC1+sc.Soil.PC2+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=belowground) 
summary(root_CN.wind.2013)

## WIND 2013 - SEPARATE TRAIT ANALYSES ----

leafCN.w.2013 <- general_brm(scale(log(leaf_C_N))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(leafCN.w.2013) 

leafWC.w.2013 <- general_brm(scale(log(leaf_WC))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(leafWC.w.2013) 

SLA.w.2013 <- general_brm(scale(log(SLA))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(SLA.w.2013) 

height.w.2013 <- general_brm(scale(Height)~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(height.w.2013) 

shoot.count.w.2013 <- general_brm(scale(log(all.shoot.count))~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(shoot.count.w.2013) 

shoot.length.w.2013 <- general_brm(scale(all.shoot.avg.length)~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013)
summary(shoot.length.w.2013) 

## WIND ARTHROPOD RICHNESS 2012 ANALYSIS ----

arth.rich.wind.2012.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2012)
summary(arth.rich.wind.2012.brm)

## ANT-APHID ARTHROPOD RICHNESS 2012 ANALYSIS ----

arth.rich.aa.2012.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.arth.df)
summary(arth.rich.aa.2012.brm)

## ANT-APHID APHIS FARINOSA 2012 ANALYSIS ----

aphis.aa.2012.brm <- general_brm(scale(log(aphid_Aphis+1))~sc.Ant.mound.dist+(1+sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.arth.df)
summary(aphis.aa.2012.brm)

## ANT-APHID FORMICA OBSCURIPES 2012 ANALYSIS ----

formica.aa.2012.brm <- general_brm(scale(log(ant_F_obscuripes+1))~sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.arth.df)
summary(formica.aa.2012.brm)

## WIND ARTHROPOD RICHNESS 2013 ANALYSIS ----

arth.rich.wind.2013.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2013)
summary(arth.rich.wind.2013.brm)

## WIND FUNGI RAREFIED-RICHNESS 2013 ANALYSIS ----

fungal.rarerich.wind.2013.brm <- general_brm(sc.Fungi.Rarerich~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=fungal.df)
summary(fungal.rarerich.wind.2013.brm)

## WIND BACTERIAL RAREFIED-RICHNESS 2013 ANALYSIS ----

bacteria.rarerich.wind.2013.brm <- general_brm(sc.Bacteria.Rarerich~sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=bacteria.df)
summary(bacteria.rarerich.wind.2013.brm)

## ANT-APHID TRAIT-ARTHROPOD 2012 ANALYSIS ----
aa.trait.arth.2012 <- left_join(aa.arth.df, select(aa.trait.2012, plant_ID, sc.Trait.PC1, sc.Trait.PC2), by="plant_ID") %>%
  mutate(fact.Ant.mound.dist = factor(Ant.mound.dist))

trait.rich.aa.2012.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Trait.PC1+sc.Trait.PC2+sc.Aphid.treatment*sc.Ant.mound.dist+(1+sc.Aphid.treatment*sc.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.arth.2012)
summary(trait.rich.aa.2012.brm) # trait.PC1 is key driver, but there is also a weak effect of trait.PC2
trait.rich.aa.12 <- posterior_samples(trait.rich.aa.2012.brm, pars = "^b") # for plotting

plot(marginal_effects(trait.rich.aa.2012.brm, effects = "sc.Trait.PC1"), points=T)
plot(marginal_effects(trait.rich.aa.2012.brm, effects = "sc.Trait.PC2"), points=T)


## WIND TRAIT-ARTHROPOD 2012 ANALYSIS ----
w.trait.arth.2012 <- left_join(w.arth.2012, select(w.trait.2012, plant_ID, sc.Trait.PC1, sc.Trait.PC2 = sc.log.trans.Trait.PC2), by="plant_ID") 

trait.rich.wind.2012.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Trait.PC1+sc.Trait.PC2+sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.arth.2012)
summary(trait.rich.wind.2012.brm) # trait.PC1 is primary effect
trait.rich.wind.12 <- posterior_samples(trait.rich.wind.2012.brm, pars = "^b") # for plotting

plot(marginal_effects(trait.rich.wind.2012.brm, effects = "sc.Trait.PC1"), points=T)


## WIND TRAIT-ARTHROPOD 2013 ANALYSIS ----
w.trait.arth.2013 <- left_join(w.arth.2013, select(w.trait.2013, plant_ID, sc.Trait.PC1, sc.Trait.PC2), by="plant_ID") 

trait.rich.wind.2013.brm <- general_brm(sc.log1.Arthropod.Rich~sc.Trait.PC1+sc.Trait.PC2+sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.arth.2013)
summary(trait.rich.wind.2013.brm) # trait.PC1 is the primary effect
trait.rich.wind.13 <- posterior_samples(trait.rich.wind.2013.brm, pars = "^b") # for plotting

plot(marginal_effects(trait.rich.wind.2013.brm, effects = "sc.Trait.PC1"), points=T)


## WIND TRAIT/SOIL-FUNGAL 2013 ANALYSIS ----
w.trait.fung.2013 <- left_join(fungal.df, select(w.trait.2013, plant_ID, sc.Trait.PC1, sc.Trait.PC2, sc.log.Root.CN), by="plant_ID") %>%
  left_join(., select(w.soil, sc.Soil.PC1 = sc.log.trans.Soil.PC1, sc.Soil.PC2, Plot_code)) 

trait.rarerich.wind.2013.brm <- general_brm(sc.Fungi.Rarerich~sc.log.Root.CN+sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.fung.2013) # sc.log.trans.Soil.PC1+sc.Soil.PC2+
summary(trait.rarerich.wind.2013.brm) # log_root_CN 
trait.rarerich.wind.13 <- posterior_samples(trait.rarerich.wind.2013.brm, pars = "^b") # for plotting

plot(marginal_effects(trait.rarerich.wind.2013.brm, effects = "sc.log.Root.CN"), points=T)


## WIND TRAIT/SOIL-BACTERIA 2013 ANALYSIS ----
w.trait.bact.2013 <- left_join(bacteria.df, select(w.trait.2013, plant_ID, sc.Trait.PC1, sc.Trait.PC2, sc.log.Root.CN), by="plant_ID") %>%
  left_join(., select(w.soil, sc.Soil.PC1 = sc.log.trans.Soil.PC1, sc.Soil.PC2, Plot_code)) 

bact.trait.rarerich.wind.2013.brm <- general_brm(sc.Bacteria.Rarerich~sc.log.Root.CN+sc.Wind.Exposure+(1+sc.Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.bact.2013) # sc.log.trans.Soil.PC1+sc.Soil.PC2+
summary(bact.trait.rarerich.wind.2013.brm) # strong effect of soil.PC2
bact.trait.rarerich.wind.13 <- posterior_samples(bact.trait.rarerich.wind.2013.brm, pars = "^b") # for plotting


## TIDY AND SAVE OUTPUT ----
wind_SDs <- bind_rows(
  mutate(posterior_SDs(trait.PC1.wind.2012.brm, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Trait PC1"), #Response="scale(Trait PC1)"),
  mutate(posterior_SDs(trait.PC2.wind.2012.brm, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Trait PC2"), #Response="scale(log(Trait PC2 + min(Trait PC2) + 1))"),
  mutate(posterior_SDs(trichomes.w.2012, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Trichome Density"),
  mutate(posterior_SDs(leafWC.w.2012, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Leaf Water Content"),
  mutate(posterior_SDs(height.w.2012, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Plant Height"),
  mutate(posterior_SDs(shoot.count.w.2012, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Shoot Count"),
  mutate(posterior_SDs(shoot.length.w.2012, df=w.trait.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Shoot Length"),
  mutate(posterior_SDs(trait.PC1.wind.2013.brm, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Trait PC1"), #Response="scale(Trait PC1)"),
  mutate(posterior_SDs(trait.PC2.wind.2013.brm, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Trait PC2"), #Response="scale(Trait PC2)"),
  mutate(posterior_SDs(leafCN.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Leaf C:N"),
  mutate(posterior_SDs(leafWC.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Leaf Water Content"),
  mutate(posterior_SDs(SLA.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="SLA"),
  mutate(posterior_SDs(height.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Plant Height"),
  mutate(posterior_SDs(shoot.count.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Shoot Count"),
  mutate(posterior_SDs(shoot.length.w.2013, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Shoot Length"),
  #mutate(posterior_SDs(soil.PC1.wind.brm, df=w.soil, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Soil PC1"), #Response="scale(log(Soil PC1 + min(Soil PC1) + 1))"), 
  #mutate(posterior_SDs(soil.PC2.wind.brm, df=w.soil, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Soil PC2"), #Response="scale(Soil PC2)"),
  mutate(posterior_SDs(root_CN.wind.2013.brm, df=w.trait.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Root C:N"), #Response="scale(log(Root C:N))"),
  mutate(posterior_SDs(arth.rich.wind.2012.brm, df=w.arth.2012, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))"),
  mutate(posterior_SDs(arth.rich.wind.2013.brm, df=w.arth.2013, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))"),
  mutate(posterior_SDs(fungal.rarerich.wind.2013.brm, df=fungal.df, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Fungi Rarefied Richness"), #Response="scale(Fungi Rarefied Richness)"),
  mutate(posterior_SDs(bacteria.rarerich.wind.2013.brm, df=bacteria.df, FE_formula = "~sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Bacteria Rarefied Richness") #Response="scale(Bacteria Rarefied Richness)")
)
write_csv(wind_SDs, path="output_brms/wind_SDs.csv")

ant.aphid_SDs <- bind_rows(
  mutate(posterior_SDs(trait.PC1.aa.2012.brm, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Trait PC1"), #Response="scale(Trait PC1)"),
  mutate(posterior_SDs(trait.PC2.aa.2012.brm, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Trait PC2"), #Response="scale(Trait PC2)"),
  mutate(posterior_SDs(trichomes.aa.2012, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Trichome Density"),
  mutate(posterior_SDs(leafWC.aa.2012, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Leaf Water Content"),
  mutate(posterior_SDs(height.aa.2012, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Plant Height"),
  mutate(posterior_SDs(shoot.count.aa.2012, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Shoot Count"),
  mutate(posterior_SDs(shoot.length.aa.2012, df=aa.trait.2012, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Shoot Length"),
  mutate(posterior_SDs(arth.rich.aa.2012.brm, df=aa.arth.df, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))")
  mutate(posterior_SDs(aphis.aa.2012.brm, df=aa.arth.df, FE_formula = "~sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Aphis farinosa"),
  mutate(posterior_SDs(formica.aa.2012.brm, df=aa.arth.df, FE_formula = "~sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Formica obscuripes")
)
write_csv(ant.aphid_SDs, path="output_brms/ant.aphid_SDs.csv")

lanphere_trait_regs <- bind_rows(
  mutate(posterior_SDs(trait.rich.aa.2012.brm, df=aa.trait.arth.2012, FE_formula = "~sc.Trait.PC1+sc.Trait.PC2+sc.Aphid.treatment*sc.Ant.mound.dist"), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))"),
  mutate(posterior_SDs(trait.rich.wind.2012.brm, df=w.trait.arth.2012, FE_formula = "~sc.Trait.PC1+sc.Trait.PC2+sc.Wind.Exposure"), Experiment="Wind", Year="2012", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))"),
  mutate(posterior_SDs(trait.rich.wind.2013.brm, df=w.trait.arth.2013, FE_formula = "~sc.Trait.PC1+sc.Trait.PC2+sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Arthropod Richness"), #Response="scale(log(Arthropod Richness + 1))"),
  mutate(posterior_SDs(trait.rarerich.wind.2013.brm, df=w.trait.fung.2013, FE_formula = "~sc.log.Root.CN+sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Fungi Rarefied Richness"), #Response="scale(Fungi Rarefied Richness)"),
  mutate(posterior_SDs(bact.trait.rarerich.wind.2013.brm, df=w.trait.bact.2013, FE_formula = "~sc.log.Root.CN+sc.Wind.Exposure"), Experiment="Wind", Year="2013", Response="Bacteria Rarefied Richness") #Response="scale(Bacteria Rarefied Richness)")
)
write_csv(lanphere_trait_regs, path="output_brms/lanphere_trait_regs.csv")


## PLOT TRAIT EFFECTS
library(cowplot)
library(coda)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

trait.effects.df <- bind_rows(mutate(gather(select(trait.rich.aa.12, b_sc.Trait.PC1, b_sc.Trait.PC2), Trait, Effect_Size), Response = "Arthropod Richness", Experiment_Year = "Ant-Aphid 2012"), 
          mutate(gather(select(trait.rich.wind.12, b_sc.Trait.PC1, b_sc.Trait.PC2), Trait, Effect_Size), Response = "Arthropod Richness", Experiment_Year = "Wind 2012"), 
          mutate(gather(select(trait.rich.wind.13, b_sc.Trait.PC1, b_sc.Trait.PC2), Trait, Effect_Size), Response = "Arthropod Richness", Experiment_Year = "Wind 2013"), 
          mutate(gather(select(trait.rarerich.wind.13, b_sc.log.Root.CN), Trait, Effect_Size), Response = "Fungi Rarefied Richness", Experiment_Year = "Wind 2013"), 
          mutate(gather(select(bact.trait.rarerich.wind.13, b_sc.log.Root.CN), Trait, Effect_Size), Response = "Bacteria Rarefied Richness", Experiment_Year = "Wind 2013")) %>%
  separate(Experiment_Year, into=c("Experiment","Year"), sep = " ", remove=F) %>%
  group_by(Experiment_Year, Experiment, Year, Response, Trait) %>%
  summarise(mode = Mode(round(Effect_Size,2)), 
            HPDI_lower_50 = HPDinterval(as.mcmc(Effect_Size), prob=0.5)[ ,1],
            HPDI_upper_50 = HPDinterval(as.mcmc(Effect_Size), prob=0.5)[ ,2],
            HPDI_lower_95 = HPDinterval(as.mcmc(Effect_Size), prob=0.95)[ ,1],
            HPDI_upper_95 = HPDinterval(as.mcmc(Effect_Size), prob=0.95)[ ,2])
trait.effects.df$Trait <- factor(trait.effects.df$Trait, 
                                 levels = rev(c("b_sc.Trait.PC1","b_sc.Trait.PC2","b_sc.log.Root.CN")),
                                 labels = rev(c("Trait PC1", "Trait PC2", "Root C:N")))

plot_trait.effects <- ggplot(trait.effects.df, aes(x = Trait, y=mode, shape=Response)) +
  geom_linerange(aes(ymin=HPDI_lower_95, ymax=HPDI_upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=HPDI_lower_50, ymax=HPDI_upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=3.5, fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() + 
  scale_shape_manual(values = c(21,23,24)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  xlab("") +
  ylab("Standardized Effect Size") +
  facet_wrap(~Experiment_Year, ncol=1, scales = "free_y")
save_plot(filename = "fig_trait_effects.png", plot = plot_trait.effects, base_height = 6, base_width = 8.5)


