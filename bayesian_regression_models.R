
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
w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block), 
         Wind.Exposure = C(as.factor(Wind.Exposure), "contr.sum"),
         Plot_code=paste(Block, Wind.Exposure, sep="_"))
w.trait.2012 <- filter(w.trait.df, Year=="2012")
w.trait.2013 <- filter(w.trait.df, Year=="2013")


## ANT-APHID TRAIT DATA ----
aa.trait.df <- read.csv('final_data/ant_aphid_trait_df.csv') %>% tbl_df() %>% 
  mutate(Block=as.factor(Block), 
         Aphid.treatment=C(Aphid.treatment, "contr.sum"),        # set sum to zero contrasts
         c.Ant.mound.dist=Ant.mound.dist - mean(Ant.mound.dist), # center
         Plot_code=paste(Block, Ant.mound.dist, sep="_"))
aa.trait.2012 <- filter(aa.trait.df, Year=="2012")


## ANT-APHID ARTHROPOD COMMUNITY DATA ----
aa.arth.df <- read.csv('final_data/ant_aphid_arthropod_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block), 
         Aphid.treatment=C(Aphid.treatment, "contr.sum"),        # set sum to zero contrasts
         c.Ant.mound.dist=Ant.mound.dist - mean(Ant.mound.dist)) # center 

mean.aa.rich.2012 <- mean(aa.arth.df$total.rich)

## WIND ARTHROPOD COMMUNITY DATA ----
# (dead plants have already been removed)
wind.arth.df <- read.csv('final_data/wind_arthropod_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block),
         Wind.Exposure = C(as.factor(Wind.Exposure), "contr.sum"),
         Plot_code = interaction(Block, Wind.Exposure)) 
w.arth.2012 <- wind.arth.df %>% filter(Year == "2012")
w.arth.2013 <- wind.arth.df %>% filter(Year == "2013")

mean.w.rich.2012 <- mean(w.arth.2012$total.rich)
mean.w.rich.2013 <- mean(w.arth.2013$total.rich)

## WIND FUNGAL COMMUNITY DATA ----
fungal.df <- read.csv("final_data/fungal.df.csv") %>% tbl_df() %>% 
  mutate(Block = as.factor(Block),
         Wind.Exposure = C(as.factor(Wind.Exposure), "contr.sum"),
         X = as.factor(X))


## WIND BACTERIA COMMUNITY DATA ----
bacteria.df <- read.csv("final_data/bacteria.df.csv") %>% tbl_df() %>% 
  mutate(Block = as.factor(Block), 
         Wind.Exposure = C(as.factor(Wind.Exposure), "contr.sum"),
         X = as.factor(X))


## WIND SOIL DATA ----
w.soil <- read.csv('final_data/wind_soil_df.csv') %>% tbl_df() %>% 
  mutate(Block = as.factor(Block),
         Wind.Exposure = C(as.factor(Wind.Exposure), "contr.sum")) 


## FUNCTIONS FOR ANALYSIS ----

general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=family, 
      prior=c(prior(normal(0,1), class=b),
              prior(normal(0,1), class=sd)),
      control=list(adapt_delta=0.99, max_treedepth=20),
      chains=4)
  # all other brm parameters correspond to the defaults
}

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

wind_posterior_SDs <- function(brm_model, df, FE_formula="~Wind.Exposure"){
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
  SD_df <- data.frame(sd_Wind.Exposure=FE_SD_list[[2]],                # fixed effects
                      posterior_samples(brm_model, pars = "^sd")) %>%  # random effects
    gather(key=term, value=posterior_SD) 
  return(SD_df)
}

ant.aphid_posterior_SDs <- function(brm_model, df, FE_formula="~Aphid.treatment*c.Ant.mound.dist"){
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
  SD_df <- data.frame(sd_Aphid.treatment=FE_SD_list[[2]], sd_c.Ant.mound.dist=FE_SD_list[[3]], sd_Aphid.x.Ant=FE_SD_list[[4]], # fixed effects
                      posterior_samples(brm_model, pars = "^sd")) %>%                                                          # random effects
    gather(key=term, value=posterior_SD)
  return(SD_df)
}

# still thinking about using coda package to calculate intervals

## WIND TRAIT PC1 2012 ANALYSIS ----
hist(w.trait.2012$trait.PC1)
trait.PC1.wind.2012.brm <- general_brm(trait.PC1~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012, family=gaussian(link="identity"))
summary(trait.PC1.wind.2012.brm)

y_w.trait.PC1.2012 <- w.trait.2012$trait.PC1
yrep_w.trait.PC1.2012 <- posterior_predict(trait.PC1.wind.2012.brm, nsamples=100)
#launch_shinystan(trait.PC1.wind.2012.brm)


## WIND TRAIT PC2 2012 ANALYSIS ----
hist(w.trait.2012$trait.PC2)
w.trait.2012$trait.PC2.trans <- w.trait.2012$trait.PC2-min(w.trait.2012$trait.PC2)+1 # make so minimum value is 1
hist(log(w.trait.2012$trait.PC2.trans))
trait.PC2.wind.2012.brm <- brm(log(trait.PC2.trans)~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2012, family=gaussian(link="identity"))
summary(trait.PC2.wind.2012.brm) 

y_w.trait.PC2.2012 <- log(w.trait.2012$trait.PC2.trans)
yrep_w.trait.PC2.2012 <- posterior_predict(trait.PC2.wind.2012.brm, nsamples=100)
#launch_shinystan(trait.PC2.wind.2012.brm)


## ANT-APHID TRAIT PC1 2012 ANALYSIS ----

hist(aa.trait.2012$trait.PC1)
trait.PC1.aa.2012.brm <- general_brm(trait.PC1~Aphid.treatment*c.Ant.mound.dist+(1+Aphid.treatment*c.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012, family=gaussian(link="identity"))
summary(trait.PC1.aa.2012.brm) 

y_aa.trait.PC1.2012 <- aa.trait.2012$trait.PC1
yrep_aa.trait.PC1.2012 <- posterior_predict(trait.PC1.aa.2012.brm, nsamples=100)
#launch_shinystan(trait.PC1.aa.2012.brm)


## ANT-APHID TRAIT PC2 2012 ANALYSIS ----
hist(aa.trait.2012$trait.PC2)
trait.PC2.aa.2012.brm <- general_brm(trait.PC2~Aphid.treatment*c.Ant.mound.dist+(1+Aphid.treatment*c.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.trait.2012, family=gaussian(link="identity"))
summary(trait.PC2.aa.2012.brm) 

y_aa.trait.PC2.2012 <- aa.trait.2012$trait.PC2
yrep_aa.trait.PC2.2012 <- posterior_predict(trait.PC2.aa.2012.brm, nsamples=100)
#launch_shinystan(trait.PC2.aa.2012.brm)


## WIND TRAIT PC1 2013 ANALYSIS ----
hist(w.trait.2013$trait.PC1)
trait.PC1.wind.2013.brm <- general_brm(trait.PC1~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013, family=gaussian(link="identity"))
summary(trait.PC1.wind.2013.brm)

y_w.trait.PC1.2013 <- w.trait.2013$trait.PC1
yrep_w.trait.PC1.2013 <- posterior_predict(trait.PC1.wind.2013.brm, nsamples=100)
#launch_shinystan(trait.PC1.wind.2013.brm)


## WIND TRAIT PC2 2013 ANALYSIS ----
hist(w.trait.2013$trait.PC2)
trait.PC2.wind.2013.brm <- general_brm(trait.PC2~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.trait.2013, family=gaussian(link="identity"))
summary(trait.PC2.wind.2013.brm) 

y_w.trait.PC2.2013 <- w.trait.2013$trait.PC2
yrep_w.trait.PC2.2013 <- posterior_predict(trait.PC2.wind.2013.brm, nsamples=100)
#launch_shinystan(trait.PC2.wind.2013.brm)

## WIND ROOT C:N 2013 ANALYSIS ----

# note that the model improves if I exclude root_CN > 100, but the results are qualitatively the same even if I retain these data.
hist(log(w.trait.2013$root_CN))
root_CN.wind.2013.brm <- general_brm(log(root_CN)~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=filter(w.trait.2013, root_CN>0), family=gaussian(link="identity")) 
summary(root_CN.wind.2013.brm) 

y_root_CN.2013 <- log(filter(w.trait.2013, root_CN>0)$root_CN)
yrep_root_CN.2013 <- posterior_predict(root_CN.wind.2013.brm, nsamples=100)
#launch_shinystan(root_CN.wind.2013.brm)

## WIND ARTHROPOD RICHNESS 2012 ANALYSIS ----
get_prior(total.rich~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2012, family=poisson(link="log"))
hist(w.arth.2012$total.rich)
arth.rich.wind.2012.brm <- general_brm(total.rich~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2012, family=poisson(link="log"))
summary(arth.rich.wind.2012.brm)

y_w.arth.rich.2012 <- w.arth.2012$total.rich
yrep_w.arth.rich.2012 <- posterior_predict(arth.rich.wind.2012.brm, nsamples=100)
#launch_shinystan(arth.rich.wind.2012.brm)

## WIND ARTHROPOD COMPOSITION 2012 ANALYSIS ----

## RETHINKING THIS BASED ON MULTIVARIATE ANALYSIS, SHOULD BE ABLE TO KEEP IN SAME FORMAT

#w.arth.2012.comp <- select(w.arth.2012, X:spider_Larionoides, Plot_code) %>% #, GxE
#  gather(key=Species, value=Abundance, ant_F_obscuripes:spider_Larionoides) %>%
#  mutate(Occurrence=ifelse(Abundance>0, 1, 0)) %>%
#  select(-Abundance) %>%
#  spread(Species, Occurrence)

#composition_plot(w.arth.2012.comp, term="Wind.Exposure")
#composition_plot(w.arth.2012.comp, term="Genotype")
#composition_plot(w.arth.2012.comp, term="GxE")

#get_prior(cbind(ant_black, ant_F_obscuripes)~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2012.comp, family=bernoulli(link="logit"))

#wind.community <- cbind(ant_black, ant_F_obscuripes, aphid_Aphis, aphid_LG, aphid_Tuberolachnus, caterpillar_LB,
#                        caterpillar_looper, caterpillar_unk, gall_Aculus, gall_Iteomyia, gall_Pontania, gall_R_rigidae,
#                        gall_R_salicisbattatus, gall_R_salicisbrassicoides, grasshopper, leafhopper_camo, leafhopper_C_reductus,
#                        leafhopper_green, leafhopper_nymph_unk, leafhopper_unk, leafhopper_YK, leaftier_Tortricid, 
#                        LTF_Caloptilia, psyllid, red_scale, sawfly_larva, spider_BY, spider_CS, spider_Larionoides, 
#                        spider_NW, spider_Tetragnathid, spider_Theridion, spider_unk, stinkbug, tentmine_Phyllonorycter)
#hist(w.arth.2012.comp$Occurrence)
#arth.comp.wind.2012.brm <- brm(cbind(ant_black, ant_F_obscuripes, aphid_Aphis, aphid_LG, aphid_Tuberolachnus, caterpillar_LB,
#                                     caterpillar_looper, caterpillar_unk, gall_Aculus, gall_Iteomyia, gall_Pontania, gall_R_rigidae,
#                                     gall_R_salicisbattatus, gall_R_salicisbrassicoides, grasshopper, leafhopper_camo, leafhopper_C_reductus,
#                                     leafhopper_green, leafhopper_nymph_unk, leafhopper_unk, leafhopper_YK, leaftier_Tortricid, 
#                                     LTF_Caloptilia, psyllid, red_scale, sawfly_larva, spider_BY, spider_CS, spider_Larionoides, 
#                                     spider_NW, spider_Tetragnathid, spider_Theridion, spider_unk, stinkbug, tentmine_Phyllonorycter)
#                               ~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2012.comp, family=bernoulli(link="logit"), control=list(adapt_delta=0.95, max_treedepth=20))
#summary(arth.comp.wind.2012.brm)
#round(fixef(arth.comp.wind.2012.brm),1)
#tidy(arth.comp.wind.2012.brm, par_type = "hierarchical")

#y_w.arth.comp.2012 <- w.arth.2012.comp$Occurrence
#yrep_w.arth.comp.2012 <- posterior_predict(arth.comp.wind.2012.brm, nsamples=100)
#launch_shinystan(arth.comp.wind.2012.brm)


## ANT-APHID ARTHROPOD RICHNESS 2012 ANALYSIS ----

# unable to fit the most complex random effects structure to the model
hist(aa.arth.df$total.rich)
arth.rich.aa.2012.brm <- general_brm(total.rich~Aphid.treatment*c.Ant.mound.dist+(1+Aphid.treatment|Genotype)+(1|Block)+(1|Plot_code), data=aa.arth.df, family=poisson(link="log"))
summary(arth.rich.aa.2012.brm)

y_aa.arth.rich.2012 <- aa.arth.df$total.rich
yrep_aa.arth.rich.2012 <- posterior_predict(arth.rich.aa.2012.brm, nsamples=100)
#launch_shinystan(arth.rich.aa.2012.brm)


## ANT-APHID ARTHROPOD COMPOSITION 2012 ANALYSIS ----
#aa.arth.2012.comp <- select(aa.arth.df, X:LTF_Caloptilia) %>%
#  gather(key=Species, value=Abundance, aphid_Aphis:LTF_Caloptilia) %>%
#  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

## EVALUATE WHETHER CBIND FUNCTION WILL WORK BEST

#hist(aa.arth.2012.comp$Occurrence)
#arth.comp.aa.2012.brm <- general_brm(Occurrence~Aphid.treatment*c.Ant.mound.dist+(1+Aphid.treatment*c.Ant.mound.dist|Genotype)+(1|Block)+(1|Plot_code), data=aa.arth.2012, family=bernoulli(link="logit"))
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

## WIND ARTHROPOD RICHNESS 2013 ANALYSIS ----

hist(w.arth.2013$total.rich)
arth.rich.wind.2013.brm <- general_brm(total.rich~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=w.arth.2013, family=poisson(link="log"))
summary(arth.rich.wind.2013.brm)

y_w.arth.rich.2013 <- w.arth.2013$total.rich
yrep_w.arth.rich.2013 <- posterior_predict(arth.rich.wind.2013.brm, nsamples=100)
#launch_shinystan(arth.rich.wind.2013.brm)


## WIND FUNGAL RAREFIED-RICHNESS 2013 ANALYSIS ----

# note that scaling the response variable makes my general priors (normal(mean=0, sd=1)) appropriate.
hist(scale(fungal.df$fungal.rarerich))

fungal.rarerich.wind.2013.brm <- general_brm(scale(fungal.rarerich)~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=fungal.df, family=gaussian(link="identity"))
summary(fungal.rarerich.wind.2013.brm)

y_w.fungal.rarerich.2013 <- as.numeric(scale(fungal.df$fungal.rarerich))
yrep_w.fungal.rarerich.2013 <- posterior_predict(fungal.rarerich.wind.2013.brm, nsamples=100)
#launch_shinystan(fungal.rarerich.wind.2013.brm)

## WIND BACTERIAL RAREFIED-RICHNESS 2013 ANALYSIS ----
hist(scale(bacteria.df$bacteria.rarerich))
bacteria.rarerich.wind.2013.brm <- general_brm(scale(bacteria.rarerich)~Wind.Exposure+(1+Wind.Exposure|Genotype)+(1|Block)+(1|Plot_code), data=bacteria.df, family=gaussian(link="identity"))
summary(bacteria.rarerich.wind.2013.brm)

y_w.bacteria.rarerich.2013 <- scale(bacteria.df$bacteria.rarerich)
yrep_w.bacteria.rarerich.2013 <- posterior_predict(bacteria.rarerich.wind.2013.brm, nsamples=100)
#launch_shinystan(bacteria.rarerich.wind.2013.brm)

## WIND ARTHROPOD COMPOSITION 2013 ANALYSIS ----
#w.arth.2013.comp <- select(w.arth.2013, X:spider_Larionoides, Plot_code) %>% # , GxE
#  gather(key=Species, value=Abundance, ant_F_obscuripes:spider_Larionoides) %>%
#  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

#composition_plot(w.arth.2013.comp, term="Wind.Exposure")
#composition_plot(w.arth.2013.comp, term="Genotype")
#composition_plot(w.arth.2013.comp, term="GxE")

## RETHINK WHETHER LONG DATA FORMAT IS NEEDED. SHOULDN'T NEED IT.
#hist(w.arth.2013.comp$Occurrence)
#arth.comp.wind.2013.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=w.arth.2013.comp, family=bernoulli(link="logit"))
#summary(arth.comp.wind.2013.brm)

#y_w.arth.comp.2013 <- w.arth.2013.comp$Occurrence
#yrep_w.arth.comp.2013 <- posterior_predict(arth.comp.wind.2013.brm, nsamples=100)
#launch_shinystan(arth.comp.wind.2013.brm)

## WIND SOIL PC1 ANALYSIS ----
hist(w.soil$soil.PC1)
hist(log(w.soil$soil.PC1-min(w.soil$soil.PC1)+1))
w.soil$soil.PC1.trans <- w.soil$soil.PC1-min(w.soil$soil.PC1)+1
soil.PC1.wind.brm <- general_brm(log(soil.PC1.trans)~Wind.Exposure+(1|Block), data=w.soil, family=gaussian(link="identity"))
summary(soil.PC1.wind.brm) 

y_w.soil.PC1 <- log(w.soil$soil.PC1.trans)
yrep_w.soil.PC1 <- posterior_predict(soil.PC1.wind.brm, nsamples=100)
#launch_shinystan(soil.PC1.wind.brm)

## WIND SOIL PC2 ANALYSIS ----
hist(w.soil$soil.PC2)
soil.PC2.wind.brm <- general_brm(soil.PC2~Wind.Exposure+(1|Block), data=w.soil, family=gaussian(link="identity"))
summary(soil.PC2.wind.brm)

y_w.soil.PC2 <- w.soil$soil.PC2
yrep_w.soil.PC2 <- posterior_predict(soil.PC2.wind.brm, nsamples=100)
#launch_shinystan(soil.PC2.wind.brm)

## ANT-APHID TRAIT-ARTHROPOD 2012 ANALYSIS ----
aa.trait.arth.2012 <- left_join(aa.arth.df, select(aa.trait.2012, plant_ID, trait.PC1, trait.PC2), by="plant_ID") 

trait.rich.aa.2012.brm <- general_brm(total.rich~scale(trait.PC1)+scale(trait.PC2)+(1|Block)+(1|Plot_code), data=aa.trait.arth.2012, family=poisson(link="log"))
summary(trait.rich.aa.2012.brm) # trait.PC1 is key driver, but there is also a weak effect of trait.PC2

y_aa.trait.rich.2012 <- aa.trait.arth.2012$total.rich
yrep_aa.trait.rich.2012 <- posterior_predict(trait.rich.aa.2012.brm, nsamples=100)
#launch_shinystan(trait.rich.aa.2012.brm)

plot(marginal_effects(trait.rich.aa.2012.brm, effects = "trait.PC1"), points=T)
plot(marginal_effects(trait.rich.aa.2012.brm, effects = "trait.PC2"), points=T)

## WIND TRAIT-ARTHROPOD 2012 ANALYSIS ----
w.trait.arth.2012 <- left_join(w.arth.2012, select(w.trait.2012, plant_ID, trait.PC1, trait.PC2), by="plant_ID") 

trait.rich.wind.2012.brm <- general_brm(total.rich~scale(trait.PC1)+scale(trait.PC2)+(1|Block)+(1|Plot_code), data=w.trait.arth.2012, family=poisson(link="log"))
summary(trait.rich.wind.2012.brm) # trait.PC1 is primary effect

y_w.trait.rich.2012 <- w.trait.arth.2012$total.rich
yrep_w.trait.rich.2012 <- posterior_predict(trait.rich.wind.2012.brm, nsamples=100)
#launch_shinystan(trait.rich.wind.2012.brm)

plot(marginal_effects(trait.rich.wind.2012.brm, effects = "trait.PC1"), points=T)


## WIND TRAIT-ARTHROPOD 2013 ANALYSIS ----
w.trait.arth.2013 <- left_join(w.arth.2013, select(w.trait.2013, plant_ID, trait.PC1, trait.PC2), by="plant_ID") 

trait.rich.wind.2013.brm <- general_brm(total.rich~scale(trait.PC1)+scale(trait.PC2)+(1|Block)+(1|Plot_code), data=w.trait.arth.2013, family=poisson(link="log"))
summary(trait.rich.wind.2013.brm) # trait.PC1 is the primary effect

plot(marginal_effects(trait.rich.wind.2013.brm, effects = "trait.PC1"), points=T)

## WIND TRAIT/SOIL-FUNGAL 2013 ANALYSIS ----
w.trait.fung.2013 <- left_join(fungal.df, select(w.trait.2013, plant_ID, trait.PC1, trait.PC2, root_CN), by="plant_ID") %>%
  left_join(., select(w.soil, soil.PC1, soil.PC2, Plot_code)) %>%
  mutate(log_root_CN = log(root_CN))

trait.rarerich.wind.2013.brm <- general_brm(scale(fungal.rarerich)~scale(trait.PC1)+scale(trait.PC2)+scale(soil.PC1)+scale(soil.PC2)+scale(log_root_CN)+(1|Block)+(1|Plot_code), data=w.trait.fung.2013, family=gaussian(link="identity"))
summary(trait.rarerich.wind.2013.brm) # log_root_CN and soil.PC1 to a lesser extent

plot(marginal_effects(trait.rarerich.wind.2013.brm, effects = "log_root_CN"), points=T)
plot(marginal_effects(trait.rarerich.wind.2013.brm, effects = "soil.PC1"), points=T)

## WIND TRAIT/SOIL-BACTERIA 2013 ANALYSIS ----
w.trait.bact.2013 <- left_join(bacteria.df, select(w.trait.2013, plant_ID, trait.PC1, trait.PC2, root_CN), by="plant_ID") %>%
  left_join(., select(w.soil, soil.PC1, soil.PC2, Plot_code)) %>%
  mutate(log_root_CN = log(root_CN))

bact.trait.rarerich.wind.2013.brm <- general_brm(scale(bacteria.rarerich)~scale(trait.PC1)+scale(trait.PC2)+scale(soil.PC1)+scale(soil.PC2)+scale(log_root_CN)+(1|Block)+(1|Plot_code), data=w.trait.bact.2013, family=gaussian(link="identity"))
summary(bact.trait.rarerich.wind.2013.brm) # strong effect of soil.PC2

plot(marginal_effects(bact.trait.rarerich.wind.2013.brm, effects = "soil.PC2"), points=T)

## TIDY AND SAVE OUTPUT ----

lanphere_SDs <- bind_rows(
  mutate(wind_posterior_SDs(trait.PC1.wind.2012.brm, df=w.trait.2012), Experiment="Wind", Year="2012", Response="Trait PC1"),
  mutate(wind_posterior_SDs(trait.PC2.wind.2012.brm, df=w.trait.2012), Experiment="Wind", Year="2012", Response="log(Trait PC2 + min(Trait PC2) + 1)"),
  mutate(wind_posterior_SDs(trait.PC1.wind.2013.brm, df=w.trait.2013), Experiment="Wind", Year="2013", Response="Trait PC1"),
  mutate(wind_posterior_SDs(trait.PC2.wind.2013.brm, df=w.trait.2013), Experiment="Wind", Year="2013", Response="Trait PC2"),
  mutate(wind_posterior_SDs(soil.PC1.wind.brm, df=w.soil), Experiment="Wind", Year="2013", Response="log(Soil PC1 + min(Soil PC1) + 1)"),
  mutate(wind_posterior_SDs(soil.PC2.wind.brm, df=w.soil), Experiment="Wind", Year="2013", Response="Soil PC2"),
  mutate(wind_posterior_SDs(root_CN.wind.2013.brm, df=w.trait.2013), Experiment="Wind", Year="2013", Response="log(Root C:N)"),
  mutate(ant.aphid_posterior_SDs(trait.PC1.aa.2012.brm, df=aa.trait.2012), Experiment="Ant-Aphid", Year="2012", Response="Trait PC1"),
  mutate(ant.aphid_posterior_SDs(trait.PC2.aa.2012.brm, df=aa.trait.2012), Experiment="Ant-Aphid", Year="2012", Response="Trait PC2"),
  mutate(wind_posterior_SDs(arth.rich.wind.2012.brm, df=w.arth.2012), Experiment="Wind", Year="2012", Response="Arthropod Richness"),
  #mutate(wind_posterior_SDs(arth.comp.wind.2012.brm), Experiment="Wind", Year="2012", Response="Arthropod Composition"),
  mutate(wind_posterior_SDs(arth.rich.wind.2013.brm, df=w.arth.2013), Experiment="Wind", Year="2013", Response="Arthropod Richness"),
  #mutate(wind_posterior_SDs(arth.comp.wind.2013.brm), Experiment="Wind", Year="2013", Response="Arthropod Composition"),
  mutate(wind_posterior_SDs(fungal.rarerich.wind.2013.brm, df=fungal.df), Experiment="Wind", Year="2013", Response="scale(Fungi Rarefied Richness)"),
  mutate(wind_posterior_SDs(bacteria.rarerich.wind.2013.brm, df=bacteria.df), Experiment="Wind", Year="2013", Response="scale(Bacteria Rarefied Richness)"),
  mutate(ant.aphid_posterior_SDs(arth.rich.aa.2012.brm, df=aa.arth.df), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness")#,
  #mutate(ant.aphid_posterior_SDs(arth.comp.aa.2012.brm), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Composition")
)
write_csv(lanphere_SDs, path="output_brms/lanphere_SDs.csv")

## CONSIDER CHANGING BACTERIAL TO BACTERIA AND FUNGAL TO FUNGI TO SHORTEN EVERYTHING

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

lanphere_trait_regs <- bind_rows(
  mutate(posterior_samples(trait.rich.aa.2012.brm, pars = "^b"), Experiment="Ant-Aphid", Year="2012", Response="Arthropod Richness"),
  mutate(posterior_samples(trait.rich.wind.2012.brm, pars = "^b"), Experiment="Wind", Year="2012", Response="Arthropod Richness"),
  mutate(posterior_samples(trait.rich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="Arthropod Richness"),
  mutate(posterior_samples(trait.rarerich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="scale(Fungi Rarefied Richness)"),
  mutate(posterior_samples(bact.trait.rarerich.wind.2013.brm, pars = "^b"), Experiment="Wind", Year="2013", Response="scale(Bacteria Rarefied Richness)")
)
write_csv(lanphere_trait_regs, path="output_brms/lanphere_trait_regs.csv")

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
