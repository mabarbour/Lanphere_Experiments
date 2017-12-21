
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(rstan)
library(broom)
library(parallel)

## SET OPTIONS ----

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## LOAD BACTERIA COMMUNITY DATA ----

microbial.df <- read.csv('bacteria.df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))

# comment out appropriate OTUs
#OTUs <- colnames(select(microbial.df, -(X:fungal.rarerich))) # Fungal Composition
OTUs <- colnames(select(microbial.df, -(X:bacteria.rarerich))) # Bacteria Composition

microbial.comp <- microbial.df[ ,c("Block", "Plot_code", "Wind.Exposure", "Genotype", OTUs)] %>%
  gather_(key="Species", value="Abundance", OTUs) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

## FUNCTIONS FOR ANALYSIS ----

general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=family, 
      prior=prior(normal(0,1), class=sd),
      control=list(adapt_delta=0.99),
      chains=1)
  # all other brm parameters correspond to the defaults
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


## WIND BACTERIA COMPOSITION 2013 ANALYSIS ----

microbial.comp.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=microbial.comp, family=bernoulli(link="logit"))
summary(microbial.comp.brm)

## TIDY AND SAVE OUTPUT ----

microbial_VarComps <- mutate(
  get_VarComps(microbial.comp.brm, Distrib_Var=pi^2/3), 
               Experiment="Wind", Year="2013", Trait="Bacteria Composition" # "Bacteria or Fungal Composition"
)
write_csv(microbial_VarComps, path="bacteria_composition_VarComps.csv")
