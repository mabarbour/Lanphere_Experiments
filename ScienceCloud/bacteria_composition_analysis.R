# args object will contain all the arguments given from the command line e.g. Rscript example.R xyz.csv
#args <- commandArgs(trailingOnly = TRUE)

# filname will contain the csv file
#filename <- args[1]

## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(broom)
library(parallel)

#options(mc.cores=parallel::detectCores ())

## LOAD BACTERIA COMMUNITY DATA ----

bacteria.df <- read.csv(file = "final_data/bacteria.df.csv", header=TRUE) %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))

bacteria.comp <- select(bacteria.df, X:plant_ID, OTU_1347:OTU_713) %>%
  gather(key=Species, value=Abundance, OTU_1347:OTU_713) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0)) #%>%
  #sample_n(size=100)

## FUNCTIONS FOR ANALYSIS ----

general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=family, 
      prior=prior(normal(0,1), class=sd),
      control=list(adapt_delta=0.99),
      chains=1,
      cores=10)
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

bacteria.comp.wind.2013.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=bacteria.comp, family=bernoulli(link="logit"))


## TIDY AND SAVE OUTPUT ----

lanphere_bacteria_VarComps <- mutate(
  get_VarComps(bacteria.comp.wind.2013.brm, Distrib_Var=pi^2/3), 
               Experiment="Wind", Year="2013", Trait="Bacterial Composition"
)
write_csv(lanphere_bacteria_VarComps, path="results.csv")
