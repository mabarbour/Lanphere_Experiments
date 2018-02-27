
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(rstan)
library(broom)
library(parallel)

## SET OPTIONS ----

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## LOAD MICROBIAL COMMUNITY DATA ----

fungal.df <- read.csv(file = "fungal.df.csv") %>% # change to appropriate directory if running in Rstudio
  tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
bacteria.df <- read.csv(file = "bacteria.df.csv") %>% # change to appropriate directory if running in Rstudio
  tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))

# get OTUs
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich))) # Fungal Composition
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich))) # Bacteria Composition

# divide up fungi
fungi.occurrs <- colSums(fungal.df[ ,f.OTUs]>0)
summary(fungi.occurrs) # Fungi: 100-75% = 128 to 10; 75-50% = 10 to 3; 50-25% = 3 to 1; 25-0% = 1 to 0
abund.Fungi <- f.OTUs[which(fungi.occurrs>10)]
mid.Fungi <- f.OTUs[which(fungi.occurrs<11 & fungi.occurrs>2)]

fungal.comp <- fungal.df[ ,c("Block", "Plot_code", "Wind.Exposure", "Genotype", f.OTUs)] %>%
  gather_(key="Species", value="Abundance", f.OTUs) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0)) 
fungal.comp.abund <- filter(fungal.comp, Species %in% abund.Fungi)
fungal.comp.mid <- filter(fungal.comp, Species %in% mid.Fungi)

# divide up bactieria
bacteria.occurrs <- colSums(bacteria.df[ ,b.OTUs]>0)
summary(bacteria.occurrs)
abund.Bacteria <- b.OTUs[which(bacteria.occurrs>28)]
mid.Bacteria <- b.OTUs[which(bacteria.occurrs<29 & bacteria.occurrs>8)]
low.Bacteria <- b.OTUs[which(bacteria.occurrs<9 & bacteria.occurrs>3)]

bacteria.comp <- bacteria.df[ ,c("Block", "Plot_code", "Wind.Exposure", "Genotype", b.OTUs)] %>%
  gather_(key="Species", value="Abundance", b.OTUs) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0)) 
bacteria.comp.abund <- filter(bacteria.comp, Species %in% abund.Bacteria)
bacteria.comp.mid <- filter(bacteria.comp, Species %in% mid.Bacteria)
bacteria.comp.low <- filter(bacteria.comp, Species %in% low.Bacteria)


## FUNCTIONS FOR ANALYSIS ----

general_brm <- function(formula, family, data, ...) {
  brm(formula=formula, data=data, family=family, 
      prior=prior(normal(0,1), class=sd),
      #algorithm="sampling",
      control=list(adapt_delta=0.999, max_treedepth=20),
      chains=4)
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


## WIND FUNGAL COMPOSITION 2013 ANALYSIS ----

fungal.comp.abund.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=fungal.comp.abund, family=bernoulli(link="logit"))
summary(fungal.comp.abund.brm)

fungal.comp.mid.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=fungal.comp.mid, family=bernoulli(link="logit"))
summary(fungal.comp.mid.brm)

## WIND BACTERIAL COMPOSITION 2013 ANALYSIS ----

bacteria.comp.abund.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=bacteria.comp.abund, family=bernoulli(link="logit"))
summary(bacteria.comp.abund.brm)

bacteria.comp.mid.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=bacteria.comp.mid, family=bernoulli(link="logit"))
summary(bacteria.comp.mid.brm)

bacteria.comp.low.brm <- general_brm(Occurrence~(1|Genotype*Wind.Exposure*Species)+(1|Block)+(1|Plot_code)+(1|Block:Species), data=bacteria.comp.low, family=bernoulli(link="logit"))
summary(bacteria.comp.low.brm)

## TIDY AND SAVE OUTPUT ----

microbial_VarComps <- bind_rows(
  mutate(get_VarComps(fungal.comp.abund.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Abundant Fungi Composition"),
  mutate(get_VarComps(fungal.comp.mid.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Mid Fungi Composition"),
  mutate(get_VarComps(bacteria.comp.abund.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Abundant Bacteria Composition"),
  mutate(get_VarComps(bacteria.comp.mid.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Mid Bacteria Composition"),
  mutate(get_VarComps(bacteria.comp.low.brm, Distrib_Var=pi^2/3), 
         Experiment="Wind", Year="2013", Trait="Low Bacteria Composition")
)
write_csv(microbial_VarComps, path="microbial_composition_VarComps.csv")
