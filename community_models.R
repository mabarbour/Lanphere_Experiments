############################################
## Description: This script analyzes arthropod community data for the Lanphere experiments.
## Code author: Matt Barbour
## Email: matthew.a.barbour@gmail.com
############################################

## load required libraries ----
source('scripts_for_analysis/required_libraries.R')

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
contrasts(wind.arth.df$Wind.Exposure) <- "contr.sum" # important for calculating variance explained
wind.arth.names <- colnames(select(wind.arth.df, Gracilliaridae_miner:Spider)) # for subsetting community data


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

## COMMUNITY MODELS ----
library(brms)
library(bayesplot)
library(broom)
source('scripts_for_analysis/variance_partitioning_brms.R')

wind.arth.2012 <- filter(wind.arth.df, Year=="2012") %>%
  select(X:spider_Larionoides) %>%
  gather(key=Species, value=Abundance, ant_F_obscuripes:spider_Larionoides) %>%
  mutate(Occurrence=ifelse(Abundance>0, 1, 0))

species_guild_info <- data.frame(Species = names(table(wind.arth.2012$Species)),
  Guild = c(rep("Formica_ant",2), rep("Aphididae",3), rep("Caterpillar_other",3),
            "Erophyidae_gall", "Cecidomyiidae_gall", "Tenthredinidae_sawfly", 
            rep("Cecidomyiidae_gall",3), "Orthoptera", rep("Cicadellidae",6), 
            "Tortricidae_leaftier", "Gracilliaridae_leafminer", "Psyllidae", 
            "Coccoidea", "Tenthredinidae_sawfly", rep("Spider",7), 
            "Pentatomidae", "Gracilliaridae_leafminer"))
table(species_guild_info) # everything corresponds  
species_guild_info$Key_Guilds <- c(rep("low",2), rep("Aphididae",3),
                                   rep("low",16), "Tortricidae_leaftier",
                                   "Gracilliaridae_leafminer", rep("low",3),
                                   rep("Spider", 7), "low", "Gracilliaridae_leafminer")
with(species_guild_info, table(Guild, Key_Guilds))

wind.arth.2012.df <- left_join(wind.arth.2012, species_guild_info)
contrasts(wind.arth.2012.df$Wind.Exposure) <- "contr.sum"
relevel(as.factor(wind.arth.2012.df$Key_Guilds), ref="low")

wind.arth.model <- brm(Occurrence ~ 
                         Wind.Exposure*Key_Guilds + (Key_Guilds + Wind.Exposure|Genotype) + (1|Species) +
                         (1|plant_ID) + (1|Block) + (1|Block:Wind.Exposure),
                       data=wind.arth.2012.df,
                       family=bernoulli(link="logit"),
                       algorithm="sampling",
                       chains=1,
                       control=list(adapt_delta=0.9))
plot(wind.arth.model)
summary(wind.arth.model)
plot(marginal_effects(wind.arth.model))
ppc_bars_grouped(y=wind.arth.2012.df$Occurrence, 
                 yrep=posterior_predict(wind.arth.model, nsamples=100),
                 group=wind.arth.2012.df$Guild)

binomial_variance <- pi^2/3
var.comp <- get_variance_wind_brms(wind.arth.model, data=wind.arth.2012.df, RE_rows=3:9)
var.comp$perc_var <- var.comp$estimate^2/(sum(var.comp$estimate^2)+binomial_variance)

pps <- data.frame(t(posterior_predict(wind.arth.model, nsamples=100)))

## EXAMINE RICHNESS PREDICTIONS
wind.arth.2012.df_pps <- cbind.data.frame(wind.arth.2012.df, pps)

rich.df <- wind.arth.2012.df_pps %>%
  group_by(Block, Wind.Exposure, Genotype, plant_ID) %>%
  summarise_at(vars(Occurrence, X1:X100), sum) %>%
  ungroup()

rich.df$predict.rich <- rowMeans(select(rich.df, X1:X100))
rich.df$MSE <- (rich.df$Occurrence - rich.df$predict.rich)^2
ggplot(rich.df, aes(x = Occurrence, y=predict.rich)) + geom_point() + geom_abline(intercept=0, slope=1)

ggplot(rich.df, aes(x = Occurrence, y=MSE)) + geom_point() 
## seems that the model overestimates richness when it is observed to be zero,
## and underestimates richness when it is high.

table(filter(wind.arth.2012.df, Occurrence>0)$Guild)
table(filter(wind.arth.2012.df, Occurrence>0)$Guild)/sum(table(filter(wind.arth.2012.df, Occurrence>0)$Guild))
which(table(filter(wind.arth.2012.df, Occurrence>0)$Guild)/sum(table(filter(wind.arth.2012.df, Occurrence>0)$Guild)) > 0.05)
28+16.5+28+12.9 # 85% of the arthropods
