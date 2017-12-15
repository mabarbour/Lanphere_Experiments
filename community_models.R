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
#contrasts(wind.arth.df$Wind.Exposure) <- "contr.sum" # important for calculating variance explained
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
species_guild_info$Key_Guilds <- relevel(as.factor(c(rep("low",2), rep("Aphididae",3),
                                   rep("low",16), "Tortricidae_leaftier",
                                   "Gracilliaridae_leafminer", rep("low",3),
                                   rep("Spider", 7), "low", "Gracilliaridae_leafminer")), ref="low")
with(species_guild_info, table(Guild, Key_Guilds))

wind.arth.2012.df <- left_join(wind.arth.2012, species_guild_info)
#contrasts(wind.arth.2012.df$Wind.Exposure) <- "contr.sum"
#relevel(as.factor(wind.arth.2012.df$Key_Guilds), ref="low")
sum(wind.arth.2012.df$Occurrence)/length(wind.arth.2012.df$Occurrence) # only 2% are ones...

## CONSIDER MODELING GENOTYPE AND SPECIES AS FIXED-EFFECTS. THIS MAKES SENSE IF THERE
## EFFECTS ARE NOT NORMALLY DISTRIBUTED, WHICH I CAN BASICALLY TEST FOR BY LOOKING AT
## THE DISTRIBUTION OF COEFFICIENTS WHEN THEY ARE TREATED AS FIXED-EFFECTS.

wind.arth.model.base <- brm(Occurrence ~ Wind.Exposure + (Wind.Exposure|Genotype) + (1|Species) +
                         (1|plant_ID) + (1|Block) + (1|Block:Wind.Exposure),
                       data=wind.arth.2012.df, family=bernoulli(link="logit"),
                       algorithm="sampling",
                       prior=c(prior(normal(0,2), class=b),
                               prior(normal(0,2), class=sd)),
                       chains=4,
                       control=list(adapt_delta=0.99))

wind.arth.model <- brm(Occurrence ~ 
                         Wind.Exposure + (Wind.Exposure|Genotype*Species) +
                         (1|plant_ID) + (1|Block) + (1|Block:Wind.Exposure),
                       data=wind.arth.2012.df,
                       family=bernoulli(link="logit"),
                       algorithm="sampling",
                       prior=c(prior(normal(0,2), class=b),
                               prior(normal(0,2), class=sd)),
                       chains=4,
                       control=list(adapt_delta=0.99))
tidy(wind.arth.model)
shinystan::launch_shinystan(wind.arth.model)

wind.arth.2012.df$Species <- as.factor(wind.arth.2012.df$Species)
contrasts(wind.arth.2012.df$Species) <- "contr.sum"

## VERY POOR MODEL, STICK WITH RANDOM EFFECTS AND JUST ESTIMATE THE PARAMETERS.
wind.arth.model.FE <- brm(Occurrence ~ 
                         (Wind.Exposure + Species + Genotype)^2 +
                         (1|plant_ID) + (1|Block) + (1|Block:Wind.Exposure),
                       data=wind.arth.2012.df,
                       family=bernoulli(link="logit"),
                       algorithm="sampling",
                       prior=c(prior(normal(0,2), class=b),
                               prior(normal(0,2), class=sd)),
                       chains=1,
                       control=list(adapt_delta=0.99))
summary(wind.arth.model.FE)
tidy(wind.arth.model.FE)
FE.df <- data.frame(fixef(wind.arth.model.FE))
hist(FE.df[3:36,"Estimate"])
#hist(FE.df[12:45,"Estimate"])

LOO(wind.arth.model.base, wind.arth.model, wind.arth.model.FE)

#plot(wind.arth.model)
summary(wind.arth.model)
tidy(wind.arth.model)[1:15,]
ranef(wind.arth.model)
pps <- data.frame(t(posterior_predict(wind.arth.model, nsamples=100)))

## EXAMINE RICHNESS PREDICTIONS
wind.arth.2012.df_pps <- cbind.data.frame(wind.arth.2012.df, pps)

rich.df <- wind.arth.2012.df_pps %>%
  group_by(Block, Wind.Exposure, Genotype, plant_ID) %>%
  summarise_at(vars(Occurrence, X1:X100), sum) %>%
  ungroup()
rich.df$predict.rich <- rowMeans(select(rich.df, X1:X100))
rich.df$error <- rich.df$Occurrence - rich.df$predict.rich
rich.df$MSE <- rich.df$error^2

rich.df.gather <- select(rich.df, -Occurrence, -predict.rich, -error, -MSE) %>%
  gather(key=sim, value=Occurrence.predict, X1:X100)

ggplot(rich.df.gather, aes(x=Genotype, y=Occurrence.predict)) + geom_boxplot()

# model appears to over predict richness when it is zero, and under predict it when it is higher
#ggplot(rich.df, aes(x = Abundance, y=predict.rich)) + geom_point() + geom_abline(intercept=0, slope=1) + geom_smooth(method="lm") + scale_y_continuous(limits=c(0,50))
#ggplot(rich.df, aes(x = Abundance, y=error)) + geom_point() + geom_hline(yintercept=0)

rich.df.genos <- rich.df %>%
  group_by(Genotype) %>%
  summarise_at(vars(Occurrence, predict.rich, X1:X100), mean) %>%
  ungroup() %>%
  gather(key=sim, value=Occurrence.predict, X1:X100)
ggplot(rich.df.genos, aes(x=Genotype, y=Occurrence.predict)) + 
  stat_sum_df("mean_cl_boot") +
  #geom_point(shape=1, color="grey") +
  geom_point(aes(y=Occurrence))

rich.df.wind <- rich.df %>%
  group_by(Wind.Exposure) %>%
  summarise_at(vars(Occurrence, predict.rich, X1:X100), mean) %>%
  ungroup() %>%
  gather(key=sim, value=Occurrence.predict, X1:X100)
ggplot(rich.df.wind, aes(x=Wind.Exposure, y=Occurrence.predict)) + 
  stat_sum_df("mean_cl_boot") +
  #geom_boxplot() + 
  geom_point(aes(y=Occurrence))

rich.df.block <- rich.df %>%
  group_by(Block) %>%
  summarise_at(vars(Occurrence, predict.rich, X1:X100), mean) %>%
  ungroup() %>%
  gather(key=sim, value=Occurrence.predict, X1:X100)
ggplot(rich.df.block, aes(x=Block, y=Occurrence.predict)) + 
  stat_sum_df("mean_cl_boot") +
  #geom_boxplot() + 
  geom_point(aes(y=Occurrence))

# A set of useful summary functions is provided from the Hmisc package:
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, colour = "red", geom = geom, width = 0.2, ...)
}

table(filter(wind.arth.2012.df, Occurrence>0)$Guild)
table(filter(wind.arth.2012.df, Occurrence>0)$Guild)/sum(table(filter(wind.arth.2012.df, Occurrence>0)$Guild))
which(table(filter(wind.arth.2012.df, Occurrence>0)$Guild)/sum(table(filter(wind.arth.2012.df, Occurrence>0)$Guild)) > 0.05)
28+16.5+28+12.9 # 85% of the arthropods
