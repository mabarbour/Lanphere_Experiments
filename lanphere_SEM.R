## load required libraries ----
#library(devtools)
#install_github("gavinsimpson/ggvegan")
library(ggvegan) # requires devtools and install_github()
library(cowplot)
#source('~/Documents/miscellaneous_R/autoplot.custom.R')
#library(devtools)
#install_github("jslefche/piecewiseSEM")
library(piecewiseSEM) # dev version as of Nov. 10
#library(missMDA) # for imputing missing values in PCA
library(dplyr)
library(tidyr) # for separate()
library(psych) # for correlations
library(lme4)
library(lmerTest)
library(vegan) # community composition analysis
library(effects)
library(car)
library(RCurl) # get code directly from github

## source and parse github code
script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/autoplot.custom.R", ssl.verifypeer = FALSE)
eval(parse(text = script))


## Wind SEM data management ----
# combine required datasets. First load datasets for lanphere_plant_trait_analysis.R and lanphere_arthropod_analysis.R and lanphere_root_analysis.R

## load and analyze soil properties at plot-level ----
w.soil <- read.csv('final_data/wind_soil_df.csv') %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block),
         Plot_code = paste(Block, Wind.Exposure, sep = "."),
         num.Wind = ifelse(Wind.Exposure == "Exposed", 2, 1))

# Soil PC1
soil.PC1.lmer <- lmerTest::lmer(scale(soil.PC1) ~ scale(num.Wind) + (1|Block), w.soil)
summary(soil.PC1.lmer) # effect is not significant
fixef(soil.PC1.lmer) # SD = -0.26, wind exposure results in reduced soil moisture, organic matter and microelements, but more available N03 and NH4.
sem.model.fits(soil.PC1.lmer) # marginal R2 = 0.06

# Soil PC2
soil.PC2.lmer <- lmerTest::lmer(scale(soil.PC2) ~ scale(num.Wind) + (1|Block), w.soil)
summary(soil.PC2.lmer) # effect is not significant
fixef(soil.PC2.lmer) # SD = 0.20, wind exposure does influence soil chemistry
sem.model.fits(soil.PC2.lmer) # marginal R2 = 0.04

## load fungal and bacteria community data ----
fungal.df <- read.csv("fungal.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code = paste(Block, Wind.Exposure, sep = "."))
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich)))

bacteria.df <- read.csv("bacteria.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code = paste(Block, Wind.Exposure, sep = "."))
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich)))

## load arthropod community data ----
wind.arth.df <- read.csv('~/Lanphere_Experiments/final_data/wind_arthropod_df.csv') %>%
  tbl_df() %>% #rename(plant_ID = plant_code) %>%
  mutate(Block = as.factor(Block), Year = as.factor(Year), Plot_code = paste(Block, Wind.Exposure, sep = "."))

w.arth.names <- colnames(select(wind.arth.df, Gracilliaridae_miner:Spider))

## load plant trait data ----
w.trait.df <- read.csv('~/Lanphere_Experiments/final_data/wind_trait_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Year = as.factor(Year), Plot_code = paste(Block, Wind.Exposure, sep = "."))
glimpse(w.trait.df)

## combine datasets into data frame ----
tmp.combo <- left_join(select(wind.arth.df, -X), select(w.trait.df, -X)) %>% left_join(., select(w.soil, Plot_code, soil.PC1, soil.PC2)) 
w.SEM.2012 <- tmp.combo %>% filter(Year == "2012") %>% select(Block, Wind.Exposure, Genotype, plant_ID, total.abund:total.rarerich, trait.PC1, trait.PC2, soil.PC1, soil.PC2) %>% mutate(num.Wind = ifelse(Wind.Exposure == "Exposed", 2, 1)) 

w.SEM.2013 <- tmp.combo %>% 
  filter(Year == "2013") %>% 
  select(Block, Wind.Exposure, Genotype, plant_ID, total.abund:total.rarerich, trait.PC1, trait.PC2, soil.PC1, soil.PC2, root_CN) %>% 
  left_join(., select(fungal.df, plant_ID, fungal.rarerich)) %>% # fungal.abund:
  left_join(., select(bacteria.df, plant_ID, bacteria.rarerich)) %>% # bacteria.abund:
  mutate(num.Wind = ifelse(Wind.Exposure == "Exposed", 2, 1)) 

## look at correlation among univariate communtity responses in 2013 ----
corr.test(dplyr::select(w.SEM.2013, total.abund:total.rarerich, bacteria.rarerich, fungal.rarerich)) # bacteria.abund: fungal.abund:

## Scale variables to mean = 0 and SD = 1 so coefficients are standardized ----
vars.2013 <- colnames(select(w.SEM.2013, num.Wind, total.rich, total.abund, total.rarerich, fungal.rarerich, bacteria.rarerich, trait.PC1, trait.PC2, root_CN, soil.PC1, soil.PC2)) # fungal.rich, fungal.abund, bacteria.rich, bacteria.abund, 
w.SEM.2013[ ,vars.2013] <- scale(w.SEM.2013[ ,vars.2013])

## Wind: richness models ----
# fit separate models
arth.rich.2013 <- lmerTest::lmer(total.rich ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013) # removed Genotype to avoid modelling direct effect
fixef(arth.rich.2013)
sem.model.fits(arth.rich.2013)

#fung.rich.2013 <- lmerTest::lmer(fungal.rich ~ soil.PC1 + soil.PC2 + root_CN + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013)
#sem.model.fits(fung.rich.2013)

#bact.rich.2013 <- lmerTest::lmer(bacteria.rich ~ soil.PC1 + soil.PC2 + root_CN + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013)
#sem.model.fits(bact.rich.2013)

soilPC1.2013 <- lmerTest::lmer(soil.PC1 ~ num.Wind + (1|Block), data = w.SEM.2013)

soilPC2.2013 <- lmerTest::lmer(soil.PC2 ~ num.Wind + (1|Block), data = w.SEM.2013) # note that the coefficient is approximately correct, however, the p-value for is not based on the appropriate residual degrees of freedom (way too liberal). This means I will have to manually adjust the p-value for this path-coefficient in the structural equation model.

rootCN.2013 <- lmerTest::lmer(root_CN ~ soil.PC1 + soil.PC2 + Genotype + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013, contrasts = list(Genotype = "contr.sum"))
sd(fixef(rootCN.2013)[-(1:3)]) # Genotype SD = 0.23
sd(fixef(rootCN.2013)[-(1:3)])*25 # plotting scale = 5.6
sem.model.fits(rootCN.2013) # marginal R2 = 0.05

traitPC2.2013 <- lmerTest::lmer(trait.PC2 ~ soil.PC1 + soil.PC2 + num.Wind + Genotype + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013, contrasts = list(Genotype = "contr.sum"))
sd(fixef(traitPC2.2013)[-(1:4)]) # Genotype SD = 0.37,
sd(fixef(traitPC2.2013)[-(1:4)])*25 # plotting scale = 9.3
sem.model.fits(traitPC2.2013) # marginal R2 = 0.15

traitPC1.2013 <- lmerTest::lmer(trait.PC1 ~ soil.PC1 + soil.PC2 + num.Wind + Genotype + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013, contrasts = list(Genotype = "contr.sum"))
sd(fixef(traitPC1.2013)[-(1:4)]) # Genotype SD = 0.38,
sd(fixef(traitPC1.2013)[-(1:4)])*25 # plotting scale = 9.4
sem.model.fits(traitPC1.2013) # marginal R2 = 0.19

# generate model list
mods <- list(traitPC1.2013,traitPC2.2013, soilPC1.2013, soilPC2.2013, arth.rich.2013)#,rootCN.2013, fung.rich.2013, bact.rich.2013)

# fit SEM
sem.fit(mods, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2")) # adequate fit to the data: Fisher.C = 2.94, df = 10, P = 0.983; 

mods.coefs <- sem.coefs(mods, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), standardize = "none") #"fungal.rich~~bacteria.rich",
mods.coefs <- mutate(mods.coefs, plot.scale = estimate*25)
filter(mods.coefs, predictor %in% c("num.Wind","soil.PC1","soil.PC2","trait.PC1","trait.PC2","root_CN","~~ trait.PC2","~~ bacteria.rich", "~~ soil.PC2")) # note that I used original models for standardized coefficient estimates on num.Wind -> soil.PC1 and soil.PC2 because they have the appropriate degrees of freedom (unlike the ones in this model). # why are the correlations so high?

## Wind: abundance models 2013 ----
arth.abund.2013 <- lmerTest::lmer(total.abund ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013) # removed Genotype to avoid modelling direct effect
fixef(arth.abund.2013)
sem.model.fits(arth.abund.2013)

#fung.abund.2013 <- lmerTest::lmer(fungal.abund ~ soil.PC1 + soil.PC2 + root_CN + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013)
#sem.model.fits(fung.abund.2013)

#bact.abund.2013 <- lmerTest::lmer(bacteria.abund ~ soil.PC1 + soil.PC2 + root_CN + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013)
#sem.model.fits(bact.abund.2013)

mods.abund <- list(traitPC1.2013,traitPC2.2013, soilPC1.2013, soilPC2.2013, arth.abund.2013) #, rootCN.2013, fung.abund.2013, bact.abund.2013)

sem.fit(mods.abund, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2")) # adequate fit to the data: Fisher.C = 6.73, df = 10, P = 0.751

mods.abund.coefs <- sem.coefs(mods.abund, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), standardize = "none")
mods.abund.coefs <- mutate(mods.abund.coefs, plot.scale = estimate*25)
filter(mods.abund.coefs, predictor %in% c("num.Wind","soil.PC1","soil.PC2","trait.PC1","trait.PC2","root_CN"))

## Wind: rarefied richness models 2013 ----
arth.rarerich.2013 <- lmerTest::lmer(total.rarerich ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block),# + (1|Block:Wind.Exposure), 
                                     data = w.SEM.2013) 
sem.model.fits(arth.rarerich.2013)
fixef(arth.rarerich.2013)
-0.26*25

fung.rarerich.2013 <- lmerTest::lmer(fungal.rarerich ~ soil.PC1 + soil.PC2 + root_CN + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2013)
fixef(fung.rarerich.2013)
sem.model.fits(fung.rarerich.2013)

bact.rarerich.2013 <- lmerTest::lmer(bacteria.rarerich ~ soil.PC1 + soil.PC2 + root_CN + (1|Block),# + (1|Block:Wind.Exposure), 
                                     data = w.SEM.2013)
fixef(bact.rarerich.2013)
sem.model.fits(bact.rarerich.2013)

# for an unknown reason, I receive an error when I try to include all of the models together. I'm going to work around this by manually calculating Fisher's C, df, and P-value.
mods.rarerich <- list(traitPC1.2013, traitPC2.2013, soilPC1.2013, soilPC2.2013, rootCN.2013, fung.rarerich.2013, bact.rarerich.2013)#, arth.rarerich.2013)
mods.rarerich2 <- list(traitPC1.2013, traitPC2.2013, soilPC1.2013, soilPC2.2013, arth.rarerich.2013) 

sem.fit(mods.rarerich, data = w.SEM.2013, corr.errors = c("fungal.rarerich~~bacteria.rarerich", "trait.PC1~~trait.PC2","soil.PC1~~soil.PC2")) # adequate fit to the data: Fisher.C = 22.26, df = 38, P = 0.98; AIC = 174.26, AICc = 6026.26, K = 76, n = 79
sem.fit(mods.rarerich2, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"))

rarerich.sem <- sem.missing.paths(mods.rarerich, w.SEM.2013, corr.errors = c("fungal.rarerich~~bacteria.rarerich","trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), conditional = TRUE) %>% 
  rbind.data.frame(., sem.missing.paths(mods.rarerich2, data = w.SEM.2013, corr.errors = c("trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), conditional = TRUE)[-c(1,2), ]) # remove overlapping missing paths

rarerich.sem.C <- -2*sum(log(rarerich.sem$p.value)) # Fisher.C = 22.32
rarerich.df <- 2 * length(rarerich.sem$p.value)*2 # df = 32
1 - pchisq(rarerich.sem.C, rarerich.df) # P = 0.899

#mods.rarerich.coefs <- sem.coefs(mods.rarerich2, data = w.SEM.2013, corr.errors = c("fungal.rarerich~~bacteria.rarerich","trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), standardize = "none")
#mods.rarerich.coefs <- mutate(mods.rarerich.coefs, plot.scale = estimate*25)
#filter(mods.rarerich.coefs, predictor %in% c("num.Wind","soil.PC1","soil.PC2","trait.PC1","trait.PC2","root_CN"))

## Wind: arthropods 2012 ----
vars.2012 <- colnames(select(w.SEM.2012, num.Wind, total.rich, total.abund, total.rarerich, trait.PC1, trait.PC2,soil.PC1, soil.PC2))
w.SEM.2012[ ,vars.2012] <- scale(w.SEM.2012[ ,vars.2012])

# fit separate models
arth.rich.2012 <- lmerTest::lmer(total.rich ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2012) # removed Genotype to avoid modelling direct effect
summary(arth.rich.2012)
sem.model.fits(arth.rich.2012)

arth.abund.2012 <- lmerTest::lmer(total.abund ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2012) # removed Genotype to avoid modelling direct effect
summary(arth.abund.2012)
sem.model.fits(arth.abund.2012)

arth.rarerich.2012 <- lmerTest::lmer(total.rarerich ~ trait.PC1 + trait.PC2 + num.Wind + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2012) # removed Genotype to avoid modelling direct effect
summary(arth.rarerich.2012)
sem.model.fits(arth.rarerich.2012)

soilPC1.2012 <- lmerTest::lmer(soil.PC1 ~ num.Wind + (1|Block), data = w.SEM.2012)

soilPC2.2012 <- lmerTest::lmer(soil.PC2 ~ num.Wind + (1|Block), data = w.SEM.2012) # note that the coefficient is approximately correct, however, the p-value for is not based on the appropriate residual degrees of freedom (way too liberal). This means I will have to manually adjust the p-value for this path-coefficient in the structural equation model.

traitPC2.2012 <- lmerTest::lmer(trait.PC2 ~ soil.PC1 + soil.PC2 + num.Wind + Genotype + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2012, contrasts = list(Genotype = "contr.sum"))
sd(fixef(traitPC2.2012)[-(1:4)]) # Genotype SD = 0.31,
sd(fixef(traitPC2.2012)[-(1:4)])*25 # plotting scale = 7.6
sem.model.fits(traitPC2.2012) # marginal R2 = 0.20

traitPC1.2012 <- lmerTest::lmer(trait.PC1 ~ soil.PC1 + soil.PC2 + num.Wind + Genotype + (1|Block) + (1|Block:Wind.Exposure), data = w.SEM.2012, contrasts = list(Genotype = "contr.sum"))
sd(fixef(traitPC1.2012)[-(1:4)]) # Genotype SD = 0.53,
sd(fixef(traitPC1.2012)[-(1:4)])*25 # plotting scale = 13.2
sem.model.fits(traitPC1.2012) # marginal R2 = 0.36

# generate model list
mods.12 <- list(traitPC1.2012,traitPC2.2012, soilPC1.2012, soilPC2.2012, arth.rich.2012, arth.abund.2012, arth.rarerich.2012)

# fit SEM
sem.fit(mods.12, data = w.SEM.2012, corr.errors = c("total.abund~~total.rich","total.rich~~total.rarerich", "total.rarerich~~total.abund", "trait.PC1~~trait.PC2","soil.PC1~~soil.PC2")) # adequate fit to the data: Fisher.C = 26.52, df = 22, P = 0.23; AIC = 148.52, AICc = Inf, K = 61, n = 62

mods.12.coefs <- sem.coefs(mods.12, data = w.SEM.2012, corr.errors = c("total.abund~~total.rich","total.rich~~total.rarerich", "total.rarerich~~total.abund", "trait.PC1~~trait.PC2","soil.PC1~~soil.PC2"), standardize = "none")
mods.12.coefs <- mutate(mods.12.coefs, plot.scale = estimate*25)
filter(mods.12.coefs, predictor %in% c("num.Wind","soil.PC1","soil.PC2","trait.PC1","trait.PC2"))

## Wind: arthropod community assembly ----

# only analyzing 2013, because community composition only differed in this year.
w.arth.hell.mech.df <- left_join(filter(wind.arth.df, Year == "2013", total.abund > 1), select(filter(w.trait.df, Year == "2013"), plant_ID, trait.PC1, trait.PC2)) %>% rename(Trait_PC1 = trait.PC1, Trait_PC2 = trait.PC2) %>% mutate(num.Wind = ifelse(Wind.Exposure == "Exposed",1,0))

w.arth.hell <- decostand(w.arth.hell.mech.df[ ,w.arth.names], method = "hellinger")

# test effect of trait PCs
w.arth.hell.rda <- rda(w.arth.hell ~ Trait_PC1 + Trait_PC2, data = w.arth.hell.mech.df)
summary(w.arth.hell.rda)
vif(w.arth.hell.rda)

anova(w.arth.hell.rda, by = "margin", permutations = how(block = w.arth.hell.mech.df$Block, nperm = 999))

w.arth.hell.mech.p <- autoplot.custom(w.arth.hell.rda, scaling = 3, color = "grey") + scale_shape_manual(values = 21) + scale_color_manual(values = "grey") + theme(legend.position = "none") + xlab("RDA 1 (15%)") + ylab("RDA 2 (1%)") + scale_x_continuous(limits = c(-1.6,1.2)); w.arth.hell.mech.p

save_plot("fig_w_arth_comm_mech.png", w.arth.hell.mech.p, base_height = 5, base_width = 5)

# get plot-level centroids after controlling for trait effects
w.arth.hell.plots <- rda(w.arth.hell ~ Condition(Trait_PC1) + Condition(Trait_PC2) + Plot_code, data = w.arth.hell.mech.df)

w.arth.hell.plots.df <- data.frame(scores(w.arth.hell.plots, display = "cn", choices = 1:10), Plot_code = levels(as.factor(w.arth.hell.mech.df$Plot_code))) %>% separate(Plot_code, into = c("Block","Wind.Exposure"))

# test effect of wind exposure after accounting for trait PC1
w.arth.hell.rda.wind <- rda(w.arth.hell.plots.df[ ,1:10] ~ Condition(Block) + Wind.Exposure, data = w.arth.hell.plots.df)
anova(w.arth.hell.rda.wind, by = "margin",permuations = how(block = w.arth.hell.plots.df$Block, nperm = 999))  # suggest that wind exposure did not have a direct effect on community composition (after accounting for indirect effects on plant traits).

## Wind: fungal community assembly ----
fungal.df.mech <- left_join(fungal.df, select(filter(w.trait.df, Year == "2013"), plant_ID, root_CN)) %>% left_join(., select(w.soil, Plot_code, soil.PC1, soil.PC2)) 

fungal.df.mech.rootCNsub <- filter(fungal.df.mech, root_CN > 0)

f.hell.mech <- decostand(fungal.df.mech[ ,f.OTUs], method = "hellinger")
f.hell.mech.rootCNsub <- decostand(fungal.df.mech.rootCNsub[ ,f.OTUs], method = "hellinger")

# to test soil effects, we need to first get the centroid at the plot level
w.f.plots.hell <- betadisper(vegdist(f.hell.mech, method = "euclidean"),  fungal.df.mech$Plot_code, bias.adjust = TRUE)

w.f.plots.centr.hell <- data.frame(w.f.plots.hell$centroids, 
                                   id = rownames(w.f.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))
w.f.plots.centr.hell <- left_join(w.f.plots.centr.hell, select(w.soil, Block, Wind.Exposure, soil.PC1, soil.PC2)) %>% mutate(Block = as.factor(Block), Wind.Exposure = as.factor(Wind.Exposure))

# test effect of soil PCs
w.f.mech.soil.rda <- rda(w.f.plots.hell$centroids ~ Condition(Block) + soil.PC1 + soil.PC2, data = w.f.plots.centr.hell) # no effect of soil.PC2
anova(w.f.mech.soil.rda, by = "margin", permutations = how(block = w.f.plots.centr.hell$Block, nperm = 999)) # no effect of soil PCs

adonis(w.f.plots.hell$centroids ~ Block + soil.PC1 + soil.PC2, data = w.f.plots.centr.hell, permutations = how(block = w.f.plots.centr.hell$Block, nperm = 999), method = "euclidean") # this is proof of concept that the F-test on soil.PC2 is the exact same as the RDA.

# test effect of root_CN (type 2) and effect of Genotype (type 1). Used adonis because it was much faster than RDA and should give the same results with Euclidean distance. Put root_CN last though so the test would effectively be the same as type 2 SS for this test.
#w.f.mech.rootCN.rda <- rda(f.hell.mech.rootCNsub ~ Condition(soil.PC1) + Condition(soil.PC2) + root_CN + Genotype, data = fungal.df.mech.rootCNsub) 
adonis(f.hell.mech.rootCNsub ~ soil.PC1 + soil.PC2 + root_CN, data = fungal.df.mech.rootCNsub, permutations = how(block = fungal.df.mech.rootCNsub$Block, nperm = 999), method = "euclidean") # important to have root_CN in last position to make sure the test preserves marginality. Note that the tests for soil.PC1 and soil.PC2 are bogous because they aren't based on the appropriate degrees of freedom. Effect of root C:N is marginal.

adonis(f.hell.mech.rootCNsub ~ soil.PC1 + soil.PC2 + root_CN + Genotype, data = fungal.df.mech.rootCNsub, permutations = how(block = fungal.df.mech.rootCNsub$Block, nperm = 999), method = "euclidean") # test whether Genotype still has an effect after accounting for these factors.nGenotype effect is significant after accounting for all of the factors, suggesting we don't know what is mediating the assembly of the fungal community.

#anova(w.f.mech.rootCN.rda, by = "term", permutations = how(block = fungal.df.mech.rootCNsub$Block, nperm = 999)) # 

#w.f.mech.rootCN.p <- autoplot.custom(w.f.mech.rootCN.rda, scaling = 3, color = "grey") + scale_shape_manual(values = 21) + scale_color_manual(values = "grey") + theme(legend.position = "none") + xlab("RDA 1 (15%)") + ylab("RDA 2 (1%)") + scale_x_continuous(limits = c(-1.6,1.2)); w.f.mech.rootCN.p 

## Wind: Bacteria community assembly ----
bacteria.df.mech <- left_join(bacteria.df, select(filter(w.trait.df, Year == "2013"), plant_ID, root_CN)) %>% left_join(., select(w.soil, Plot_code, soil.PC1, soil.PC2)) 

bacteria.df.mech.rootCNsub <- filter(bacteria.df.mech, root_CN > 0)

b.hell.mech <- decostand(bacteria.df.mech[ ,b.OTUs], method = "hellinger")
b.hell.mech.rootCNsub <- decostand(bacteria.df.mech.rootCNsub[ ,b.OTUs], method = "hellinger")

# to test soil effects, we need to first get the centroid at the plot level
w.b.plots.hell <- betadisper(vegdist(b.hell.mech, method = "euclidean"),  bacteria.df.mech$Plot_code, bias.adjust = TRUE)

w.b.plots.centr.hell <- data.frame(w.b.plots.hell$centroids, 
                                   id = rownames(w.b.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))
w.b.plots.centr.hell <- left_join(w.b.plots.centr.hell, select(w.soil, Block, Wind.Exposure, soil.PC1, soil.PC2)) %>% mutate(Block = as.factor(Block), Wind.Exposure = as.factor(Wind.Exposure))

# test effect of soil PCs
w.b.soil.rda <- rda(w.b.plots.hell$centroids ~ Condition(Block) + soil.PC1 + soil.PC2, data = w.b.plots.centr.hell) 

anova(w.b.soil.rda, by = "margin",  permutations = how(block = w.b.plots.centr.hell$Block, nperm = 999)) # nothing sig.

# test effect of root_CN (type 2) and effect of Genotype (type 1). Used adonis because it was much faster than RDA and should give the same results with Euclidean distance. Put root_CN last though so the test would effectively be the same as type 2 SS for this test. 
adonis(b.hell.mech.rootCNsub ~ soil.PC1 + soil.PC2 + root_CN, data = bacteria.df.mech.rootCNsub, method = "euclidean", permutations = how(block = bacteria.df.mech.rootCNsub$Block, nperm = 999)) # effect of root C:N is marginal (explains 1% of variance)

## Ant-aphid: SEM data management ----
aa.arth.df <- read.csv('~/Lanphere_Experiments/final_data/ant_aphid_arthropod_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), fact.Ant.mound.dist = as.factor(Ant.mound.dist)) %>% select(-X)
aa.arth.names <- colnames(select(aa.arth.df, Gracilliaridae_miner:Spider))

aa.trait.df <- read.csv('~/Lanphere_Experiments/final_data/ant_aphid_trait_df.csv') %>% tbl_df() %>% #mutate(fact.Ant.Mound.Dist = as.factor(Ant.mound.dist), Plot_code = paste(Block, Ant.mound.dist, sep = "_"), Block = as.factor(Block)) %>% 
  filter(Year == "2012") %>% select(plant_ID, trait.PC1, trait.PC2)

aa.mech.df <- left_join(aa.arth.df, aa.trait.df)

aa.SEM.df <- aa.mech.df %>% mutate(num.Aphid = ifelse(Aphid.treatment == "aphid",2,1)) %>% select(Block, Genotype, fact.Ant.mound.dist, Ant.mound.dist, num.Aphid, aphid_Aphis, ant_F_obscuripes, total.abund, total.rarerich, total.rich, trait.PC1, trait.PC2)

aa.vars <- colnames(select(aa.SEM.df, Ant.mound.dist:trait.PC2))
aa.SEM.df[ ,aa.vars] <- scale(aa.SEM.df[ ,aa.vars])

p.scale <- 15

## Ant-aphid SEM: richness ----
# We only observed an effect of willow genotype on arthropod richness, therefore, we only modelled the indirect effects of willow genotype on richness mediated through plant traits.

# fit separate models
aa.arth.rich.2012 <- lmerTest::lmer(total.rich ~ trait.PC1 + trait.PC2 + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df)
sem.model.fits(aa.arth.rich.2012) # marginal R2 = 0.14
fixef(aa.arth.rich.2012)
fixef(aa.arth.rich.2012)*p.scale

aa.arth.rich.2012.Gmiss <- lmerTest::lmer(total.rich ~ trait.PC1 + trait.PC2 + Genotype + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
sd(fixef(aa.arth.rich.2012.Gmiss)[-(1:3)]) # missing Genotype SD = 0.38
sd(fixef(aa.arth.rich.2012.Gmiss)[-(1:3)])*p.scale # plot scale = 5.6

aa.traitPC1.2012 <- lmerTest::lmer(trait.PC1 ~ Genotype + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
sd(fixef(aa.traitPC1.2012)[-1]) # Genotype SD = 0.39
sd(fixef(aa.traitPC1.2012)[-1])*p.scale # plot scale = 9.8
sem.model.fits(aa.traitPC1.2012) # marginal R2 = 0.12

aa.traitPC2.2012 <- lmerTest::lmer(trait.PC2 ~ Genotype + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
sd(fixef(aa.traitPC2.2012)[-1]) # Genotype SD = 0.47
sd(fixef(aa.traitPC2.2012)[-1])*p.scale # plot scale = 11.6
sem.model.fits(aa.traitPC2.2012) # marginal R2 = 0.23

aa.rich.mods <- list(aa.arth.rich.2012, aa.traitPC1.2012, aa.traitPC2.2012)

sem.fit(aa.rich.mods, aa.SEM.df, corr.errors = c("trait.PC1~~trait.PC2")) # missing path between genotype and total richness. 
aa.rich.mods.coefs <- sem.coefs(aa.rich.mods, aa.SEM.df, corr.errors = c("trait.PC1~~trait.PC2"), standardize = "none")
aa.rich.mods.coefs <- mutate(aa.rich.mods.coefs, plot.scale = estimate*p.scale)
filter(aa.rich.mods.coefs, predictor %in% c("trait.PC1","trait.PC2","~~ trait.PC2")) # mainly effect of trait PC1 rather than PC2

## Ant-aphid SEM: rarefied richness ----
# only observed an effect of aphid treatment so we modelled the direct and indirect (via F obscuripes) on rarefied richness
aa.arth.rarerich.2012 <- lmerTest::lmer(total.rarerich ~ ant_F_obscuripes + aphid_Aphis  + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df) 
sem.model.fits(aa.arth.rarerich.2012) # marginal R2 = 0.01
fixef(aa.arth.rarerich.2012)

aa.arth.rarerich.2012.treat <- lmerTest::lmer(total.rarerich ~ ant_F_obscuripes + aphid_Aphis + num.Aphid  + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df) 
fixef(aa.arth.rarerich.2012.treat) # -0.18
fixef(aa.arth.rarerich.2012.treat)*p.scale # plot scale = 2.8

aa.aphis.2012 <- lmerTest::lmer(aphid_Aphis ~ num.Aphid + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df)
sem.model.fits(aa.aphis.2012) # marginal R2 = 0.19

aa.Fobs.2012 <- lmerTest::lmer(ant_F_obscuripes ~ aphid_Aphis + (1|Block),# + (1|Block:fact.Ant.mound.dist), zero variance
                               aa.SEM.df)
sem.model.fits(aa.Fobs.2012) # marginal R2 = 0.18. Note that adding in aphid_Aphis removed Genotype*num.Aphid effect.

aa.rarerich.mods <- list(aa.arth.rarerich.2012, aa.aphis.2012, aa.Fobs.2012)

sem.fit(aa.rarerich.mods, aa.SEM.df) # missing path between aphid treatment and rarefied richness
aa.rarerich.mods.coefs <- data.frame(sem.coefs(aa.rarerich.mods, aa.SEM.df, standardize = "none"))
aa.rarerich.mods.coefs <- mutate(aa.rarerich.mods.coefs, plot.scale = estimate*p.scale); aa.rarerich.mods.coefs

## Ant-aphid SEM: abundance ----
# modelled direct and indirect effects of both aphid treatment and plant traits.

aa.arth.abund.2012 <- lmerTest::lmer(total.abund ~ ant_F_obscuripes + aphid_Aphis + trait.PC1 + trait.PC2 + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df)
sem.model.fits(aa.arth.abund.2012) # marginal R2 = 0.06

aa.aphis.2012.treats <- lmerTest::lmer(aphid_Aphis ~ num.Aphid*Ant.mound.dist + trait.PC1 + trait.PC2 + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df) # added Ant.mound.dist because this was one of the interactive effects driving total abundance
sem.model.fits(aa.aphis.2012.treats) # marginal R2 = 0.19

aa.aphis.2012.treats.Gmiss <- lmerTest::lmer(aphid_Aphis ~ num.Aphid*Ant.mound.dist + trait.PC1 + trait.PC2 + Genotype + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
sd(fixef(aa.aphis.2012.treats.Gmiss)[-(c(1:5,15))]) # Genotype SD = 0.28
sd(fixef(aa.aphis.2012.treats.Gmiss)[-(c(1:5,15))])*p.scale # plot scale = 4.3

aa.Fobs.2012.dist <- lmerTest::lmer(ant_F_obscuripes ~ aphid_Aphis + Ant.mound.dist + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df)
sem.model.fits(aa.Fobs.2012) # marginal R2 = 0.18

# list of mods
aa.abund.mods <- list(aa.arth.abund.2012, aa.aphis.2012.treats, aa.Fobs.2012.dist, aa.traitPC1.2012, aa.traitPC2.2012)

sem.fit(aa.abund.mods, aa.SEM.df, corr.errors = c("trait.PC1~~trait.PC2")) # missing path between genotype and aphid abundance as well as ant.mound.dist:num.Aphid on total abundance
# Note however, that this Fisher test is inaccurate because the test of the interaction num.Aphid*Ant.mound.dist does not incorporate the main effects (therefore, doesn't preserve marginality). To account for this, I'm going to calculate the p-values of the interaction effect after accounting for both main effects in the models for trait.PC1, trait.PC2, and total.abund

aa.arth.abund.2012.add <- lmerTest::lmer(total.abund ~ ant_F_obscuripes + aphid_Aphis + trait.PC1 + trait.PC2 + num.Aphid*Ant.mound.dist + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df)
summary(aa.arth.abund.2012.add) # p-values of 0.000944 after accounting for significance of main effects
0.188881*p.scale # 2.83 plot scale

aa.traitPC1.2012.add <- lmerTest::lmer(trait.PC1 ~ Genotype + num.Aphid*Ant.mound.dist + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
summary(aa.traitPC1.2012.add) # 0.854631, interactive effect is not significant after accounting for main effects (maintaining marginality)

aa.traitPC2.2012.add <- lmerTest::lmer(trait.PC2 ~ Genotype + num.Aphid*Ant.mound.dist + (1|Block) + (1|Block:fact.Ant.mound.dist), aa.SEM.df, contrasts = list(Genotype = "contr.sum"))
summary(aa.traitPC2.2012.add) # 0.782435, interaction effect is not significant after accounting for main effects (maintaining marginality)

aa.miss <- sem.missing.paths(aa.abund.mods, aa.SEM.df, corr.errors = c("trait.PC1~~trait.PC2"), conditional = TRUE)
aa.miss$p.value[c(8:9,11)] <- c(0.854631, 0.782435, 0.000944) # replacing with p-values that maintain marginality of interaction effect
aa.SEM.C <- -2*sum(log(aa.miss$p.value)) # Fisher C = 48.87655 
1 - pchisq(aa.SEM.C, 2 * length(aa.miss$p.value)) # df = 32, P = 0.02851434

aa.abund.mods.coefs <- sem.coefs(aa.abund.mods, aa.SEM.df, corr.errors = c("trait.PC1~~trait.PC2"), standardize = "none")
aa.abund.mods.coefs <- mutate(aa.abund.mods.coefs, plot.scale = estimate*p.scale)
filter(aa.abund.mods.coefs, predictor %in% c("trait.PC1","trait.PC2","aphid_Aphis","Ant.mound.dist","num.Aphid","ant_F_obscuripes","num.Aphid:Ant.mound.dist"))

## Ant-aphid: miscellaneous stats to aid interpretation ----
round(colSums(aa.mech.df[ ,aa.arth.names])/sum(colSums(aa.mech.df[ ,aa.arth.names]))*100,0)
corr.test(aa.mech.df[ ,c("aphid_Aphis", "ant_F_obscuripes", aa.arth.names, "total.abund","total.rich","total.rarerich")])
with(aa.mech.df, cor.test(Formica_ant, aphid_Aphis))
with(aa.mech.df, cor.test((Syrphidae + Formica_ant + Spider), aphid_Aphis)) # higher correlation with just Formica_ants

## Ant-aphid: arthropod composition mechanisms ----
aa.arth.hell.mech.df <- filter(aa.mech.df, total.abund > 1) %>% rename(A_farinosa_abund = aphid_Aphis, F_obscuripes_abund = ant_F_obscuripes, Trait_PC1 = trait.PC1, Trait_PC2 = trait.PC2)

aa.arth.hell <- decostand(aa.arth.hell.mech.df[ ,aa.arth.names], method = "hellinger")

# test effect of F obscuripes 
aa.hell.mech.rda <- rda(aa.arth.hell ~ A_farinosa_abund + F_obscuripes_abund + Trait_PC1 + Trait_PC2, data = aa.arth.hell.mech.df) # nonsig.
vif(aa.hell.mech.rda) # multicollinearity shouldn't be a problem
summary(aa.hell.mech.rda)
anova(aa.hell.mech.rda, by = "margin", permutations = how(block = aa.arth.hell.mech.df$Block, nperm = 999))

aa.hell.mech.p <- autoplot.custom(aa.hell.mech.rda, scaling = 3, color = "grey") + scale_shape_manual(values = 21) + scale_color_manual(values = "grey") + theme(legend.position = "none") + xlab("RDA 1 (2%)") + ylab("RDA 2 (1%)") + scale_x_continuous(limits = c(-0.8,1.3))

save_plot("fig_aa_comm_mech.png", aa.hell.mech.p, base_height = 5, base_width = 5)

# test remainder effect of genotype*aphid interaction
RsquareAdj(update(aa.hell.mech.rda, .~. + Genotype*Aphid.treatment))
anova(update(aa.hell.mech.rda, .~. + Genotype*Aphid.treatment), by = "margin", permutations = how(block = aa.arth.hell.mech.df$Block, nperm = 999)) # only interpreting GxE effect. 
anova(update(aa.hell.mech.rda, .~. + Genotype + Aphid.treatment), by = "margin", permutations = how(block = aa.arth.hell.mech.df$Block, nperm = 999)) # only interpreting G and Aphid.treatment effects.

