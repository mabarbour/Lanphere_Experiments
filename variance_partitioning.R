## VARIANCE PARTITIONING GXE AND TRAIT CONTRIBUTIONS

library(tidyverse)
library(brms)

## WIND DATA

wind.arth.df <- read.csv('final_data/wind_arthropod_df.csv') %>%
  tbl_df() %>% 
  mutate(Block = as.factor(Block), Year = as.factor(Year), Plot_code = paste(Block, Wind.Exposure, sep = "."))

w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Year = as.factor(Year), Plot_code = paste(Block, Wind.Exposure, sep = "."))
glimpse(w.trait.df)

wind.data.2012 <- left_join(select(wind.arth.df, -X), select(w.trait.df, -X)) %>%
  filter(Year == "2012")
contrasts(wind.data.2012$Wind.Exposure) <- "contr.sum"
wind.data.2012$trait.PC1 <- scale(wind.data.2012$trait.PC1)
wind.data.2012$trait.PC2 <- scale(wind.data.2012$trait.PC2)

## ARTHROPOD RICHNESS 2012

w.rich.2012 <- brm(total.rich ~ Wind.Exposure + (Wind.Exposure | Genotype) + (1|Block) + (1|Block:Wind.Exposure), 
                   data = wind.data.2012, family = "poisson", control = list(adapt_delta = 0.95))
bayes_R2(w.rich.2012)

fixef.rich.2012 <- as.data.frame(fixef(w.rich.2012))
fixef.model <- model.matrix(total.rich ~ Wind.Exposure, data = wind.data.2012)
fixef.var <- var(as.vector(fixef.rich.2012$Estimate %*% t(fixef.model)))

re.var <- as.data.frame(VarCorr(w.rich.2012, old = TRUE))
re.var$Std.Dev^2
dist.var <- log(1/exp(fixef.rich.2012$Estimate[1])+1)

(fixef.var + sum(re.var$Std.Dev^2))/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # conditional = 29%
fixef.var /(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 7% wind
re.var$Std.Dev[1]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 3% Block
re.var$Std.Dev[2]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 4% Block:WindExposure
re.var$Std.Dev[3]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 14% Genotype
re.var$Std.Dev[4]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 1% GxE



w.traits.rich.2012 <- brm(total.rich ~ trait.PC1 + trait.PC2 + Wind.Exposure + (Wind.Exposure | Genotype) + (1|Block) + (1|Block:Wind.Exposure), 
                          data = wind.data.2012, family = "poisson", control = list(adapt_delta = 0.95))
bayes_R2(w.traits.rich.2012)

fixef.traits.rich.2012 <- as.data.frame(fixef(w.traits.rich.2012))
fixef.traits.model <- model.matrix(total.rich ~ trait.PC1 + trait.PC2 + Wind.Exposure, data = wind.data.2012)
fixef.traits.var <- var(as.vector(fixef.traits.rich.2012$Estimate %*% t(fixef.traits.model)))
fixef.traits.var

wind.var <- var(as.vector(fixef.traits.rich.2012$Estimate[c(1,4)] %*% t(fixef.traits.model[,c(1,4)])))
trait.PC1.var <- var(as.vector(fixef.traits.rich.2012$Estimate[c(1,2)] %*% t(fixef.traits.model[,c(1,2)])))
trait.PC2.var <- var(as.vector(fixef.traits.rich.2012$Estimate[c(1,3)] %*% t(fixef.traits.model[,c(1,3)])))
wind.var + trait.PC1.var + trait.PC2.var

re.var.traits <- as.data.frame(VarCorr(w.traits.rich.2012, old = TRUE))
re.var.traits$Std.Dev^2
dist.var.traits <- log(1/exp(fixef.traits.rich.2012$Estimate[1])+1)

(fixef.traits.var + sum(re.var.traits$Std.Dev^2))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # conditional = 29%
(wind.var + trait.PC1.var + trait.PC2.var + sum(re.var.traits$Std.Dev^2))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # conditional for partitioning = 26%
wind.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 4% wind
trait.PC1.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 4% trait PC1
trait.PC2.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # < 1% trait PC2
re.var.traits$Std.Dev[1]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 2% Block
re.var.traits$Std.Dev[2]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 3% Block:WindExposure
re.var.traits$Std.Dev[3]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 12% Genotype
re.var.traits$Std.Dev[4]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 1% GxE



## ARTHROPOD RICHNESS 2013

wind.data.2013 <- left_join(select(wind.arth.df, -X), select(w.trait.df, -X)) %>%
  filter(Year == "2013")
contrasts(wind.data.2013$Wind.Exposure) <- "contr.sum"
wind.data.2013$trait.PC1 <- scale(wind.data.2013$trait.PC1)
wind.data.2013$trait.PC2 <- scale(wind.data.2013$trait.PC2)

w.rich.2013 <- brm(total.rich ~ Wind.Exposure + (Wind.Exposure | Genotype) + (1|Block) + (1|Block:Wind.Exposure), 
                   data = wind.data.2013, family = "poisson", control = list(adapt_delta = 0.95))
bayes_R2(w.rich.2013)

fixef.rich.2013 <- as.data.frame(fixef(w.rich.2013))
fixef.model <- model.matrix(total.rich ~ Wind.Exposure, data = wind.data.2013)
fixef.var <- var(as.vector(fixef.rich.2013$Estimate %*% t(fixef.model)))

re.var <- as.data.frame(VarCorr(w.rich.2013, old = TRUE))
re.var$Std.Dev^2
dist.var <- log(1/exp(fixef.rich.2013$Estimate[1])+1)

(fixef.var + sum(re.var$Std.Dev^2))/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # conditional = 33%
fixef.var /(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 14% wind
re.var$Std.Dev[1]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 9% Block
re.var$Std.Dev[2]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 5% Block:WindExposure
re.var$Std.Dev[3]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 4% Genotype
re.var$Std.Dev[4]^2/(fixef.var + sum(re.var$Std.Dev^2) + dist.var) # 1% GxE



w.traits.rich.2013 <- brm(total.rich ~ trait.PC1 + trait.PC2 + Wind.Exposure + (Wind.Exposure | Genotype) + (1|Block) + (1|Block:Wind.Exposure), 
                          data = wind.data.2013, family = "poisson", control = list(adapt_delta = 0.95))
bayes_R2(w.traits.rich.2013)

fixef.traits.rich.2013 <- as.data.frame(fixef(w.traits.rich.2013))
fixef.traits.model <- model.matrix(total.rich ~ trait.PC1 + trait.PC2 + Wind.Exposure, data = wind.data.2013)
fixef.traits.var <- var(as.vector(fixef.traits.rich.2013$Estimate %*% t(fixef.traits.model)))
fixef.traits.var

wind.var <- var(as.vector(fixef.traits.rich.2013$Estimate[c(1,4)] %*% t(fixef.traits.model[,c(1,4)])))
trait.PC1.var <- var(as.vector(fixef.traits.rich.2013$Estimate[c(1,2)] %*% t(fixef.traits.model[,c(1,2)])))
trait.PC2.var <- var(as.vector(fixef.traits.rich.2013$Estimate[c(1,3)] %*% t(fixef.traits.model[,c(1,3)])))
wind.var + trait.PC1.var + trait.PC2.var

re.var.traits <- as.data.frame(VarCorr(w.traits.rich.2013, old = TRUE))
re.var.traits$Std.Dev^2
dist.var.traits <- log(1/exp(fixef.traits.rich.2013$Estimate[1])+1)

(fixef.traits.var + sum(re.var.traits$Std.Dev^2))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # conditional = 33%
(wind.var + trait.PC1.var + trait.PC2.var + sum(re.var.traits$Std.Dev^2))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # conditional for partitioning = 29%

all.fix.R2 <- (fixef.traits.var)/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits)
wind.R2 <- wind.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 6% wind
trait.PC1.R2 <- trait.PC1.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 11% trait PC1
trait.PC2.R2 <- trait.PC2.var/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 1% trait PC2

re.var.traits$Std.Dev[1]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 5% Block
re.var.traits$Std.Dev[2]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 3% Block:WindExposure
re.var.traits$Std.Dev[3]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 2% Genotype
re.var.traits$Std.Dev[4]^2/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits) # 1% GxE

## COMMONALITY PRACTICE ANALYSIS
U_trait.PC2 <- all.fix.R2 - var(as.vector(fixef.traits.rich.2013$Estimate[c(1,2,4)] %*% t(fixef.traits.model[,c(1,2,4)])))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits)
U_wind <- all.fix.R2 - var(as.vector(fixef.traits.rich.2013$Estimate[c(1,2,3)] %*% t(fixef.traits.model[,c(1,2,3)])))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits)
U_trait.PC1 <- all.fix.R2 - var(as.vector(fixef.traits.rich.2013$Estimate[c(1,3,4)] %*% t(fixef.traits.model[,c(1,3,4)])))/(fixef.traits.var + sum(re.var.traits$Std.Dev^2) + dist.var.traits)
Cij <- 
Cik
Cjk
Cijk

