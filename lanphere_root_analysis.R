## Analysis of root-associated communities

## load required libraries ----
source('~/miscellaneous_R/model_diagnostic_functions.R')
source('~/miscellaneous_R/veganCovEllipse.R')
library(adespatial) #source('~/miscellaneous_R/beta.div.R')
library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(dplyr)
library(tidyr)
library(vegan)
library(lme4)
library(ggplot2)
library(cowplot) # throwing an error for unknown reason
library(effects)

## load require data ----

## Fungal community (all)
f.sub.OTUs <- read.csv("normalized_fungi_otu_table.csv") %>% select(-X) %>% colnames()

f.taxa <- read.csv("fungi_taxa_table.csv") %>% tbl_df() %>% rename(OTU_ID = X)

fungal.df <- read.csv("fungal.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich)))

fungal.df$fungal.sub.abund <- rowSums(round(fungal.df[ ,f.sub.OTUs]),0)
#fungal.df$fungal.sub.PIE <- rarefy(round(fungal.df[ ,f.sub.OTUs],0), 2) - 1 # Hurlbert's pie. 
#fungal.df$fungal.sub.rarerich <- rarefy(round(fungal.df[ ,f.sub.OTUs],0), min(fungal.df$fungal.sub.abund))

# arbuscular mycorrhizal OTUs
abm <- read.csv("~/Lanphere_Experiments/fungi_FUNGUILD_AMF.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(abm)[c(1,137:141)] <- c("OTU_ID","taxon","taxon_level","trophic_mode","guild","confidence_ranking")
abm.char <- select(abm, OTU_ID, taxonomy:confidence_ranking)

# ectomycorrhizal OTUs
ecm <- read.csv("~/Lanphere_Experiments/fungi_FUNGUILD_ECM.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(ecm)[c(1,137:142)] <- c("OTU_ID","taxon","taxon_level","trophic_mode","guild","confidence_ranking","growth_morphology")
ecm.char <- select(ecm, OTU_ID, taxonomy:growth_morphology)

# pathogen OTUs
path <- read.csv("~/Lanphere_Experiments/fungi_FUNGUILD_pathogens.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(path)[c(1,137:143)] <- c("OTU_ID","taxon","taxon_level","trophic_mode","guild","confidence_ranking","growth_morphology","trait")
path.char <- select(path, OTU_ID, taxonomy:trait)


# bacteria community
bacteria.df <- read.csv("bacteria.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
#b.OTUs <- colnames(select(bacteria.df, Crenarchaeota:MVP.21))
#b.OTUs <- colnames(select(bacteria.df, Thaumarchaeota:BD1.5))
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich)))


## Fungal abundance analysis ----

#fung.abund.lmer <- lmer(fungal.abund ~ Wind.Exposure*Genotype +
 #                         (1|Block) + (1|Block:Wind.Exposure),
  #                      data = fungal.df,
   #                     contrasts = list(Wind.Exposure = "contr.sum",
    #                                     Genotype = "contr.sum"))
#summary(fung.abund.lmer)
#plot(fung.abund.lmer) # looks out without transformation

# no effects on fungal abundance
#(fungal.anova <- anova.table(fung.abund.lmer, type = 2, experiment = "wind"))

#fungal.abund.up <- update(fung.abund.lmer, .~. -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))
#anova(fungal.abund.up, update(fungal.abund.up, .~. -(1|Block:Wind.Exposure)))

#(fungal.abund.R2 <- var.table(fungal.abund.up, experiment = "wind"))

## Fungal richness analysis ----
#fung.rich.lmer <- lmer(fungal.rich ~ Wind.Exposure*Genotype +
 #                        (1|Block) + (1|Block:Wind.Exposure),
  #                     data = fungal.df,
   #                    contrasts = list(Wind.Exposure = "contr.sum",
    #                                    Genotype = "contr.sum"))
#summary(fung.rich.lmer)
#plot(fung.rich.lmer) # looks out without transformation

# no effects on fungal richness
#(fungal.rich.anova <- anova.table(fung.rich.lmer, type = 2, experiment = "wind"))

#fungal.rich.up <- update(fung.rich.lmer, .~. -Wind.Exposure*Genotype -Genotype + Wind.Exposure + (1|Genotype))
#plotREsim(REsim(fungal.rich.up)) # not working...

#(fungal.rich.R2 <- var.table(fungal.rich.up, experiment = "wind"))

## Fungal rarefied richness analysis ----
# not that I get the same results if I round everything and use Poisson that accounts for overdispersion. Residuals look a bit better with linear model...

hist(fungal.df$fungal.rarerich)
fung.rarerich.lmer <- lmer(fungal.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = fungal.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))
summary(fung.rarerich.lmer)
plot(fung.rarerich.lmer) # looks out without transformation
#overdisp_fun(fung.rarerich.lmer)

# no effects on fungal rarerichness
#drop1(fung.rarerich.lmer, test = "Chisq")
#drop1(update(fung.rarerich.lmer, .~. -Wind.Exposure:Genotype), test = "Chisq")

(fungal.rarerich.anova <- anova.table(fung.rarerich.lmer, type = 2, experiment = "wind"))

fungal.rarerich.up <- update(fung.rarerich.lmer, .~. -Wind.Exposure*Genotype -Genotype + Wind.Exposure + (1|Genotype))
#plotREsim(REsim(fungal.rarerich.up)) # not working...

(fungal.rarerich.R2 <- var.table(fungal.rarerich.up, experiment = "wind"))

## Bacteria abundance analysis ----
#bact.abund.lmer <- lmer(log(bacteria.abund) ~ Wind.Exposure*Genotype +
 #                         (1|Block) + (1|Block:Wind.Exposure),
  #                      data = bacteria.df,
   #                     contrasts = list(Wind.Exposure = "contr.sum",
    #                                     Genotype = "contr.sum"))
#summary(bact.abund.lmer)
#plot(bact.abund.lmer) # looks out without transformation

# no effects on bacteria abundance
#anova.table(bact.abund.lmer, type = 2, experiment = "wind")

## Bacteria richness analysis ----
#bact.rich.lmer <- lmer(log(bacteria.rich) ~ Wind.Exposure*Genotype +
 #                        (1|Block) + (1|Block:Wind.Exposure),
  #                     data = bacteria.df,
   #                    contrasts = list(Wind.Exposure = "contr.sum",
    #                                    Genotype = "contr.sum"))
#summary(bact.rich.lmer)
#plot(bact.rich.lmer) # looks better with transformation

#plot(Effect("Wind.Exposure", bact.rich.lmer, transformation = list(link = log, inverse = exp))) # 1.1 fold increase or about 10%

# no effects on bacteria richness. Marginal effect of wind exposure
#anova.table(bact.rich.lmer, type = 2, experiment = "wind")

## Bacteria rarefied richness analysis ----
bact.rarerich.lmer <- lmer(bacteria.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = bacteria.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))
summary(bact.rarerich.lmer)
plot(bact.rarerich.lmer) 

# no effects on bacteria rarerichness. Marginal effect of wind exposure
anova.table(bact.rarerich.lmer, type = 2, experiment = "wind")

plot(Effect("Wind.Exposure", bact.rarerich.lmer)) 
as.data.frame(Effect("Wind.Exposure", bact.rarerich.lmer))
1153.492/1085.967 # 6% higher rarefied richness on exposed plants.

## Total Fungal Community analysis: hellinger ----
key.OTUs <- which(f.sub.OTUs %in% names(which(f.beta$SCBD > 0.004)))
other.OTUs <- which(f.sub.OTUs %in% c("OTU_54","OTU_68","OTU_134","OTU_66","OTU_119","OTU_262","OTU_1606","OTU_113"))
names(which(f.beta$SCBD > 0.004))
sum(f.beta$SCBD[which(f.beta$SCBD > 0.004)]) # together these species explain about 6% of the beta-diversity among genotypes.

car::scatterplotMatrix(log(mean.fcomm[ ,f.sub.OTUs[key.OTUs]]+1))

f.hell <- decostand(fungal.df[ ,f.sub.OTUs[-other.OTUs]], method = "hellinger") # -key.OTUs


# test G and GxE
adonis(f.hell ~ Wind.Exposure*Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no GxE

adonis(f.hell ~ Wind.Exposure + Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # but a G

f.hell.geno <- betadisper(vegdist(f.hell, method = "euclidean"), fungal.df$Genotype, bias.adjust = TRUE)
boxplot(f.hell.geno) # doesn't seem like there is too much variation in dispersion
plot(f.hell.geno)
permutest(f.hell.geno, permutations = how(block = fungal.df$Plot_code, nperm = 999)) # marginaly significant overdispersion, suggesting genotype effect could be due to the variance in dispersion.

# test wind effect
w.f.plots.hell <- betadisper(vegdist(f.hell, method = "euclidean"),  fungal.df$Plot_code, bias.adjust = TRUE)

w.f.plots.centr.hell <- data.frame(w.f.plots.hell$centroids, 
                                   id = rownames(w.f.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.f.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.f.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.f.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Which species are driving fungal community differences among genotypes?
mean.fcomm <- fungal.df %>% 
  group_by(Block, Genotype) %>% # summarise by Genotype at block-level first
  summarise_each(funs(mean)) %>% 
  group_by(Genotype) %>%        # then by Genotype
  summarise_each(funs(mean))

f.beta <- beta.div(mean.fcomm[ ,f.sub.OTUs])
hist(f.beta$SCBD)
names(which(f.beta$SCBD > 0.002))

OTU.abund <- colSums(fungal.df[ ,f.sub.OTUs]) # total abundance
OTU.mean.beta <- f.beta$SCBD
plot(log(OTU.mean.beta) ~ log(OTU.abund), type = 'n')
text(log(OTU.mean.beta) ~ log(OTU.abund), label = names(OTU.abund))

#plot(OTU.all.beta ~ OTU.abund, type = 'n')
#text(OTU.all.beta ~ OTU.abund, label = names(OTU.abund))

as.matrix(mean.fcomm[ ,names(which(f.beta$SCBD > 0.004))])
colSums(fungal.df[ ,names(which(f.beta$SCBD > 0.004))])

taxa.sub <- which(as.character(f.taxa$OTU_ID) %in% c("OTU_54","OTU_68","OTU_134","OTU_66","OTU_119","OTU_262","OTU_1606","OTU_113"))#names(which(f.beta$SCBD > 0.004)) == TRUE)

table(f.taxa$Guild)
as.data.frame(f.taxa[taxa.sub, ]) %>% select(OTU_ID, Species)
write.csv(as.data.frame(f.taxa[taxa.sub, ]), "key_fungal_taxa.csv")

hist(fungal.df$OTU_54)
library(glmmADMB)

hist(colSums(fungal.df[ ,f.sub.OTUs]))
summary(colSums(fungal.df[ ,f.sub.OTUs]))
sort(colSums(fungal.df[ ,f.sub.OTUs]))

# def. sig: 54,
# sig: 313, 146, 135, 104 (but errors)
# not sig: 171, 199, 115, 82
OTU.lmer <- glmer(round(OTU_113,0) ~ offset(log(fungal.sub.abund)) + Wind.Exposure*Genotype + (1|X) + (1|Block) + (1|Block:Wind.Exposure), data = fungal.df, family = "poisson", contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
overdisp_fun(OTU.lmer)
plot(OTU.lmer)
summary(OTU.lmer)
drop1(OTU.lmer, test = "Chisq")
drop1(update(OTU.lmer, .~. -Wind.Exposure:Genotype), test = "Chisq")
# def. sig: 54, 68, 
#anova.table(OTU.lmer, type = 2, experiment = "wind")

ectos <- f.taxa$Guild %in% "Ectomycorrhizal"
ecto.names <- as.character(f.taxa$OTU_ID[which(ectos == TRUE)])

plot(OTU_10/fungal.sub.abund ~ Wind.Exposure, fungal.df)
plot(OTU_68/fungal.sub.abund ~ Genotype, fungal.df)
plot(OTU_134/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_54/fungal.sub.abund  ~ Genotype, fungal.df)
plot(f.hell$OTU_54  ~ Genotype, fungal.df)
plot(OTU_119/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_262/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_113/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_82/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_66/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_11/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_68/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_171/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_199/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_6/fungal.sub.abund  ~ Genotype, fungal.df)
plot(OTU_5/fungal.sub.abund  ~ Genotype, fungal.df)

library(glmmADMB)
library(lme4)
hist(fungal.df$OTU_199)
fungal.df$X <- as.factor(fungal.df$X)
plot(OTU_54 ~ Genotype, fungal.df)
OTU.glmer <- glmer(round(OTU_134,0) ~ Wind.Exposure + Genotype + 
                        (1|X) +
                        (1|Block) + (1|Block:Wind.Exposure),
                      data = fungal.df,
                   family = "poisson",
                      contrasts = list(Wind.Exposure = "contr.sum",
                                       Genotype = "contr.sum"),
                      control=glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e5)))
summary(OTU.glmer)
#overdisp_fun(OTU.glmer)
#plot(OTU.glmer)
drop1(OTU.glmer, test = "Chisq")


## Arbuscular mycorrhizal Community analysis: hellinger ----
abm.hell <- decostand(fungal.df[ ,abm.char$OTU_ID], method = "hellinger")

abm.rda <- rda(abm.hell ~ Condition(Block) + Condition(Wind.Exposure) + Genotype, data = fungal.df)
plot(abm.rda, display = "sp")
OTU.lmer <- glmer(round(OTU_107,0) ~ offset(log(fungal.sub.abund)) + Wind.Exposure  + Genotype + (1|X) + (1|Block) + (1|Block:Wind.Exposure), data = fungal.df, family = "poisson", contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
overdisp_fun(OTU.lmer)
plot(OTU.lmer)
summary(OTU.lmer)
drop1(OTU.lmer, test = "Chisq")
fungal.df$OTU_107



# test G and GxE
adonis(abm.hell ~ Wind.Exposure*Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no GxE

adonis(abm.hell ~ Wind.Exposure + Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no G

abm.hell.geno <- betadisper(vegdist(abm.hell, method = "euclidean"), fungal.df$Genotype, bias.adjust = TRUE)
boxplot(abm.hell.geno) 
plot(abm.hell.geno)
permutest(abm.hell.geno, permutations = how(block = fungal.df$Plot_code, nperm = 999)) # no clear evidence of overdispersion.

# test wind effect
w.abm.plots.hell <- betadisper(vegdist(abm.hell, method = "euclidean"),  fungal.df$Plot_code, bias.adjust = TRUE)

w.abm.plots.centr.hell <- data.frame(w.abm.plots.hell$centroids, 
                                   id = rownames(w.abm.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.abm.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.abm.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.abm.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Ectomycorrhizal Community analysis: hellinger ----
ecm.hell <- decostand(fungal.df[ ,ecm.char$OTU_ID], method = "hellinger")

ecm.rda <- rda(ecm.hell ~ Condition(Block) + Condition(Wind.Exposure) + Genotype, data = fungal.df)
plot(ecm.rda, display = "sp", type = "t")
OTU.lmer <- glmer(round(OTU_39,0) ~ offset(log(fungal.sub.abund)) + Wind.Exposure  + Genotype + (1|X) + (1|Block) + (1|Block:Wind.Exposure), data = fungal.df, family = "poisson", contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
overdisp_fun(OTU.lmer)
plot(OTU.lmer)
summary(OTU.lmer)
drop1(OTU.lmer, test = "Chisq")
fungal.df$OTU_1819 # marginal effect of Genotype
fungal.df$OTU_311 # marginal effect of Genotype
fungal.df$OTU_8 # strong effect of wind
fungal.df$OTU_39 # effect of Genotype
plot(OTU_39/fungal.sub.abund ~ Genotype, fungal.df) # this doesn't look good.

# test G and GxE
adonis(ecm.hell ~ Wind.Exposure*Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no GxE

adonis(ecm.hell ~ Wind.Exposure + Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no G

ecm.hell.geno <- betadisper(vegdist(ecm.hell, method = "euclidean"), fungal.df$Genotype, bias.adjust = TRUE)
boxplot(ecm.hell.geno) 
plot(ecm.hell.geno)
permutest(ecm.hell.geno, permutations = how(block = fungal.df$Plot_code, nperm = 999)) # marginal evidence of overdispersion.

# test wind effect
w.ecm.plots.hell <- betadisper(vegdist(ecm.hell, method = "euclidean"),  fungal.df$Plot_code, bias.adjust = TRUE)

w.ecm.plots.centr.hell <- data.frame(w.ecm.plots.hell$centroids, 
                                     id = rownames(w.ecm.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.ecm.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.ecm.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.ecm.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Pathogen Community analysis: hellinger ----
path.hell <- decostand(fungal.df[ ,path.char$OTU_ID], method = "hellinger")

# test G and GxE
adonis(path.hell ~ Wind.Exposure*Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no GxE

adonis(path.hell ~ Wind.Exposure + Genotype, data = fungal.df, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999)) # no G

path.hell.geno <- betadisper(vegdist(path.hell, method = "euclidean"), fungal.df$Genotype, bias.adjust = TRUE)
boxplot(path.hell.geno) 
plot(path.hell.geno)
permutest(path.hell.geno, permutations = how(block = fungal.df$Plot_code, nperm = 999)) # no evidence of overdispersion.

# test wind effect
w.path.plots.hell <- betadisper(vegdist(path.hell, method = "euclidean"),  fungal.df$Plot_code, bias.adjust = TRUE)

w.path.plots.centr.hell <- data.frame(w.path.plots.hell$centroids, 
                                     id = rownames(w.path.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.path.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.path.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.path.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Bacteria Community analysis: hellinger distance ----
b.hell <- decostand(bacteria.df[ ,b.OTUs], method = "hellinger")
#b.comm.weff <- select(bacteria.df.weff, -Block, -Wind.Exposure, -Plot_code)

# test GxE
adonis(b.hell ~ Wind.Exposure*Genotype, 
       data = bacteria.df, 
       method = "euclidean", 
       permutations = how(block = bacteria.df$Block, nperm = 999)) # no GxE

# test G
adonis(b.hell ~ Wind.Exposure + Genotype, 
       data = bacteria.df, 
       method = "euclidean", 
       permutations = how(block = bacteria.df$Block, nperm = 999)) # no GxE


# test wind effect
w.b.plots.hell <- betadisper(vegdist(b.hell, method = "euclidean"),  bacteria.df$Plot_code, bias.adjust = TRUE)

w.b.plots.centr.hell <- data.frame(w.b.plots.hell$centroids, 
                                   id = rownames(w.b.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.b.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.b.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.b.plots.centr.hell$Block, nperm = 999)) # marginal effect of wind exposure

b.hell.wind <- betadisper(vegdist(w.b.plots.hell$centroids, method = "euclidean"), w.b.plots.centr.hell$Wind.Exposure, bias.adjust = TRUE)
boxplot(b.hell.wind)
plot(b.hell.wind)
permutest(b.hell.wind, permutations = how(block = w.b.plots.centr.hell$Block, nperm = 999)) # no difference in disperion 

## Mantel tests for correlation in community composition ----
common <- bacteria.df$plant_ID[bacteria.df$plant_ID %in% fungal.df$plant_ID] %>% as.character()

b.hell.sub <- bacteria.df %>% filter(plant_ID %in% common) %>% .[ ,b.OTUs] %>% decostand(., method = "hellinger")

f.hell.sub <- fungal.df %>% filter(plant_ID %in% common) %>% .[ ,f.OTUs] %>% decostand(., method = "hellinger")

sub.strata <- fungal.df %>% filter(plant_ID %in% common) %>% select(Plot_code)
sub.strata = as.numeric(sub.strata$Plot_code)

# strong, positive correlation in community composition
mantel(vegdist(f.hell.sub, method = "euclidean"), vegdist(decostand(b.hell.sub, method = "hellinger"), method = "euclidean"), strata = sub.strata)

## Plots: Fungal and Bacteria community ----

# general settings and required functions
pd <- 0.2
cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")
p.size <- 5
l.size <- 1.5

## Richness plots
f.rich.df <- data.frame(Effect(c("Wind.Exposure","Genotype"), fung.rich.lmer))
levels(f.rich.df$Wind.Exposure) <- c("Exposed","Unexposed","NA")
f.rich.p <- ggplot(f.rich.df, aes(x = Wind.Exposure, y = fit, color = Genotype)) + geom_line(aes(group = Genotype), size = l.size) + geom_point(size = p.size)  + scale_color_manual(values = cbPal.10) + ylab("No. of OTUs") + xlab("Wind exposure") + theme(legend.position = "none") + ggtitle("Mycorrhiza")

b.rich.df <- data.frame(Effect(c("Wind.Exposure","Genotype"), bact.rich.lmer, transformation = list(link = log, inverse = exp)))
levels(b.rich.df$Wind.Exposure) <- c("Exposed","Unexposed","NA")
b.rich.p <- ggplot(b.rich.df, aes(x = Wind.Exposure, y = fit, color = Genotype))+ geom_line(aes(group = Genotype), size = l.size) + geom_point(size = p.size)  + scale_color_manual(values = cbPal.10) + ylab("No. of OTUs") + xlab("Wind exposure")+ theme(legend.position = "none") + ggtitle("Bacteria") 

## Genotype effect on fungal community
f.hell.geno.rda <- rda(f.hell ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = fungal.df)
summary(f.hell.geno.rda)$cont$importance[ ,1:2] # first 2 axes only explain 2% and 1% of the variance in fungal community composition.

library(ggvegan)
autoplot(f.hell.geno.rda, display = "species")#, display = "sp", scaling = 2, type = "t")
scrs <- data.frame(scores(f.hell.geno.rda, display = "species", scaling = 3))
scrs$guild <- NA
scrs$guild[which(row.names(scrs) %in% path.char$OTU_ID)] <- "pathogen"
scrs$guild[which(row.names(scrs) %in% abm.char$OTU_ID)] <- "abm"
scrs$guild[which(row.names(scrs) %in% ecm.char$OTU_ID)] <- "ecm"
scrs$guild[which(is.na(scrs$guild))] <- "other"
scrs$guild <- as.factor(scrs$guild)

which(ecm.char$OTU_ID == "OTU_54")

#plot.new()
#plot.window(xlim = c(-0.05,0.05), ylim = c(-0.05,0.05), asp = 1)
plot(f.hell.geno.rda, scaling = 3, display = c("cn"))
#abline(h = 0, lty = "dotted")
#abline(v = 0, lty = "dotted")
colvec <- c("red2", "green4","grey", "mediumblue")
with(scrs, points(x = scrs$RDA1, y = scrs$RDA2, col = colvec[scrs$guild], pch = 1))

ellip <- ordiellipse(f.hell.geno.rda, groups = fungal.df$Genotype, col = "gray50",draw = "polygon")

centroids.cap.geno <- data.frame(scores(f.hell.geno.rda, display = "cn"), Genotype = levels(fungal.df$Genotype))

sites.cap.geno <- data.frame(scores(f.hell.geno.rda, choices = c(1,2), display = "sites"), droplevels(fungal.df$Genotype))
colnames(sites.cap.geno)[3] <- "Genotype"

# get data for ellipse. 
df_ell.cap.geno <- data.frame() #sets up a data frame before running the function.
for(g in levels(sites.cap.geno$Genotype)){
  df_ell.cap.geno <- rbind(df_ell.cap.geno, 
                           cbind(as.data.frame(
                             with(sites.cap.geno[sites.cap.geno$Genotype == g, ], 
                                  veganCovEllipse(ellip[[g]]$cov, ellip[[g]]$center, ellip[[g]]$scale))), Genotype = g))
}

f.ord <- ggplot(data = df_ell.cap.geno, aes(x = RDA1, y = RDA2, group = Genotype)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_point(data = sites.cap.geno, color = "gray", shape = 1) +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap.geno, 
            aes(x = RDA1, y = RDA2, label = Genotype), size = 5) +
  ylab("RDA 2 (1%)") + xlab("RDA 1 (2%)") + ggtitle("Mycorrhiza")
f.ord

## Wind effect on bacteria community
levels(w.b.plots.centr.hell$Wind.Exposure) <- c("Exp.","Unexp.")
levels(w.b.plots.centr.hell$Block)
b.hell.wind.rda <- rda(w.b.plots.hell$centroids ~ Block + Wind.Exposure, data = w.b.plots.centr.hell)
summary(b.hell.wind.rda)$cont$importance[ ,1:2] # first axis explains 14% of the variation, second axis explains 10%

ellip <- ordiellipse(b.hell.wind.rda, groups = w.b.plots.centr.hell$Wind.Exposure, draw = "polygon")

centroids.cap <- data.frame(scores(b.hell.wind.rda, display = "cn"), Wind.Exposure = levels(w.b.plots.centr.hell$Wind.Exposure))

sites.cap <- data.frame(scores(b.hell.wind.rda, choices = c(1,2), display = "sites"), droplevels(w.b.plots.centr.hell$Wind.Exposure), droplevels(w.b.plots.centr.hell$Block))
colnames(sites.cap)[3] <- "Wind.Exposure"
colnames(sites.cap)[4] <- "Block"

# get data for ellipse. 
df_ell.cap <- data.frame() #sets up a data frame before running the function.
for(g in levels(sites.cap$Wind.Exposure)){
  df_ell.cap <- rbind(df_ell.cap, 
                      cbind(as.data.frame(
                        with(sites.cap[sites.cap$Wind.Exposure == g, ], 
                             veganCovEllipse(ellip[[g]]$cov, ellip[[g]]$center, ellip[[g]]$scale))), Wind.Exposure = g))
}

b.ord.wind <- ggplot(data = df_ell.cap, aes(x = RDA1, y = RDA2, group = Wind.Exposure)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_text(data = sites.cap, aes(x = RDA1, y = RDA2, label = Block), color = "gray") +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap[11:12, ], 
            aes(x = RDA1, y = RDA2, label = Wind.Exposure), size = 5) +
  ylab("RDA 2 (10%)") + xlab("RDA 1 (14%)") + ggtitle("Bacteria")
b.ord.wind

# create multipanel plot. Run plot scripts for arthropod analysis first though #aa.rGE, w.rGE, f.rich.p, b.rich.p,
(comm.p <- plot_grid(w.arth.ord.wind, f.ord, b.ord.wind, labels = "AUTO", ncol = 3, align = 'h'))

save_plot("fig_5_wind_comps.png", comm.p, base_height = 3, base_width = 8.5)
#save_plot("fig_3_root_comm.png", root.comm.p, base_height = 8, base_width = 8)

################### Likely don't need below ----

## Root C:N analysis ----
rootCN.df$Carbon_mg
rootCN.df$Carbon_percent
rootCN.df$Nitrogen_mg
rootCN.df$Nitrogen_percent
plot(Nitrogen_percent ~ root_CN, filter(rootCN.df, root_CN < 100))

plot(Carbon_percent ~ Genotype, rootCN.df)
plot(Nitrogen_percent ~ Genotype, rootCN.df)
plot(root_CN ~ Genotype, filter(rootCN.df, root_CN < 80))

rootCN.lmer <- lmer(root_CN ~ Wind.Exposure*Genotype + (1|Block) + (1|Block:Wind.Exposure),
                    data = filter(rootCN.df, root_CN < 100), # results are qualitatively the same with the full data set.
                    contrasts = list(Wind.Exposure = "contr.sum",
                                     Genotype = "contr.sum")) # restricting to biologically reasonable values
summary(rootCN.lmer)
plot(rootCN.lmer)

(rootCN.anova <- anova.table(rootCN.lmer, type = 2, experiment = "wind"))

rootCN.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), rootCN.lmer)) %>% mutate(response = "root_CN")
write.csv(rootCN.GxE, "rootCN.means.csv")

rootCN.up <- update(rootCN.lmer, .~. -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))

(rootCN.R2 <- var.table(rootCN.up, experiment = "wind"))
## Plot of G and E effects on bacterial community ----
b.hell.samp <- b.hell[ ,sample(colnames(b.hell), 10)] # subsample 2000 OTUs of the original hellinger distance (after hellinger distance was calculated for original samples). Need to subset the data for analyzing the genotype effect within rda.

## Genotype effect
b.hell.geno.rda <- rda(b.hell.samp ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = bacteria.df)
summary(b.hell.geno.rda)$cont$importance[ ,1:2] 

plot(b.hell.geno.rda, type = "n")
points(b.hell.geno.rda, display = "sites", col = "gray50")
ordiellipse(b.hell.geno.rda, groups = bacteria.df$Genotype, col = "gray50",draw = "polygon")
text(b.hell.geno.rda, display = "cn", labels = levels(bacteria.df$Genotype))

## Wind exposure
b.hell.wind.rda <- rda(w.b.plots.hell$centroids ~ Block + Wind.Exposure, data = w.b.plots.centr.hell)
summary(b.hell.wind.rda)$cont$importance[ ,1:2] # first axis explains 13% of the variation, second axis explains 10%

plot(b.hell.wind.rda, type = "n")
points(b.hell.wind.rda, display = "sites", col = "gray50")
ordiellipse(b.hell.wind.rda, groups = w.b.plots.centr.hell$Wind.Exposure, draw = "polygon", col = "gray50")
text(x = scores(b.hell.wind.rda, display = "cn")[11:12,1], y = scores(b.hell.wind.rda, display = "cn")[11:12,2], labels = c("Exposed","Unexposed"))

# Interaction effect ----
adonis(bacteria.df[ ,bact.comm] ~ Wind.Exposure*Genotype, 
       data = bacteria.df, 
       strata = bacteria.df$Block, 
       method = "horn")

# Genotype effect. No effect
adonis(bacteria.df[ ,bact.comm] ~ Genotype, 
       data = bacteria.df, 
       strata = bacteria.df$Plot_code, 
       method = "horn")

bacteria.permdisp <- betadisper(vegdist(bacteria.df[ ,bact.comm],
                                        method = "horn"), 
                                group = bacteria.df$Genotype, 
                                bias.adjust = TRUE)
anova(bacteria.permdisp) # okay
plot(bacteria.permdisp) # doesn't look like too much differentiation among communities...
boxplot(bacteria.permdisp)

meandist(vegdist(bacteria.df[ ,bact.comm], 
                 method = "horn"), 
         group = bacteria.df$Genotype)

# Wind exposure effect. No effect once I aggregate.
adonis(wind.effect.bact[ ,bact.comm] ~ Wind.Exposure, 
       data = wind.effect.bact, 
       strata = wind.effect.bact$Block, 
       method = "horn")

## Root CN to microbial community response
f.sh <- fungal.df %>% dplyr::select(Block, Wind.Exposure, Genotype, Plot_code, Sample, fungal.abund, fungal.rich)
b.sh <- bacteria.df %>% dplyr::select(Block, Wind.Exposure, Genotype, Plot_code, Sample, bacteria.abund, bacteria.rich)
CN.sh <- rootCN.df %>% dplyr::select(Block, Wind.Exposure, Genotype, Plot_code, Sample, root_CN)

trait.microbe.df <- plyr::join_all(list(f.sh, b.sh, CN.sh), by = c("Block","Wind.Exposure","Genotype","Plot_code","Sample")) %>%
  tbl_df()

trait.soil.microbe.df <- left_join(trait.microbe.df, 
                                   select(w.soil.df, Plot_code, Total.N:avg.moisture.vwc), by = "Plot_code")

corr.test(select(filter(trait.soil.microbe.df, root_CN < 100),
                 Total.N, nut.PC1, percent.TOM, avg.moisture.vwc,
                 root_CN, fungal.abund, fungal.rich, 
                 bacteria.abund, bacteria.rich))

plot(root_CN ~ Total.N + nut.PC1 + percent.TOM + avg.moisture.vwc, filter(trait.soil.microbe.df, root_CN < 100))

plot(fungal.abund ~ Total.N + nut.PC1 + percent.TOM + avg.moisture.vwc, filter(trait.soil.microbe.df, root_CN < 100))
plot(fungal.rich ~ Total.N + nut.PC1 + percent.TOM + avg.moisture.vwc, filter(trait.soil.microbe.df, root_CN < 100))

plot(log(bacteria.abund) ~ Total.N + nut.PC1 + percent.TOM + avg.moisture.vwc, filter(trait.soil.microbe.df, root_CN < 100))
plot(log(bacteria.rich) ~ Total.N + nut.PC1 + percent.TOM + avg.moisture.vwc, filter(trait.soil.microbe.df, root_CN < 100))

ggplot(filter(trait.microbe.df, root_CN < 100),
       aes(x = root_CN, y = log(fungal.abund))) +
  geom_point() +
  stat_smooth(method = "lm")

# effect of soil N on root_CN?
root_CN.mech <- lmer(root_CN ~ Total.N + 
                       (1|Block) + 
                       (1|Block:Wind.Exposure), 
                     data = filter(trait.soil.microbe.df, root_CN < 100))
summary(root_CN.mech)
Anova(root_CN.mech, test = "F")

# apparent effect of CN on fungal abundance and richness
microbe.CN <- lmer(fungal.abund ~ scale(root_CN) +
                     (1|Block) + (1|Block:Wind.Exposure),
                   data = filter(trait.soil.microbe.df, root_CN < 100))
summary(microbe.CN)
plot(microbe.CN)

plot(allEffects(microbe.CN))
plotFEsim(FEsim(microbe.CN))
plotREsim(REsim(microbe.CN))

Anova(microbe.CN, test = "F")

var.table(microbe.CN, experiment = "wind")


## OLD ----
# lmer
library(pbkrtest)
library(lmerTest)

otu14 <- lmer(OTU_1606 ~ Wind.Exposure*Genotype +
                (1|Block) + (1|Block:Wind.Exposure),
              data = fungal.df)
plot(otu14)
summary(otu14)
anova(otu14, ddf = "Kenward-Roger")


## load library for analysis
library(mvabund)
library(vegan)

fungal.mv <- mvabund(select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position))
sort(colSums(fungal.mv))
OTU.common <- which(colSums(fungal.mv) > 800)
fungal.com.mv <- fungal.mv[ ,OTU.common]

meanvar.plot(fungal.mv ~ fungal.df$Wind.Exposure)

plot.mvabund(fungal.com.mv, type = "bx")
plot.mvformula(fungal.com.mv ~ fungal.df$Block, legend = FALSE)

fungal.many <- manyany("lm", fungal.com.mv, 
                       data = fungal.df, 
                       log(fungal.com.mv+1) ~ Wind.Exposure*Genotype + Block, 
                       composition = TRUE)
fungal.null <- manyany("glm", fungal.com.mv, data = fungal.df, fungal.com.mv ~ 1, composition = TRUE, family = poisson(link = "log"))
anova.manyany(fungal.many, fungal.null, test = "LR")
summary(fungal.many, p.uni = "unadjusted")
plot(fungal.many)

fungal.community <- select(fungal.df, -plant_ID, -Wind.Exposure, -Genotype, -Block, -Plant.Position)
OTU.common <- which(colSums(fungal.community) > 800)
fungal.adonis <- adonis(fungal.community[ ,OTU.common] ~ Wind.Exposure*Genotype + Block, fungal.df, 
                        #strata = fungal.df$Block, 
                        method = "horn")
fungal.adonis

fungal.permdisp <- betadisper(vegdist(fungal.community[ ,OTU.common], method = "horn"), group = fungal.df$Wind.Exposure, bias.adjust = TRUE)
anova(fungal.permdisp)
plot(fungal.permdisp)
boxplot(fungal.permdisp)

meandist(vegdist(fungal.community[ ,OTU.common], method = "horn"), group = fungal.df$Block)

## Fungal community: Euclidean distance ----
f.euc <- log(select(fungal.df, -(plant_ID:GxE))+1)
#f.comm.weff <- select(fungal.df.weff, -Block, -Wind.Exposure, -Plot_code)

mantel(vegdist(f.euc, "euclidean"), vegdist(f.hell, "euclidean")) # highly correlated, so we just used the Hellinger transformation. 

# test G and GxE
adonis(f.euc ~ Wind.Exposure*Genotype, 
       data = fungal.df, 
       method = "euclidean", 
       permutations = how(block = fungal.df$Block, nperm = 999)) # no GxE

# test wind effect
w.f.plots.euc <- betadisper(vegdist(f.euc, method = "euclidean"),  fungal.df$Plot_code, bias.adjust = TRUE)

w.f.plots.centr.euc <- data.frame(w.f.plots.euc$centroids, 
                                  id = rownames(w.f.plots.euc$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.f.plots.euc$centroids ~ Block + Wind.Exposure, 
       data = w.f.plots.centr.euc, 
       method = "euclidean", 
       permutations = how(block = w.f.plots.centr.euc$Block, nperm = 999))

#adonis(log(f.comm+1) ~ Wind.Exposure*Genotype, fungal.df, permutations = how(blocks = fungal.df$Block, nperm = 999), method = "euclidean") # test G and GxE. important to put wind exposure first to account for its effects since this design is unbalanced.

#adonis(log(f.comm.weff+1) ~ Block + Wind.Exposure, fungal.df.weff, permutations = how(block = fungal.df.weff$Block, nperm = 999), method = "euclidean") # test wind effect

# Analyze model assumptions
#f.G.euc.permdisp <- betadisper(vegdist(log(f.comm+1), method = "euclidean"), group = fungal.df$Genotype, bias.adjust = TRUE) # Genotype
#boxplot(f.G.euc.permdisp)
#anova(f.G.euc.permdisp) # no difference in dispersion between groups

#f.E.euc.permdisp <- betadisper(vegdist(log(f.comm.weff+1), method = "euclidean"), group = fungal.df.weff$Wind.Exposure, bias.adjust = TRUE) # wind exposure
#boxplot(f.E.euc.permdisp)
#anova(f.E.euc.permdisp)

## Hellinger distance
#adonis(decostand(f.comm, method = "hellinger") ~ Wind.Exposure*Genotype, fungal.df, permutations = how(blocks = fungal.df$Block, nperm = 999), method = "euclidean") # test G and GxE. important to put wind exposure first to account for its effects since this design is unbalanced.

#adonis(decostand(f.comm.weff, method = "hellinger") ~ Block + Wind.Exposure, fungal.df.weff, permutations = how(block = fungal.df.weff$Block, nperm = 999), method = "euclidean") # test wind effect

## Analyze model assumptions

# Genotype
#f.G.permdisp <- betadisper(vegdist(decostand(f.comm, method = "hellinger"), method = "euclidean"), group = fungal.df$Genotype, bias.adjust = TRUE)
#boxplot(f.G.permdisp)
#anova(f.G.permdisp) # no difference in dispersion between groups

# Wind Exposure
#f.E.permdisp <- betadisper(vegdist(decostand(f.comm.weff, method = "hellinger"), method = "euclidean"), group = fungal.df.weff$Wind.Exposure, bias.adjust = TRUE)
#boxplot(f.E.permdisp)
#anova(f.E.permdisp) # no difference in dispersion between groups


## plot Genotype effect
#f.G.rda <- update(f.GE.rda, .~. -Wind.Exposure + Condition(Wind.Exposure))
#plot(f.G.rda, type = "n")
#ordiellipse(f.G.rda, groups = fungal.df$Genotype, draw = "polygon", col = "gray50", kind = "sd")

## plot Wind Exposure effect
#f.E.rda <- update(f.GE.rda, .~. -Genotype + Condition(Genotype))
#plot(f.E.rda, type = "n")
#points(f.E.rda, display = c("sites"), col = as.numeric(fungal.df$Wind.Exposure), pch = 21, bg = as.numeric(fungal.df$Wind.Exposure))
#ordiellipse(f.E.rda, groups = fungal.df$Wind.Exposure, draw = "polygon", col = "gray50", kind = "sd")

## Bacteria community analysis: euclidean distance ----
b.euc <- log(select(bacteria.df, -(plant_ID:GxE))+1)
#b.comm.weff <- select(bacteria.df.weff, -Block, -Wind.Exposure, -Plot_code)

# test GxE. Note how the two models are equivalent in their tests of the interaction term
adonis(b.euc ~ Wind.Exposure*Genotype, 
       data = bacteria.df, 
       method = "euclidean", 
       permutations = how(block = bacteria.df$Block, nperm = 999)) # no GxE


# test wind effect
w.b.plots.euc <- betadisper(vegdist(b.euc, method = "euclidean"),  bacteria.df$Plot_code, bias.adjust = TRUE)

w.b.plots.centr.euc <- data.frame(w.b.plots.euc$centroids, 
                                  id = rownames(w.b.plots.euc$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.b.plots.euc$centroids ~ Block + Wind.Exposure, 
       data = w.b.plots.centr.euc, 
       method = "euclidean", 
       permutations = how(block = w.b.plots.centr.euc$Block, nperm = 999)) # marginal effect of wind exposure

b.euc.wind <- betadisper(vegdist(w.b.plots.euc$centroids, method = "euclidean"), w.b.plots.centr.euc$Wind.Exposure, bias.adjust = TRUE)
boxplot(b.euc.wind)
plot(b.euc.wind)
anova(b.euc.wind)



## Euclidean distance
#adonis(log(b.comm+1) ~ Wind.Exposure*Genotype, bacteria.df, permutations = how(blocks = bacteria.df$Block, nperm = 999), method = "euclidean") # test G and GxE. important to put wind exposure first to account for its effects since this design is unbalanced.

#adonis(log(b.comm.weff+1) ~ Block + Wind.Exposure, bacteria.df.weff, permutations = how(block = bacteria.df.weff$Block, nperm = 999), method = "euclidean") # test wind effect

# Analyze model assumptions
#b.G.euc.permdisp <- betadisper(vegdist(log(b.comm+1), method = "euclidean"), group = bacteria.df$Genotype, bias.adjust = TRUE) # Genotype
#boxplot(b.G.euc.permdisp)
#anova(b.G.euc.permdisp) # marginal difference in dispersion between groups

#b.E.euc.permdisp <- betadisper(vegdist(log(b.comm.weff+1), method = "euclidean"), group = bacteria.df.weff$Wind.Exposure, bias.adjust = TRUE) # wind exposure
#boxplot(b.E.euc.permdisp)
#anova(b.E.euc.permdisp)

## Hellinger distance
#adonis(decostand(b.comm, method = "hellinger") ~ Wind.Exposure*Genotype, bacteria.df, permutations = how(blocks = bacteria.df$Block, nperm = 999), method = "euclidean") # test G and GxE. important to put wind exposure first to account for its effects since this design is unbalanced.

#adonis(decostand(b.comm.weff, method = "hellinger") ~ Block + Wind.Exposure, bacteria.df.weff, permutations = how(block = bacteria.df.weff$Block, nperm = 999), method = "euclidean") # test wind effect

## Analyze model assumptions

# Genotype
#b.G.permdisp <- betadisper(vegdist(decostand(b.comm, method = "hellinger"), method = "euclidean"), group = bacteria.df$Genotype, bias.adjust = TRUE)
#boxplot(b.G.permdisp)
#anova(b.G.permdisp) # no difference in dispersion between groups

# Wind Exposure
#b.E.permdisp <- betadisper(vegdist(decostand(b.comm.weff, method = "hellinger"), method = "euclidean"), group = bacteria.df.weff$Wind.Exposure, bias.adjust = TRUE)
#boxplot(b.E.permdisp)
#anova(b.E.permdisp) # no difference in dispersion between groups

## plot Genotype effect
#b.G.rda <- update(b.GE.rda, .~. -Wind.Exposure + Condition(Wind.Exposure))
#plot(b.G.rda, type = "n")
#ordiellipse(b.G.rda, groups = bacteria.df$Genotype, draw = "polygon", col = "gray50", kind = "sd")

## plot Wind Exposure effect
#b.E.rda <- update(b.GE.rda, .~. -Genotype + Condition(Genotype))
#plot(b.E.rda, type = "n")
#points(b.E.rda, display = c("sites"), col = as.numeric(bacteria.df$Wind.Exposure), pch = 21, bg = as.numeric(bacteria.df$Wind.Exposure))
#ordiellipse(b.E.rda, groups = bacteria.df$Wind.Exposure, draw = "polygon", col = "gray50", kind = "sd")

## Upload and manage root C:N data ----

rootCN <- read.csv('Output/RootCNdataLanphereDunes.csv') %>%
  tbl_df() 

# need to split up plant ID into block, treatment, genotype and position
rootCN$TreatGeno <- gsub("[[:digit:]]", "", rootCN$Sample)
rootCN$Block__Position <- gsub("[^[:digit:]]", "_", rootCN$Sample)

rootCN.df <- rootCN %>%
  separate(TreatGeno, into = c("Wind.Exposure","Genotype"), sep = 1) %>%
  filter(Genotype %in% c("F","G","I","J","L","S","T","U","W","X")) %>% 
  separate(Block__Position, into = c("Block","Plant.Position")) %>%
  mutate(Genotype = as.factor(Genotype),
         Wind.Exposure = as.factor(Wind.Exposure),
         Block = as.factor(Block),
         Plot_code = interaction(Block, Wind.Exposure),
         Sample = interaction(Block, Wind.Exposure, Genotype, Plant.Position),
         root_CN = Carbon_mg/Nitrogen_mg) %>%
  filter(Nitrogen_mg > 0 & Carbon_mg > 0)
