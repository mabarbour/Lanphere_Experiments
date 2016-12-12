## Update Figure 1

## load required libraries ----
source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(psych)
library(lme4)
library(car)
library(vegan)
library(effects)
library(ggplot2)
library(cowplot)

## upload datasets ----

## ant-aphid: arthropod community
aa.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/ant_aphid_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(
    X = as.factor(X),
    Block = as.factor(Block),
    fact.Ant.mound.dist = as.factor(Ant.mound.dist),
    ord.Ant.mound.dist = ordered(Ant.mound.dist),
    GxE = C(interaction(Genotype,Aphid.treatment), "contr.sum")) 

aa.arth.names <- colnames(select(aa.arth.df, Gracilliaridae_miner:Spider))

# subset of data where plants had at least one arthropod individual
aa.arth.12.pos <- aa.arth.df %>%
  filter(total.abund > 1) 

## wind: arthropod community
# dead plants have already been removed
wind.arth.df <- read.csv('~/Documents/Lanphere_Experiments/final_data/wind_arthropod_df.csv') %>%
  tbl_df() %>%
  mutate(X = as.factor(X),
         Block = as.factor(Block),
         Year = as.factor(Year),
         Plot_code = interaction(Block, Wind.Exposure),
         GxE = C(interaction(Genotype, Wind.Exposure), contr = "contr.sum", how.many = 9))  # create a new variable for the interaction to permit testing of main effects with type 3 sum of squares. 

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

w.13.hell.pos2 <- decostand(w.arth.13.pos2[ ,wind.arth.names], method = "hellinger")

## fungi
fungal.df <- read.csv("fungal.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich)))

## bacteria
bacteria.df <- read.csv("bacteria.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich)))

## Wind: arthropod abundance analysis ----
arth.abund.glmer <- glmer(total.abund ~ Wind.Exposure*Genotype*Year +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            (1|X) +
                            (1|plant_ID),
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

arth.abund.up <- update(arth.abund.glmer, .~. -Wind.Exposure*Genotype*Year + Wind.Exposure + Year + (1|Genotype))

(arth.abund.R2 <- var.table(arth.abund.up, experiment = "wind"))

# get useful data frames
w.w.ab <- as.data.frame(Effect("Wind.Exposure", arth.abund.glmer)) %>% mutate(Wind.Exposure = as.factor(ifelse(Wind.Exposure == "Exposed", "Exposed", "Control"))) %>% rename(Treatment = Wind.Exposure) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

w.g.ab <- as.data.frame(Effect("Genotype", arth.abund.glmer)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))

w.ab <- rbind.data.frame(w.w.ab, w.g.ab)

## Wind: rarefied arthropod richness analysis ----
# need to filter data so that there are 2 or more individuals
arth.rarerich.glmer <- lmer(total.rarerich ~ Wind.Exposure*Year+ Wind.Exposure*Genotype +
                              (1|Block) + 
                              (1|Block:Wind.Exposure) +
                              (1|plant_ID),
                            data = filter(wind.arth.df, total.abund > 1),
                            contrasts = list(Wind.Exposure = "contr.sum",
                                             Genotype = "contr.sum",
                                             Year = "contr.sum"))

arth.rarerich.up <- update(arth.rarerich.glmer, .~. -Wind.Exposure*Genotype -Wind.Exposure*Year + Wind.Exposure + Year + (1|Genotype))

(arth.rarerich.R2 <- var.table(arth.rarerich.up, experiment = "wind"))

# get useful data frames
w.w.rr <- as.data.frame(Effect("Wind.Exposure", arth.rarerich.glmer)) %>% mutate(Wind.Exposure = as.factor(ifelse(Wind.Exposure == "Exposed", "Exposed", "Control"))) %>% rename(Treatment = Wind.Exposure) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

w.g.rr <- as.data.frame(Effect("Genotype", arth.rarerich.glmer)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))

w.rr <- rbind.data.frame(w.w.rr, w.g.rr)

## Wind effect on community composition ----
w.13.plots.hell <- betadisper(vegdist(w.13.hell.pos2, method = "euclidean"),  w.arth.13.pos2$Plot_code, bias.adjust = TRUE)

w.13.plots.centr.hell <- data.frame(w.13.plots.hell$centroids, 
                                    id = rownames(w.13.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

levels(w.13.plots.centr.hell$Wind.Exposure) <- c("Exp.","Unexp.")
w.hell.wind.rda <- rda(w.13.plots.hell$centroids ~ Block + Wind.Exposure, data = w.13.plots.centr.hell) 
summary(w.hell.wind.rda)$cont$importance[ ,1:2] # first axis explains 46% of the variation, second axis explains 7%

plot.new()
ellip <- ordiellipse(w.hell.wind.rda, groups = w.13.plots.centr.hell$Wind.Exposure, draw = "polygon")

centroids.cap <- data.frame(scores(w.hell.wind.rda, display = "cn"), Wind.Exposure = levels(w.13.plots.centr.hell$Wind.Exposure))

sites.cap <- data.frame(scores(w.hell.wind.rda, choices = c(1,2), display = "sites"), droplevels(w.13.plots.centr.hell$Wind.Exposure), droplevels(w.13.plots.centr.hell$Block))
colnames(sites.cap)[3] <- "Wind.Exposure"
colnames(sites.cap)[4] <- "Block"

source('~/Documents/miscellaneous_R/veganCovEllipse.R')

# get data for ellipse. 
df_ell.cap <- data.frame() #sets up a data frame before running the function.
for(g in levels(sites.cap$Wind.Exposure)){
  df_ell.cap <- rbind(df_ell.cap, 
                      cbind(as.data.frame(
                        with(sites.cap[sites.cap$Wind.Exposure == g, ], 
                             veganCovEllipse(ellip[[g]]$cov, ellip[[g]]$center, ellip[[g]]$scale))), Wind.Exposure = g))
}

w.arth.ord.wind <- ggplot(data = df_ell.cap, aes(x = RDA1, y = RDA2, group = Wind.Exposure)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_text(data = sites.cap, aes(x = RDA1, y = RDA2, label = Block), color = "gray") +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap[11:12, ], 
            aes(x = RDA1, y = RDA2, label = Wind.Exposure), size = 5) +
  ylab("") + xlab("")
  #ylab("RDA 2 (7%)") + xlab("RDA 1 (46%)") #+ ggtitle("Arthropods")
w.arth.ord.wind

## fungus rarefied rich ----
fung.rarerich.lmer <- lmer(fungal.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = fungal.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))

# get useful data frames
w.w.frr <- as.data.frame(Effect("Wind.Exposure", fung.rarerich.lmer)) %>% mutate(Wind.Exposure = as.factor(ifelse(Wind.Exposure == "Exposed", "Exposed", "Control"))) %>% rename(Treatment = Wind.Exposure) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

w.g.frr <- as.data.frame(Effect("Genotype", fung.rarerich.lmer)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))

w.frr <- rbind.data.frame(w.w.frr, w.g.frr)

## Genotype effect on fungal community ----
f.hell <- decostand(fungal.df[ ,f.OTUs], method = "hellinger")

f.hell.geno.rda <- rda(f.hell ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = fungal.df)
summary(f.hell.geno.rda)$cont$importance[ ,1:2] # first 2 axes only explain 2% and 1% of the variance in fungal community composition.

#plot.new()
#plot.window(xlim = c(-0.05,0.05), ylim = c(-0.05,0.05), asp = 1)
#plot(f.hell.geno.rda, scaling = 3, display = c("cn"))
#abline(h = 0, lty = "dotted")
#abline(v = 0, lty = "dotted")
#colvec <- c("red2", "green4","grey", "mediumblue")
#with(scrs, points(x = scrs$RDA1, y = scrs$RDA2, col = colvec[scrs$guild], pch = 1))

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
  ylab("") + xlab("")
  #ylab("RDA 2 (1%)") + xlab("RDA 1 (2%)") + ggtitle("Rhizosphere fungi")
f.ord

## Bacteria rarefied richness analysis ----
bact.rarerich.lmer <- lmer(bacteria.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = bacteria.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))

fixef(bact.rarerich.lmer)
var.table(bact.rarerich.lmer, experiment = "wind")

# get useful data frames
w.w.brr <- as.data.frame(Effect("Wind.Exposure", bact.rarerich.lmer)) %>% mutate(Wind.Exposure = as.factor(ifelse(Wind.Exposure == "Exposed", "Exposed", "Control"))) %>% rename(Treatment = Wind.Exposure) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

w.g.brr <- as.data.frame(Effect("Genotype", bact.rarerich.lmer)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))

w.brr <- rbind.data.frame(w.w.brr, w.g.brr)

## Wind effect on bacteria community ----
b.hell <- decostand(bacteria.df[ ,b.OTUs], method = "hellinger")

w.b.plots.hell <- betadisper(vegdist(b.hell, method = "euclidean"),  bacteria.df$Plot_code, bias.adjust = TRUE)

w.b.plots.centr.hell <- data.frame(w.b.plots.hell$centroids, 
                                   id = rownames(w.b.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

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
  ylab("") + xlab("RDA Axis 1")
  #ylab("RDA 2 (10%)") + xlab("RDA 1 (14%)") #+ ggtitle("Bacteria")
b.ord.wind

## ant-aphid: arthropod abundance analysis ----

aa.arth.abund.fact <- glmer(total.abund ~ ord.Ant.mound.dist*Aphid.treatment + Genotype + (1|plant_ID) + (1|Block/fact.Ant.mound.dist), data = aa.arth.df, contrasts = list(Aphid.treatment = "contr.sum", Genotype = "contr.sum"), family = "poisson", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # model ant mound distance as an ordered factor to make plotting easier. Same qualitative result as when I model it as a continuous predictor. Also, I removed interaction with willow genotype since it was not significant and it allowed the model to converge.

# get useful data frames
a.aa.ab <- as.data.frame(Effect(c("Aphid.treatment"), aa.arth.abund.fact)) %>% mutate( Aphid.treatment = as.factor(ifelse(Aphid.treatment == "aphid","Aphids","Control"))) %>% rename(Treatment = Aphid.treatment) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)])) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

#a.aa.ab <- as.data.frame(Effect(c("Aphid.treatment","ord.Ant.mound.dist"), aa.arth.abund.fact)) %>% mutate(ord.Ant.mound.dist = ordered(as.numeric(as.character(ord.Ant.mound.dist))), Aphid.treatment = as.factor(ifelse(Aphid.treatment == "aphid","Aphids","None"))) %>% rename(Treatment = Aphid.treatment, Dummy = ord.Ant.mound.dist) %>% rbind.data.frame(data.frame(Treatment = "", Dummy = "", fit = NA, se = NA, lower = NA, upper = NA))

a.g.ab <- as.data.frame(Effect("Genotype", aa.arth.abund.fact)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))#, Dummy = "") %>% select(Treatment, Dummy, fit:upper)

a.ab <- rbind.data.frame(a.aa.ab, a.g.ab)
#aa.aExE <- as.data.frame(Effect(c("Aphid.treatment","ord.Ant.mound.dist"), aa.arth.abund.fact)) %>% mutate(ord.Ant.mound.dist = ordered(as.numeric(as.character(ord.Ant.mound.dist))), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) 

#%>% ggplot(aes(x = ord.Ant.mound.dist, y = fit, fill = Aphid.treatment, shape = Aphid.treatment)) + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd, position = position_dodge(width = 0.5)) + geom_point(size = p.size, position = position_dodge(width = 0.5)) + ylab("No. of individuals") + xlab("Distance from ant mound (m)") + scale_shape_manual(name = "Aphid treatment", values = c(23,21)) + scale_fill_manual(name = "Aphid treatment", values = c("gray50", "white")) + scale_y_continuous(name = "Arthropod abundance", breaks = c(0,2,4,6,8,10), limits = c(0,10.5)) + theme(legend.justification = c(0,1), legend.position = c(0,1)); aa.aExE

## ant-aphid: rarefied arthropod richness analysis ----
aa.arth.rarerich.glmer <- lmer(total.rarerich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype + (1|Block/fact.Ant.mound.dist), data = filter(aa.arth.df, total.abund > 1), contrasts = list(Aphid.treatment = "contr.sum", Genotype = "contr.sum"))

# get useful data frames
a.a.rr <- as.data.frame(Effect(c("Aphid.treatment"), aa.arth.rarerich.glmer)) %>% mutate(Aphid.treatment = as.factor(ifelse(Aphid.treatment == "aphid","Aphids","Control"))) %>% rename(Treatment = Aphid.treatment) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit, decreasing = TRUE)])) %>% rbind.data.frame(data.frame(Treatment = "", fit = NA, se = NA, lower = NA, upper = NA))

a.g.rr <- as.data.frame(Effect(c("Genotype"), aa.arth.rarerich.glmer)) %>% rename(Treatment = Genotype) %>% mutate(Treatment = factor(Treatment, levels = Treatment[order(fit)]))

a.rr <- rbind.data.frame(a.a.rr, a.g.rr)

## Ant-aphid: ordination of community composition ----
aa.hell <- decostand(aa.arth.12.pos[ ,aa.arth.names], method = "hellinger")

aa.hell.GxE.rda <- rda(aa.hell ~ Condition(Block) + Condition(Ant.mound.dist) +
                         GxE, 
                       data = aa.arth.12.pos)
summary(aa.hell.GxE.rda)$cont$importance[ ,1:2] 

aa.hell.GxE.df <- data.frame(scores(aa.hell.GxE.rda, display = "cn")) %>%
  mutate(factor = levels(aa.arth.12.pos$GxE)) %>%
  separate(factor, into = c("Genotype","Aphid.treatment")) %>%
  mutate(Genotype = as.factor(Genotype),
         Aphid.treatment = as.factor(Aphid.treatment),
         GxE.sig.lab = ifelse(Genotype == "J","yes","no"))

aa.hell.GxE.df.arrows <- data.frame(x = aa.hell.GxE.df$RDA1[11:20], y = aa.hell.GxE.df$RDA2[11:20], xend = aa.hell.GxE.df$RDA1[1:10], yend = aa.hell.GxE.df$RDA2[1:10], Genotype = aa.hell.GxE.df$Genotype[11:20]) %>% mutate(GxE.sig.lab = ifelse(Genotype == "J","yes","no"))

aa.hell.GxE.sites <- data.frame(scores(aa.hell.GxE.rda, display = "sites"))

(aa.arth.ord.GxE <- ggplot(data = aa.hell.GxE.sites, aes(x = RDA1, y = RDA2)) +
  #coord_fixed() + # important for maintaining aspect ratio and distance interpretation; however, this doesn't appear to be replicating the aspect ratio with the vegan plot method...
  geom_point(shape = 1, color = "gray") +
  geom_segment(data = aa.hell.GxE.df.arrows, aes(x = x, xend = xend, y = y, yend = yend, linetype = GxE.sig.lab), color = "black") + #, color = Genotype
  geom_point(data = aa.hell.GxE.df, aes(x = RDA1, y = RDA2, fill = Aphid.treatment, shape = Aphid.treatment), color = "black", size = 5) + 
  geom_text(data = aa.hell.GxE.df, 
            aes(x = RDA1, y = RDA2, label = Genotype), size = 3) +
  scale_fill_manual(values = c("gray50","white"))+ #cbPal.10) +
  scale_color_manual(values = c("gray50","white")) +# cbPal.10) +
  scale_shape_manual(values = c(23,21)) +
  scale_linetype_manual(values = c("dotted","solid"), guide = "none") +
  ggtitle("Composition") + 
  ylab("") + xlab("") +
  #ylab("RDA 2 (3%)") + xlab("RDA 1 (9%)") + 
  theme(legend.position = "none"))

## Plots ----
pd <- 0.75
ebar.w <- 0.05
l.size <- 1.5
alp <- 0.5
p.size <- 3

w.ab.p <- ggplot(w.ab, aes(x = Treatment, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Wind\nFoliar Arthropods", limits = c(0,3.1), breaks = c(0,1,2,3)) + xlab("Experimental Treatments") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)) +
  geom_vline(xintercept = 3, linetype = "dashed"); w.ab.p

w.rr.p <- w.rr %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Treatment, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0,0.75), breaks = c(0,0.25,0.5,0.75)) + xlab("") + geom_vline(xintercept = 3, linetype = "dashed") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)); w.rr.p

w.frr.p <- w.frr %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Treatment, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Wind\nRoot Fungi") + xlab("") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)) +
  geom_vline(xintercept = 3, linetype = "dashed"); w.frr.p

w.brr.p <- w.brr %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Treatment, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Wind\nRoot Bacteria") + xlab("Experimental Treatments") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)) +
  geom_vline(xintercept = 3, linetype = "dashed"); w.brr.p

a.ab.p <- ggplot(a.ab, aes(x = Treatment, y = fit)) + 
  geom_point(size = p.size) + 
  #geom_point(size = p.size, position = position_dodge(width = pd), fill = "white", shape = 21) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Ant-Aphid\nFoliar Arthropods", breaks = c(0,2,4,6,8,10), limits = c(0,10.5)) + xlab("") + geom_vline(xintercept = 3, linetype = "dashed") + 
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)) +
  #geom_text(aes(label = Dummy), position = position_dodge(width = pd)) + 
  ggtitle("Total Abundance"); a.ab.p

a.rr.p <- a.rr %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Treatment, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "") + xlab("") + geom_vline(xintercept = 3, linetype = "dashed") + ggtitle("Rarefied Richness") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11)); a.rr.p

lanphere.p <- plot_grid(a.ab.p, a.rr.p, aa.arth.ord.GxE, w.ab.p, w.rr.p, w.arth.ord.wind, NULL, w.frr.p, f.ord, NULL, w.brr.p, b.ord.wind, ncol = 3, labels = c("A","B","C","D","E","F","","G","H","","I","J"), align = 'hv'); lanphere.p

save_plot(lanphere.p, filename = "fig_1_test.png", base_height = 11, base_width = 8.5)
