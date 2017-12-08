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

## ARTHROPOD ABUNDANCE ----
library(brms)
source('scripts_for_analysis/variance_partitioning_brms.R')

# WIND 2012
hist(filter(wind.arth.df, Year == "2012")$total.abund+1)
arth.abund.brm <- brm(total.abund ~ 
                        Wind.Exposure + 
                        (Wind.Exposure|Genotype) +
                        (1|Block) + 
                        (1|Block:Wind.Exposure),
                      data=filter(wind.arth.df, Year == "2012"),
                      family="negbinomial",
                      control=list(adapt_delta=0.99)
)
pp_check(arth.abund.brm, nsamples=100)
tidy(arth.abund.brm)
bayes_R2(arth.abund.brm)
arth.abund.var <- get_variance_wind_brms(arth.abund.brm, data = filter(wind.arth.df, Year == "2012"), RE_rows = 3:7)
arth.abund.var$estimate[6]^2/sum(arth.abund.var$estimate^2)

## Wind: arthropod richness analysis ----
plot(log(total.rich+1) ~ log(total.abund+1), wind.arth.df)

# GLMM
hist(filter(wind.arth.df, Year == "2012")$total.rich)
arth.rich.brm <- brm(total.rich ~ Wind.Exposure +
                           (Wind.Exposure|Genotype) +
                           (1|Block) + 
                           (1|Block:Wind.Exposure),
                     data=filter(wind.arth.df, Year == "2012"),
                     family="poisson",
                     control=list(adapt_delta=0.9))
pp_check(arth.rich.brm, nsamples=100)
tidy(arth.rich.brm)
bayes_R2(arth.rich.brm)
plot(marginal_effects(arth.rich.brm))
pois_var <- log(1 + 1/exp(fixef(arth.rich.brm)["Intercept","Estimate"]))
arth.rich.var <- get_variance_wind_brms(arth.rich.brm, data = filter(wind.arth.df, Year == "2012"), RE_rows = 3:6)
arth.rich.var$estimate[1]^2/(sum(arth.rich.var$estimate^2) + pois_var)



## Wind: rarefied arthropod richness analysis ----
# need to filter data so that there are 2 or more individuals
hist(filter(wind.arth.df, total.abund > 1)$total.rarerich)

# GLMM
arth.rarerich.brm <- brm(total.rarerich ~ Wind.Exposure +
                              (Wind.Exposure|Genotype) +
                              (1|Block) + 
                              (1|Block:Wind.Exposure),
                            data = filter(wind.arth.df, total.abund>1, Year=="2012"),
                         )
print(summary(arth.rarerich.glmer), correlation = TRUE)

## Likelihood ratio tests
arth.rarerich.anova <- anova.table(arth.rarerich.glmer, test = "F", type = 2, experiment = "wind")

## Calculate effects
Effect("Wind.Exposure", arth.rarerich.glmer) # rarefied richness of arthropods was 60% less on exposed plants
Effect("Genotype", arth.rarerich.glmer) # 5.9-fold variation in arthropod rarefied richness among genotypes

w.rrich.GxE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.rarerich.glmer)) %>% mutate(response = "arthropod_rarerich")

# Calculate R2
arth.rarerich.up <- update(arth.rarerich.glmer, .~. -Wind.Exposure*Year -Wind.Exposure*Genotype + Wind.Exposure + (1|Genotype))

(arth.rarerich.R2 <- var.table(arth.rarerich.up, experiment = "wind"))

## Wind: hellinger distance 2012 ----
# no G, wind, or GxE effect
w.12.hell.pos2 <- decostand(w.arth.12.pos2[ ,wind.arth.names], method = "hellinger")

# test G and GxE
w.12.hell.rda <- rda(w.12.hell.pos2 ~ Wind.Exposure*Genotype, data = w.arth.12.pos2) # no G or GxE
anova(w.12.hell.rda, by = "margin", permutations = how(block = w.arth.12.pos2$Block, nperm = 999)) # no GxE
anova(update(w.12.hell.rda, .~. -Wind.Exposure:Genotype), by = "margin", permutations = how(block = w.arth.12.pos2$Block, nperm = 999)) # no G

# test wind effect
w.12.plots.hell <- betadisper(vegdist(w.12.hell.pos2, method = "euclidean"),  w.arth.12.pos2$Plot_code, bias.adjust = TRUE)

w.12.plots.centr.hell <- data.frame(w.12.plots.hell$centroids, 
                                    id = rownames(w.12.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.12.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.12.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.12.plots.centr.hell$Block, nperm = 999)) # no effect of wind exposure

## Wind: Hellinger distance 2013 ----
# effect of wind only, but could be due to overdispersion
key.arths <- which(wind.arth.names %in% c("Tortricidiae_leaftier","Cecidomyiidae_gall"))
w.13.hell.pos2 <- decostand(w.arth.13.pos2[ ,wind.arth.names], method = "hellinger")

# test G and GxE
w.13.hell.rda <- rda(w.13.hell.pos2 ~ Wind.Exposure*Genotype, data = w.arth.13.pos2) 
plot(w.13.hell.rda, display = "sp") # key species are Cecidomyiids and Tortricids

anova(w.13.hell.rda, by = "margin",permutations = how(block = w.arth.13.pos2$Block, nperm = 999)) # no GxE
anova(update(w.13.hell.rda, .~. -Wind.Exposure:Genotype), by = "margin",permutations = how(block = w.arth.13.pos2$Block, nperm = 999)) # no G

# test wind effect
w.13.plots.hell <- betadisper(vegdist(w.13.hell.pos2, method = "euclidean"),  w.arth.13.pos2$Plot_code, bias.adjust = TRUE)

w.13.plots.centr.hell <- data.frame(w.13.plots.hell$centroids, 
                                    id = rownames(w.13.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Wind.Exposure")) %>%
  mutate(Block = as.factor(Block),
         Wind.Exposure = as.factor(Wind.Exposure))

adonis(w.13.plots.hell$centroids ~ Block + Wind.Exposure, 
       data = w.13.plots.centr.hell, 
       method = "euclidean", 
       permutations = how(block = w.13.plots.centr.hell$Block, nperm = 999)) # significant effect of wind exposure

w.hell.13.wind <- betadisper(vegdist(w.13.plots.hell$centroids, method = "euclidean"), w.13.plots.centr.hell$Wind.Exposure, bias.adjust = TRUE)
boxplot(w.hell.13.wind)
plot(w.hell.13.wind)
permutest(w.hell.13.wind, permutations = how(block = w.13.plots.centr.hell$Block, nperm = 999)) # significant effect, suggesting that the effect of wind on centroid location could be due to overdispersion.


w.hell.13.wind.rda <- rda(w.13.hell.pos2 ~ Block*Wind.Exposure, 
                          data = w.arth.13.pos2)
summary(w.hell.13.wind.rda)
plot(w.hell.13.wind.rda, display = c("bp","sp"))
ordiellipse(w.hell.13.wind.rda, groups = interaction(w.arth.13.pos2$Block, w.arth.13.pos2$Wind.Exposure))

w.hell.13.wind.nmds <- metaMDS(vegdist(w.13.plots.hell$centroids, method = "euclidean"))
plot(w.hell.13.wind.nmds)
ordiellipse(w.hell.13.wind.nmds, groups = w.13.plots.centr.hell$Wind.Exposure)
points(w.hell.13.wind.nmds, col = as.numeric(w.13.plots.centr.hell$Wind.Exposure))
text(w.hell.13.wind.nmds, labels = w.13.plots.centr.hell$Block)


## Wind: Caloptilia analysis ----
hist(wind.arth.df$Gracilliaridae_miner)
with(wind.arth.df, table(Gracilliaridae_miner, Year))

# GLMM
LTF.abund.glmer <- glmer(Gracilliaridae_miner ~ Wind.Exposure*Year + 
                           Genotype + # only model without convergence
                           (1|Block) + 
                           (1|Block:Wind.Exposure) +
                           #(1|X) +
                           (1|plant_ID),
                         data = wind.arth.df,
                         contrasts = list(Wind.Exposure = "contr.sum",
                                          Genotype = "contr.sum",
                                          Year = "contr.sum"),
                         family = "poisson",
                         control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))
summary(LTF.abund.glmer)
overdisp_fun(LTF.abund.glmer) 

## Likelihood-ratio test
(LTF.2 <- drop1(LTF.abund.glmer, test = "Chisq")) # G effect
(LTF.2 <- drop1(update(LTF.abund.glmer, .~. -Wind.Exposure:Year), test = "Chisq")) # G and E effect, marginal Year effect

## calculate variance components
LTF.abund.up <- update(LTF.abund.glmer, .~. -Wind.Exposure:Year
                       -Genotype + (1|Genotype)) 
(LTF.abund.R2 <- var.table(LTF.abund.up, experiment = "wind"))



## Wind: Tortricidae leaftier analysis ----
hist(wind.arth.df$Tortricidiae_leaftier)
with(wind.arth.df, table(Tortricidiae_leaftier, Year))

# GLMM
leaftier.abund.glmer <- glmer(Tortricidiae_leaftier ~ 
                                Wind.Exposure*Year + 
                                Wind.Exposure*Genotype + # only model with convergence
                                (1|Block) + 
                                (1|Block:Wind.Exposure) +
                                #(1|X) +
                                (1|plant_ID),
                              data = wind.arth.df, 
                              contrasts = list(Wind.Exposure = "contr.sum",
                                               Genotype = "contr.sum",
                                               Year = "contr.sum"),
                              family = "poisson", 
                              control=glmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))
summary(leaftier.abund.glmer)
overdisp_fun(leaftier.abund.glmer) # no overdispersion

## Likelihood ratio test
(lt.2 <- drop1(leaftier.abund.glmer, test = "Chisq"))
(lt.1 <- drop1(update(leaftier.abund.glmer, .~. -Wind.Exposure:Year -Wind.Exposure:Genotype), test = "Chisq"))
(leaftier.abund.anova <- anova.table(leaftier.abund.glmer, test = "Chisq", experiment = "wind"))

g.tort <- as.data.frame(Effect("Genotype", leaftier.abund.glmer)) %>%
  select(Genotype, tort.fit = fit)

## calculate variance components
leaftier.abund.up <- update(leaftier.abund.glmer, .~. -Wind.Exposure*Year
                            -Wind.Exposure*Genotype + Year + (1|Genotype))
(leaftier.abund.R2 <- var.table(leaftier.abund.up, experiment = "wind"))

## Wind: Cecidomyiidae galler analysis ----
with(wind.arth.df, table(Cecidomyiidae_gall, Year))

# GLMM
gall.abund.glmer <- glmer(Cecidomyiidae_gall ~ Wind.Exposure+Genotype + Year +
                            (1|Block) + 
                            (1|Block:Wind.Exposure) +
                            #(1|X) +
                            (1|plant_ID),
                          weights = total.abund,
                          data = wind.arth.df,
                          contrasts = list(Wind.Exposure = "contr.sum",
                                           Genotype = "contr.sum",
                                           Year = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
summary(gall.abund.glmer)
overdisp_fun(gall.abund.glmer) # no overdispersion

# likelihood ratio test
(gall.1 <- drop1(gall.abund.glmer, test = "Chisq"))

Effect("Wind.Exposure", gall.abund.glmer)
g.gall <- as.data.frame(Effect("Genotype", gall.abund.glmer)) %>%
  select(Genotype, Cecid_gall.fit = fit)

gall.abund.up <- update(gall.abund.glmer, .~. -Genotype + (1|Genotype))

(gall.abund.R2 <- var.table(gall.abund.up, experiment = "wind"))


## Wind: Spider analysis ----
spider.abund.glmer <- glmer(Spider ~ Wind.Exposure*Year + 
                              Genotype +
                              (1|Block) + 
                              (1|Block:Wind.Exposure) +
                              #(1|X) +
                              (1|plant_ID),
                            data = wind.arth.df,
                            contrasts = list(Wind.Exposure = "contr.sum",
                                             Genotype = "contr.sum",
                                             Year = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
summary(spider.abund.glmer)
overdisp_fun(spider.abund.glmer) # no overdispersion
plot(spider.abund.glmer) 

## likelihood ratio tests
(sp.2 <- drop1(spider.abund.glmer, test = "Chisq"))
(sp.1 <- drop1(update(spider.abund.glmer, .~. -Wind.Exposure:Year), test = "Chisq"))

Effect("Wind.Exposure", spider.abund.glmer)

g.spid <- as.data.frame(Effect("Genotype", spider.abund.glmer)) %>%
  select(Genotype, spid.fit = fit)

spider.abund.up <- update(spider.abund.glmer, .~. -Wind.Exposure:Year -Genotype)

(spider.abund.R2 <- var.table(spider.abund.up, experiment = "wind"))

## Wind: Aphididae analysis ----
with(wind.arth.df, table(Aphididae, Year)) # no aphids in 2013, therefore, I only analyze 2012

aphid.abund.glmer <- glmer(Aphididae ~ Wind.Exposure + Genotype +
                             (1|Block) + 
                             (1|Block:Wind.Exposure) +
                             (1|plant_ID), # acts as individual-level random effect since only a single year is being analyzed
                           data = filter(wind.arth.df, Year == "2012"),
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"),
                           family = "poisson",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))
summary(aphid.abund.glmer)
overdisp_fun(aphid.abund.glmer) # sig. overdisper, model with individual-level random effect


# likelihood ratio test
(aphid.1 <- drop1(aphid.abund.glmer, test = "Chisq")) # no wind exposure or genotype effect

g.Aphids_nonAphis <- as.data.frame(Effect("Genotype", aphid.abund.glmer)) %>%
  select(Genotype, Aphids_nonAphis.fit = fit)

## Wind: plots ----
pd <- 0.2
ebar.w <- 0.05
l.size <- 1.5
alp <- 0.5
p.size <- 5
cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")

rich.Gord <- as.data.frame(Effect(c("Genotype"), arth.rich.glmer))  %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)]))

## richness
w.rE <- as.data.frame(Effect(c("Wind.Exposure"), arth.rich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Arthropod richness", limits = c(0,2), breaks = c(0,1,2)) + xlab(""); w.rE

w.rG <- as.data.frame(Effect(c("Genotype"), arth.rich.glmer))  %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0,2), breaks = c(0,1,2)) + xlab(""); w.rG

## abundance
w.aE <- as.data.frame(Effect(c("Wind.Exposure"), arth.abund.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Arthropod abundance", limits = c(0,3.1), breaks = c(0,1,2,3))+ xlab(""); w.aE

w.aG <- as.data.frame(Effect(c("Genotype"), arth.abund.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.Gord$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0, 3.1), breaks = c(0,1,2,3))+ xlab(""); w.aG

## rarefied richness
w.rrE <- as.data.frame(Effect(c("Wind.Exposure"), arth.rarerich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "Rarefied richness", limits = c(0,0.75), breaks = c(0,0.25,0.5,0.75)) + xlab("Wind treatment"); w.rrE

w.rrG <- as.data.frame(Effect(c("Genotype"), arth.rarerich.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.Gord$fit)])) %>% mutate(lower = ifelse(lower < 0, 0, lower)) %>% ggplot(aes(x = Genotype, y = fit)) + geom_point(size = p.size) + geom_errorbar(aes(ymin = lower, ymax = upper), width = ebar.w*2) + scale_y_continuous(name = "", limits = c(0,0.75), breaks = c(0,0.25,0.5,0.75)) + xlab("Willow genotype"); w.rrG # set lower confidence intrval for Genotype G to stop at 0.

w.arth.p <- plot_grid(w.rE, w.rG, w.aE, w.aG, w.rrE, w.rrG, ncol = 2, labels = "AUTO", align = 'hv')

save_plot(w.arth.p, filename = "fig_4_wind_arth_comm.png", base_height = 11, base_width = 8.5)

## composition
cp.keyE <- c( "#999999", "#E69F00","#FFFFFF", "#009E73")
w.key.E <- data.frame(Effect(c("Wind.Exposure"), LTF.abund.glmer), response = "leaf-mining moths", sig = "y") %>% 
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), leaftier.abund.glmer), response = "leaf-tiering moths", sig = "n")) %>%
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), gall.abund.glmer), response = "galling midges", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Wind.Exposure"), spider.abund.glmer), response = "spiders", sig = "y")) %>%
  #bind_rows(., data.frame(Effect(c("Wind.Exposure"), aphid.abund.glmer), response = "aphids")) %>% 
  ggplot(aes(x = Wind.Exposure, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Wind treatment") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyE) + theme(legend.justification = c(1,0), legend.position = c(1,0)); w.key.E #+ geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w, position = position_dodge(width = pd)) # na.rm = TRUE, position = position_dodge(width = pd), 

cp.keyG <- c( "#999999", "#E69F00", "#56B4E9", "#FFFFFF") 
w.key.G <- data.frame(Effect(c("Genotype"), LTF.abund.glmer), response = "leaf-mining moths", sig = "y") %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), leaftier.abund.glmer), response = "leaf-tiering moths", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), gall.abund.glmer), response = "galling midges", sig = "y")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), spider.abund.glmer), response = "spiders", sig = "n")) %>%
  #bind_rows(., data.frame(Effect(c("Genotype"), aphid.abund.glmer), response = "aphids")) %>% 
  ggplot(aes(x = Genotype, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Willow genotype") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyG) + theme(legend.justification = c(1,0), legend.position = c(1,0)); w.key.G

#w.arth.p <- plot_grid(w.rE, w.rG, w.aE, w.aG, w.rrE, w.rrG, w.key.E, w.key.G, labels = "AUTO", ncol = 2, align = 'hv'); w.arth.p

#save_plot(w.arth.p, filename = "fig_wind_arth_test.png", base_height = 11, base_width = 8.5)

## manuscript plots
w.rGE <- as.data.frame(Effect(c("Wind.Exposure","Genotype"), arth.rich.glmer)) %>% ggplot(aes(x = Wind.Exposure, y = fit, color = Genotype)) + geom_line(aes(group = Genotype), size = l.size) + geom_point(size = p.size)  + scale_color_manual(values = cbPal.10) + ylab("No. of species") + xlab("Wind exposure")+ theme(legend.position = "none") + ggtitle("Arthropods")

# Wind effect on community composition
levels(w.13.plots.centr.hell$Wind.Exposure) <- c("Exp.","Unexp.")
w.hell.wind.rda <- rda(w.13.plots.hell$centroids ~ Block + Wind.Exposure, data = w.13.plots.centr.hell) 
summary(w.hell.wind.rda)$cont$importance[ ,1:2] # first axis explains 46% of the variation, second axis explains 7%

ellip <- ordiellipse(w.hell.wind.rda, groups = w.13.plots.centr.hell$Wind.Exposure, draw = "polygon")

centroids.cap <- data.frame(scores(w.hell.wind.rda, display = "cn"), Wind.Exposure = levels(w.13.plots.centr.hell$Wind.Exposure))

sites.cap <- data.frame(scores(w.hell.wind.rda, choices = c(1,2), display = "sites"), droplevels(w.13.plots.centr.hell$Wind.Exposure), droplevels(w.13.plots.centr.hell$Block))
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

w.arth.ord.wind <- ggplot(data = df_ell.cap, aes(x = RDA1, y = RDA2, group = Wind.Exposure)) +
  geom_text(data = sites.cap, aes(x = RDA1, y = RDA2, label = Block), color = "gray") +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap[11:12, ], 
            aes(x = RDA1, y = RDA2, label = Wind.Exposure), size = 5) +
  ylab("RDA 2 (7%)") + xlab("RDA 1 (46%)") + ggtitle("Arthropods")
w.arth.ord.wind

## ant-aphid: maximum Aphis abundance analysis ----

# GLMM
aa.Aphis.glmer <- glmer(aphid_Aphis ~ scale(Ant.mound.dist)*Genotype +
                          (1|plant_ID) + 
                          (1|Block/fact.Ant.mound.dist),
                        data = filter(aa.arth.df, Aphid.treatment == "aphid"),
                        contrasts = list(Genotype = "contr.sum"),
                        family = "poisson",
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e5)))
overdisp_fun(aa.Aphis.glmer)
summary(aa.Aphis.glmer)

## Likelihood ratio tests
(aa.Aphis.2 <- drop1(aa.Aphis.glmer, test = "Chisq"))
(aa.Aphis.2 <- drop1(update(aa.Aphis.glmer, .~. -scale(Ant.mound.dist):Genotype), test = "Chisq"))

## effects
plot(Effect("Ant.mound.dist", aa.Aphis.glmer))
plot(Effect("Genotype", aa.Aphis.glmer)) # genotypes varied 153-fold in mean aphid abundance.
Effect("Genotype", aa.Aphis.glmer)

aa.Aphis.G <- as.data.frame(Effect("Genotype", aa.Aphis.glmer)) %>% mutate(response = "Aphis_farinosa")

## Calculate R2
aa.Aphis.up <- update(aa.Aphis.glmer, .~. -scale(Ant.mound.dist)*Genotype + (1|Genotype))

(aa.Aphis.R2 <- var.table(aa.Aphis.up, experiment = "ant-aphid"))

## ant-aphid: Formica obscuripes abundance analysis ----

# GLMM on abundance
aa.Fobs.abund.glmer <- glmer(ant_F_obscuripes ~ Aphid.treatment + scale(Ant.mound.dist) + (Aphid.treatment|Genotype) +
                               #(1|plant_ID) +
                               (1|Block/fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "poisson",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
summary(aa.Fobs.abund.glmer)
overdisp_fun(aa.Fobs.abund.glmer)


(aa.Fobs.1 <- drop1(aa.Fobs.abund.glmer, test = "Chisq"))
anova(aa.Fobs.abund.glmer, update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype) + (1|Genotype))) # sig GxE
anova(update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype) + (1|Genotype)), update(aa.Fobs.abund.glmer, .~. -(Aphid.treatment|Genotype))) # test G effect

var.table(update(aa.Fobs.abund.glmer, .~. +(1|plant_ID)), "ant-aphid") # need to add individual-level random effect when calculating R2

## ant-aphid: arthropod abundance analysis ----

# GLMM
aa.arth.abund.glmer <- glmer(total.abund ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                               (1|plant_ID) +
                               (1|Block/fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "poisson",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
print(summary(aa.arth.abund.glmer), correlation = TRUE)
overdisp_fun(aa.arth.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.

## Likelihood ratio tests
(aa.abund.3 <- drop1(aa.arth.abund.glmer, test = "Chisq"))
(aa.abund.2 <- drop1(update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aa.abund.1 <- drop1(update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq"))

## calculate effects
plot(Effect("Genotype", aa.arth.abund.glmer)) # arthropod abundance varied 3.8-fold among the most disparate willow genotypes.
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.arth.abund.glmer)) # arthropod abundance increased 1.9-fold in the presence of aphids at increasing distances from ant mounds, but there was little effect of ant-mound distance when there were no aphids.

plot(Effect(c("Aphid.treatment","Genotype"), aa.arth.abund.glmer)) # GxAphid marginal interaction is predominantly influenced by genotypes T and J, which had much higher arthropod abundances in the presence of aphids.

aa.abund.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.arth.abund.glmer)) %>% mutate(response = "arthropod_abund")
aa.abund.ExE <- as.data.frame(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.arth.abund.glmer)) %>% mutate(response = "arthropod_abund")

# Calculate R2
aa.arth.abund.up <- update(aa.arth.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (Aphid.treatment|Genotype))

summary(aa.arth.abund.up)

(aa.arth.abund.R2 <- var.table(aa.arth.abund.up, experiment = "ant-aphid"))

## ant-aphid: arthropod richness analysis ----
plot(total.rich ~ total.abund, aa.arth.df)

# GLMM
aa.arth.rich.glmer <- glmer(total.rich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                              #(1|plant_ID) +
                              (1|Block/fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.arth.rich.glmer), correlation = TRUE)
overdisp_fun(aa.arth.rich.glmer) # no overdispersion

## Likelihood ratio tests
(aa.rich.3 <- drop1(aa.arth.rich.glmer, test = "Chisq"))
(aa.rich.2 <- drop1(update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aa.rich.1 <- drop1(update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -Aphid.treatment:Genotype -scale(Ant.mound.dist):Genotype), test = "Chisq")) # only Genotype effect


## Effects
Effect("Genotype", aa.arth.rich.glmer) # arthropod richness varied 2.7-fold among willow genotypes.

aa.rich.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.arth.rich.glmer)) %>% mutate(response = "arthropod_rich") 

## Calculate R2
aa.arth.rich.up <- update(aa.arth.rich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + (1|Genotype))

(aa.arth.rich.R2 <- var.table(aa.arth.rich.up, experiment = "ant-aphid"))

## ant-aphid: rarefied arthropod richness analysis ----

hist(filter(aa.arth.df, total.abund > 1)$total.rarerich)

# GLMM. Note that logit-transformation (common for proportion data) did not qualitatively affect the outcome or the residuals
aa.arth.rarerich.glmer <- lmer(total.rarerich ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                 #(1|plant_ID) +
                                 (1|Block/fact.Ant.mound.dist),
                               data = filter(aa.arth.df, total.abund > 1),
                               contrasts = list(Aphid.treatment = "contr.sum",
                                                Genotype = "contr.sum"))
summary(aa.arth.rarerich.glmer)

## F tests
(aa.arth.rarerich.anova <- anova.table(aa.arth.rarerich.glmer, test = "F", type = 2, experiment = "ant-aphid")) # marginally significant aphid by genotype interactions

## Effects
plot(Effect(c("Aphid.treatment","Genotype"), aa.arth.rarerich.glmer))
plot(Effect(c("Aphid.treatment"), aa.arth.rarerich.glmer))
1-(0.5816897/0.6917055)

aa.rrich.GxE <- as.data.frame((Effect(c("Aphid.treatment","Genotype"), aa.arth.rarerich.glmer))) %>% mutate(response = "arthropod_rarerich")

## Calculate R2
aa.arth.rarerich.up <- update(aa.arth.rarerich.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist) + Aphid.treatment + (Aphid.treatment|Genotype))

(aa.arth.rarerich.R2 <- var.table(aa.arth.rarerich.up, experiment = "ant-aphid"))

## ant-aphid: Hellinger distance ----
aa.hell <- decostand(aa.arth.12.pos[ ,aa.arth.names], method = "hellinger")

# test 3-way interaction
anova(rda(aa.hell ~ Ant.mound.dist*Aphid.treatment*Genotype, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# test 2-way interactions
anova(rda(aa.hell ~ (Ant.mound.dist + Aphid.treatment + Genotype)^2, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# which GxE are driving the GxEaphid effect?
aa.GxE.filt <- filter(aa.arth.12.pos, Genotype != "J")

aa.GxE.rda <- rda(decostand(aa.GxE.filt[ ,aa.arth.names], "hellinger") ~ Condition(Ant.mound.dist) + Aphid.treatment*Genotype, data = aa.GxE.filt)
plot(aa.GxE.rda)
anova(aa.GxE.rda, by = 'margin', permutations = how(block = aa.GxE.filt$Block, nperm = 999)) # removing Genotype J nullifies GxE effect 
anova(update(aa.GxE.rda, .~. -Aphid.treatment:Genotype), by = 'margin', permutations = how(block = aa.GxE.filt$Block, nperm = 999)) # still main effects of Genotype and aphid treatment though.


# test main effects
anova(rda(aa.hell ~ Ant.mound.dist + Aphid.treatment + Genotype, data = aa.arth.12.pos), by = "margin", permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

# test Ant.mound.dist
aa.plots.hell <- betadisper(vegdist(aa.hell, method = "euclidean"), aa.arth.12.pos$Plot_code, bias.adjust = TRUE)

aa.plots.hell.centr <- data.frame(aa.plots.hell$centroids, 
                                  id = rownames(aa.plots.hell$centroids)) %>%
  separate(col = id, into = c("Block","Ant.mound.dist")) %>%
  mutate(Block = as.factor(Block),
         Ant.mound.dist = as.numeric(Ant.mound.dist))

anova(rda(aa.plots.hell$centroids ~ Block + Ant.mound.dist, data = aa.plots.hell.centr), by = "margin", permutations = how(block = aa.plots.hell.centr$Block, nperm = 999))

## assess assumptions for significant models. No evidence of overdispersion
aa.hell.geno <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$Genotype, bias.adjust = TRUE)
boxplot(aa.hell.geno)
plot(aa.hell.geno)
permutest(aa.hell.geno, permutations = how(block = aa.arth.12.pos$fact.Ant.mound.dist, nperm = 999))

aa.hell.GxE <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$GxE, bias.adjust = TRUE)
boxplot(aa.hell.GxE)
plot(aa.hell.GxE)
permutest(aa.hell.GxE, permutations = how(block = aa.arth.12.pos$Block, nperm = 999))

aa.hell.aphid <- betadisper(vegdist(aa.hell, "euclidean"), aa.arth.12.pos$Aphid.treatment, bias.adjust = TRUE)
boxplot(aa.hell.aphid)
plot(aa.hell.aphid)
permutest(aa.hell.aphid, permutations = how(block = aa.arth.12.pos$fact.Ant.mound.dist, nperm = 999))

## ant-aphid: Identify most abundant taxonomic groups ----
colSums(aa.arth.df[ ,aa.arth.names])
round(colSums(aa.arth.df[ ,aa.arth.names])/sum(colSums(aa.arth.df[ ,aa.arth.names]))*100,0)

## ant-aphid: Caloptilia analysis ----

# GLMM
aa.LTF.abund.glmer <- glmer(LTF_Caloptilia ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                              (1|plant_ID) +
                              (1|Block/fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.LTF.abund.glmer), correlation = TRUE)
overdisp_fun(aa.LTF.abund.glmer) # accounted for overdispersion by modelling individual-level random effect.

## likelihood-ratio tests
(LTF.2 <- drop1(aa.LTF.abund.glmer, test = "Chisq"))
(LTF.1 <- drop1(update(aa.LTF.abund.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))

aa.g.LTF <- as.data.frame(Effect(c("Genotype"), aa.LTF.abund.glmer)) %>%
  select(Genotype, aa.LTF.fit = fit)
plot(Effect(c("Ant.mound.dist"), aa.LTF.abund.glmer))
plot(Effect(c("Ant.mound.dist","Aphid.treatment"), aa.LTF.abund.glmer)) # expected response of LTF to distance and aphids

aa.LTF.abund.up <- update(aa.LTF.abund.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + scale(Ant.mound.dist)*Aphid.treatment + (scale(Ant.mound.dist)|Genotype))

(aa.LTF.abund.R2 <- var.table(aa.LTF.abund.up, experiment = "ant-aphid"))

## ant-aphid: Aphids - non-Aphis analysis ----

# GLMM
aa.Aphids_nonAphis.glmer <- glmer(Aphididae ~ scale(Ant.mound.dist)*Aphid.treatment*Genotype +
                                    (1|plant_ID) +
                                    (1|Block/fact.Ant.mound.dist),
                                  data = aa.arth.df,
                                  contrasts = list(Aphid.treatment = "contr.sum",
                                                   Genotype = "contr.sum"),
                                  family = "poisson",
                                  control=glmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e5)))
summary(aa.Aphids_nonAphis.glmer)
overdisp_fun(aa.Aphids_nonAphis.glmer) # accounted for overdispersion by modelling individual-level random effect.

## likelihood ratio tests
(aphid.3 <- drop1(aa.Aphids_nonAphis.glmer, test = "Chisq"))
(aphid.2 <- drop1(update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype), test = "Chisq"))
(aphid.1 <- drop1(update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment:Genotype -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))

Effect("Genotype", aa.Aphids_nonAphis.glmer)

aa.g.Aphids_nonAphis <- as.data.frame(Effect("Genotype", aa.Aphids_nonAphis.glmer)) %>% # genotypes varied 35-fold in density of non-Aphis aphids.
  select(Genotype, aa.Aphids_nonAphis.fit = fit)
plot(Effect(c("Aphid.treatment","Ant.mound.dist"), aa.Aphids_nonAphis.glmer)) # marginal effect
plot(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer))

aa.Aphids_nonAphis.up <- update(aa.Aphids_nonAphis.glmer, .~. -scale(Ant.mound.dist)*Aphid.treatment*Genotype + Aphid.treatment + (Aphid.treatment|Genotype))
summary(aa.Aphids_nonAphis.up)

(aa.Aphids_nonAphis.R2 <- var.table(aa.Aphids_nonAphis.up, experiment = "ant-aphid"))

## ant-aphid: leafhopper analysis ----
# GLMM
aa.leafhopper.glmer <- glmer(leafhopper ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                               #(1|plant_ID) +
                               (1|Block) +
                               (1|Block:fact.Ant.mound.dist),
                             data = aa.arth.df,
                             contrasts = list(Aphid.treatment = "contr.sum",
                                              Genotype = "contr.sum"),
                             family = "poisson",
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))
print(summary(aa.leafhopper.glmer), correlation = TRUE)
overdisp_fun(aa.leafhopper.glmer) # no overdispersion
plot(aa.leafhopper.glmer) 

# likelihood ratio tests
(lh.2 <- drop1(aa.leafhopper.glmer, test = "Chisq")) 
(lh.1 <- drop1(aa.leafhopper.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype, test = "Chisq")) 


plot(Effect("Genotype",aa.leafhopper.glmer))

aa.leafhopper.up <- update(aa.leafhopper.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + (1|Genotype))

(aa.leafhopper.R2 <- var.table(aa.leafhopper.up, experiment = "ant-aphid"))

## ant-aphid: spiders analysis ----
# GLMM
aa.spiders.glmer <- glmer(spiders ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                            #(1|plant_ID) +
                            (1|Block) +
                            (1|Block:fact.Ant.mound.dist),
                          data = aa.arth.df,
                          contrasts = list(Aphid.treatment = "contr.sum",
                                           Genotype = "contr.sum"),
                          family = "poisson",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))
print(summary(aa.spiders.glmer), correlation = TRUE)
overdisp_fun(aa.spiders.glmer) # no overdispersion

# likelihood ratio tests
(spid.2 <- drop1(aa.spiders.glmer, test = "Chisq"))
(spid.1 <- drop1(update(aa.spiders.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq"))

aa.g.spid <- as.data.frame(Effect("Genotype", aa.spiders.glmer)) %>%
  select(Genotype, aa.spid.fit = fit)

aa.spiders.up <- update(aa.spiders.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + Aphid.treatment + (Aphid.treatment|Genotype))

(aa.spiders.R2 <- var.table(aa.spiders.up, experiment = "ant-aphid"))

## ant-aphid: black ant analysis ----

# GLMM
aa.ant_black.glmer <- glmer(ant_black ~ (scale(Ant.mound.dist) + Aphid.treatment + Genotype)^2 +
                              #(1|plant_ID) +
                              (1|Block) +
                              (1|Block:fact.Ant.mound.dist),
                            data = aa.arth.df,
                            contrasts = list(Aphid.treatment = "contr.sum",
                                             Genotype = "contr.sum"),
                            family = "poisson",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
print(summary(aa.ant_black.glmer), correlation = TRUE)
overdisp_fun(aa.ant_black.glmer) # no evidence of overdispersion
plot(aa.ant_black.glmer) 

# likelihood ratio tests
(ant.2 <- drop1(aa.ant_black.glmer, test = "Chisq")) 
(ant.1 <- drop1(update(aa.ant_black.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment -scale(Ant.mound.dist):Genotype -Aphid.treatment:Genotype), test = "Chisq")) 

plot(Effect("Aphid.treatment",aa.ant_black.glmer))
plot(Effect("Genotype",aa.ant_black.glmer))

aa.ant_black.up <- update(aa.ant_black.glmer, .~. -(scale(Ant.mound.dist)+Aphid.treatment+Genotype)^2 + Aphid.treatment + (1|Genotype) + (1|plant_ID)) # need to add to account for overdispersion

(aa.ant_black.R2 <- var.table(aa.ant_black.up, experiment = "ant-aphid"))

## ant-aphid: leaftier Tortricidae analysis ----

# GLMM
aa.leaftier_Tortricid.glmer <- glmer(leaftier_Tortricid ~  scale(Ant.mound.dist)*Aphid.treatment + Genotype +
                                       #(1|plant_ID) +
                                       (1|Block) +
                                       (1|Block:fact.Ant.mound.dist),
                                     data = aa.arth.df,
                                     contrasts = list(Aphid.treatment = "contr.sum",
                                                      Genotype = "contr.sum"),
                                     family = "poisson",
                                     control=glmerControl(optimizer="bobyqa",
                                                          optCtrl=list(maxfun=2e5)))
summary(aa.leaftier_Tortricid.glmer)
overdisp_fun(aa.leaftier_Tortricid.glmer) # no strong evidence of overdispersion
plot(aa.leaftier_Tortricid.glmer) 

# likelihood ratio tests
(lt.2 <- drop1(aa.leaftier_Tortricid.glmer, test = "Chisq"))
(lt.1 <- drop1(update(aa.leaftier_Tortricid.glmer, .~. -scale(Ant.mound.dist):Aphid.treatment), test = "Chisq"))

plot(Effect(c("Ant.mound.dist"),aa.leaftier_Tortricid.glmer))
aa.g.tort <- as.data.frame(Effect(c("Genotype"),aa.leaftier_Tortricid.glmer)) %>%
  select(Genotype, aa.tort.fit = fit)
plot(Effect(c("Ant.mound.dist","Aphid.treatment"),aa.leaftier_Tortricid.glmer ))

aa.leaftier_Tortricid.up <- update(aa.leaftier_Tortricid.glmer, .~. -Genotype + (1|Genotype) + (1|plant_ID))

(aa.leaftier_Tortricid.R2 <- var.table(aa.leaftier_Tortricid.up, experiment = "ant-aphid"))

##### ant-aphid: manuscript plots ----
pd <- 0.2
#cbPal.10 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66", "#000000")
p.size <- 5
l.size <- 1.5

## richness
rich.sum <- as.data.frame(Effect(c("Genotype"), aa.arth.rich.glmer)) %>%
  mutate(Genotype = factor(Genotype, levels = Genotype[order(fit)])) # re-order Genotype according to increasing richness

aa.r <- rich.sum %>% ggplot(aes(x = Genotype, y = fit)) +  geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + ylab("Arthropod richness") + xlab("Willow genotype")

## rarefied richness
aa.rr <- as.data.frame(Effect(c("Aphid.treatment"), aa.arth.rarerich.glmer)) %>% mutate(Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = fit, shape = Aphid.treatment, fill = Aphid.treatment))  + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd/2)+ geom_point(size = p.size) + ylab("Rarefied richness") + xlab("Aphid treatment") + scale_shape_manual(values = c(23,21), guide = "none") + scale_fill_manual(values = c("gray50", "white"), guide = "none")

## arthropod abundance
aa.aG <- as.data.frame(Effect(c("Genotype"), aa.arth.abund.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.sum$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) +  geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + xlab("Willow genotype") + scale_y_continuous(name = "Arthropod abundance", breaks = c(0,2,4,6,8,10), limits = c(0,10.5))

aa.arth.abund.fact <- glmer(total.abund ~ ord.Ant.mound.dist*Aphid.treatment + Genotype + (1|plant_ID) + (1|Block/fact.Ant.mound.dist), data = aa.arth.df, contrasts = list(Aphid.treatment = "contr.sum", Genotype = "contr.sum"), family = "poisson", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # model ant mound distance as an ordered factor to make plotting easier. Same qualitative result as when I model it as a continuous predictor. Also, I removed interaction with willow genotype since it was not significant and it allowed the model to converge.

aa.aExE <- as.data.frame(Effect(c("Aphid.treatment","ord.Ant.mound.dist"), aa.arth.abund.fact)) %>% mutate(ord.Ant.mound.dist = ordered(as.numeric(as.character(ord.Ant.mound.dist))), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = ord.Ant.mound.dist, y = fit, fill = Aphid.treatment, shape = Aphid.treatment)) + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd, position = position_dodge(width = 0.5)) + geom_point(size = p.size, position = position_dodge(width = 0.5)) + ylab("No. of individuals") + xlab("Distance from ant mound (m)") + scale_shape_manual(name = "Aphid treatment", values = c(23,21)) + scale_fill_manual(name = "Aphid treatment", values = c("gray50", "white")) + scale_y_continuous(name = "Arthropod abundance", breaks = c(0,2,4,6,8,10), limits = c(0,10.5)) + theme(legend.justification = c(0,1), legend.position = c(0,1)); aa.aExE

## ordination of community compositino
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
    geom_point(shape = 1, color = "gray") +
    geom_segment(data = aa.hell.GxE.df.arrows, aes(x = x, xend = xend, y = y, yend = yend, linetype = GxE.sig.lab), color = "black") + #, color = Genotype
    geom_point(data = aa.hell.GxE.df, aes(x = RDA1, y = RDA2, fill = Aphid.treatment, shape = Aphid.treatment), color = "black", size = 5) + 
    geom_text(data = aa.hell.GxE.df, 
              aes(x = RDA1, y = RDA2, label = Genotype), size = 3) +
    scale_fill_manual(values = c("gray50","white"))+ #cbPal.10) +
    scale_color_manual(values = c("gray50","white")) +# cbPal.10) +
    scale_shape_manual(values = c(23,21)) +
    scale_linetype_manual(values = c("dotted","solid"), guide = "none") +
    ylab("RDA 2 (3%)") + xlab("RDA 1 (9%)") + theme(legend.position = "none"))

## visualize non-Aphis aphids driving GxE and which genotype
aa.Aphid.GxE <- as.data.frame(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer)) %>% mutate(sig.GxE = ifelse(Genotype == "J", "yes","no"), Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = fit, group = Genotype, fill = Aphid.treatment, shape = Aphid.treatment)) + geom_line(aes(linetype = sig.GxE), size = 1) + geom_point(size = p.size) + xlab("Aphid treatment") + scale_fill_manual(values = c("gray50","white"), guide = "none") + scale_shape_manual(values = c(23,21), guide = "none") + geom_text(aes(label = Genotype), size = 3) + scale_linetype_manual(values = c("dotted","solid"), guid = "none") + scale_y_continuous(name = "Abundance of other aphids", breaks = c(0,1,2), limits = c(0,2)); aa.Aphid.GxE

## establish multi-panel plot
aa.arth.p <- plot_grid(aa.r, aa.aG, aa.aExE, aa.rr, aa.arth.ord.GxE, aa.Aphid.GxE, labels = "AUTO", ncol = 3, align = 'hv'); aa.arth.p

save_plot(aa.arth.p, filename = "fig_1_aa_arth_comm.png", base_height = 8.5, base_width = 11)

##### ant-aphid: Aphis farinosa and Formica obscuripes responses ----

## Aphis farinosa
aa.Aphis.p <- as.data.frame(Effect("Genotype", aa.Aphis.glmer)) %>% mutate(Genotype = factor(Genotype, levels = Genotype[order(rich.sum$fit)])) %>% ggplot(aes(x = Genotype, y = fit)) + geom_errorbar(aes(ymax = upper, ymin = lower), width = pd) + geom_point(size = p.size, shape = 22, fill = "black") + ylab(expression(paste(italic("Aphis farinosa")," abundance"))) + xlab("Willow genotype") + scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1,10,20, 30), limits = c(0.001,30), labels = c(0.001, 0.01, 0.1, 1,10,20, 30)); aa.Aphis.p

aa.Fobs.p <- aa.arth.df %>% group_by(Genotype, Aphid.treatment) %>% summarise(mean.Fobs = mean(ant_F_obscuripes))  %>% ungroup() %>% mutate(Aphid.treatment = ifelse(Aphid.treatment == "aphid","Aphids","None")) %>% ggplot(aes(x = Aphid.treatment, y = mean.Fobs, group = Genotype, shape = Aphid.treatment, fill = Aphid.treatment)) + geom_line() + geom_point(size = p.size) + geom_text(aes(label = Genotype), size = 3) + scale_shape_manual(values = c(23,21), guide = "none") + scale_fill_manual(values = c("gray50", "white"), guide = "none") + xlab("Aphid treatment") + ylab(expression(paste(italic("Formica obscuripes")," abundance"))) + scale_y_continuous(breaks = c(0,0.2, 0.4, 0.6), limits = c(0,0.6)); aa.Fobs.p #+ scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1), limits = c(0,1))

## ant-aphid: plots for supp mat -----
# consider turning one significant genotype to black solid line and circle, while everything else is grey with white circle and different line types.
cp.keyEaa <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#009E73", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")
aa.key.E <- data.frame(Effect(c("Aphid.treatment","Genotype"), aa.Aphids_nonAphis.glmer), response = "aphids") %>% 
  ggplot(aes(x = Aphid.treatment, y = fit, group = Genotype, color = Genotype, fill = Genotype)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Genotype", values = cbPal.10) + ylab("No. of aphid individuals") + xlab("Aphid treatment")  + scale_fill_manual(name = "Genotype", values = cp.keyEaa) + scale_y_log10() + theme(legend.justification = c(1,0), legend.position = c(1,0)) #+ geom_errorbar(aes(ymax = upper, ymin = lower), width = ebar.w, position = position_dodge(width = pd)) # na.rm = TRUE, position = position_dodge(width = pd), 

cp.keyGaa <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#FFFFFF") # consider altering colors to match other wind graphs
aa.key.G <- data.frame(Effect(c("Genotype"), aa.LTF.abund.glmer), response = "leaf-mining moths") %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), aa.leaftier_Tortricid.glmer), response = "leaf-tiering moths")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), aa.ant_black.glmer), response = "ants")) %>%
  bind_rows(., data.frame(Effect(c("Genotype"), aa.spiders.glmer), response = "spiders")) %>%
  #bind_rows(., data.frame(Effect(c("Genotype"), aa.Aphids_nonAphis.glmer), response = "aphids")) %>% 
  bind_rows(., data.frame(Effect(c("Genotype"), aa.leafhopper.glmer), response = "leafhoppers")) %>% 
  ggplot(aes(x = Genotype, y = fit, group = response, color = response, fill = response)) + geom_line(size = l.size) + geom_point(size = p.size, shape = 21) + scale_color_manual(name = "Arthropod guild", values = cbPal.10) + ylab("No. of individuals") + xlab("Willow genotype") + scale_y_log10() + scale_fill_manual(name = "Arthropod guild", values = cp.keyGaa) + theme(legend.justification = c(1,0), legend.position = c(1,0))

