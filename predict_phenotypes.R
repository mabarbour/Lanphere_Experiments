## Predicting individual and community phenotypes

## load libraries ----
source('~/Documents/miscellaneous_R/model_diagnostic_functions.R') # need to load before dplyr
library(rdryad)
library(dplyr)
library(tidyr)

## upload and manage data ----
# download data from dryad
#cg.url <- download_url('10255/dryad.77493')
#dryad_fetch(cg.url, 'herb.trait.data.DRYAD.csv')

# upload data set
cg.dat <- read.csv('herb.trait.data.DRYAD.csv') 

# manage dataset for analysis
cg.df <- cg.dat %>%
  dplyr::select(X, Genotype, Height, water_content, specific_leaf_area, C_N_ratio, Trichome.No., Herbivore_abundance, vLG, SG, PG, LTF) %>%
  mutate(X = as.factor(X),
         Cecidomyiid_gall = vLG + SG + PG,
         Gracilliaridae_miner = LTF) %>%
  dplyr::select(X:Herbivore_abundance, Cecidomyiid_gall, Gracilliaridae_miner)

# calculate phenotype means
cg.df.mean <- cg.df %>%
  group_by(Genotype) %>%
  summarise_each(funs(mean.narm = mean(., na.rm = TRUE)))

# upload data from experiment
trait.H2s <- read.csv("trait.R2s.csv") %>%
  tbl_df() %>%
  filter(Factors %in% c("genotype.(Intercept)", "Genotype.(Intercept)", "genotype.Year2012","genotype.Year2013"), response %in% c("Height","SLA","leaf_WC","leaf_trichome.density","leaf_C_N")) %>%
  dplyr::select(Response = response, H2 = var_percent, experiment) %>%
  group_by(experiment, Response) %>%
  summarise_each(funs(mean)) %>% # calculating average heritability for leaf water content across both years.
  ungroup() 

# work on extracting adjusted H2s
trait.adj.H2s <- read.csv("trait.R2s.csv") %>%
  tbl_df() %>%
  filter(Factors %in% c("genotype.(Intercept)", "Genotype.(Intercept)", "genotype.Year2012","genotype.Year2013", "var_Resid"), response %in% c("Height","SLA","leaf_WC","leaf_trichome.density","leaf_C_N")) %>%
  dplyr::select(Response = response, Factors, var, experiment) %>%
  spread(Factors, var)# %>%
trait.adj.H2s$totGenVar <- rowSums(trait.adj.H2s[ ,3:6], na.rm = TRUE)

trait.adj.H2s <- trait.adj.H2s %>%
  mutate(H2 = totGenVar/(totGenVar + var_Resid)) %>%
  dplyr::select(Response, experiment, H2) 

## Calculate broad-sense heritabilities ----
# consider restricting data to genotypes from the experiment...

# plant height 
cg.height <- lmer(Height ~ (1|Genotype), data = cg.df)
cg.height.h <- var.table(cg.height, experiment = "cg")[2,"var_percent"]

# leaf water content
cg.WC <- lmer(water_content ~ (1|Genotype), data = cg.df)
cg.WC.h <- var.table(cg.WC, experiment = "cg")[2,"var_percent"]

# SLA
cg.SLA <- lmer(specific_leaf_area ~ (1|Genotype), data = cg.df)
cg.SLA.h <- var.table(cg.SLA, experiment = "cg")[2,"var_percent"]

# leaf C:N
cg.CN <- lmer(C_N_ratio ~ (1|Genotype), data = cg.df)
cg.CN.h <- var.table(cg.CN, experiment = "cg")[2,"var_percent"] 

# Trichome density
cg.tri <- lmer(log(Trichome.No.+1) ~ (1|Genotype), data = cg.df)
cg.tri.h <- var.table(cg.tri, experiment = "cg")[2,"var_percent"]

# Herbivore abundance
cg.herb <- glmer(Herbivore_abundance ~ (1|Genotype) + (1|X), data = cg.df, family = "poisson")
(cg.herb.h <- var.table(cg.herb, experiment = "cg")[3,"var_percent"])

# Cecidomyiid galls
cg.gall <- glmer(Cecidomyiid_gall ~ (1|Genotype) + (1|X), data = cg.df, family = "poisson")
(cg.gall.h <- var.table(cg.gall, experiment = "cg")[3,"var_percent"])

# Gracilliaridae miners
cg.mine <- glmer(Gracilliaridae_miner ~ (1|Genotype) + (1|X), data = cg.df, family = "poisson")
(cg.mine.h <- var.table(cg.mine, experiment = "cg")[3,"var_percent"])

# make data frame of heritabilities
cg.herit.df <- data.frame(Response = c("Height","SLA","leaf_WC", "leaf_C_N","leaf_trichome.density","Herb_abund","Cecid_gall","Grac_mine"),
                          bind_rows(cg.height.h, cg.SLA.h, cg.WC.h, cg.CN.h, cg.tri.h, cg.herb.h, cg.gall.h, cg.mine.h),
                          experiment = "common_garden") %>% # for comparison, I feel they are adjusted since they are contingent on a single common garden
  dplyr::select(Response, H2 = var_percent, experiment)

# join heritability data frames
herit.df <- bind_rows(cg.herit.df, trait.H2s) %>%
  spread(experiment, H2)

adj.herit.df <- bind_rows(cg.herit.df, trait.adj.H2s) %>%
  spread(experiment, H2)

plot(herit.df$'ant-aphid' ~ common_garden, herit.df)
plot(herit.df$wind ~ common_garden, herit.df)

plot(adj.herit.df$'ant-aphid' ~ common_garden, adj.herit.df)
plot(adj.herit.df$wind ~ common_garden, adj.herit.df)

## Calculate genetic correlations ----
# Need lanphere_plant_trait_analysis and lanphere_arthropod_analysis scripts run first.
gen.corr.df <- left_join(cg.df.mean, g.height) %>%
  left_join(., g.tri) %>%
  left_join(., g.CN) %>%
  left_join(., g.SLA) %>%
  left_join(., g.WC) %>%
  left_join(., aa.g.tri) %>%
  left_join(., aa.g.height) %>%
  left_join(., aa.g.WC) %>%
  left_join(., g.LTF) %>%
  left_join(., g.gall) %>%
  left_join(., aa.g.LTF)

# wind: height. Genotype G is an outlier
plot(height.fit ~ Height, gen.corr.df, type = "n")
text(y = gen.corr.df$height.fit, x = gen.corr.df$Height, labels = gen.corr.df$Genotype)
(height.cor <- with(gen.corr.df, cor.test(height.fit, Height, method = "pearson")))

# ant-aphid: height
plot(aa.height.fit ~ Height, gen.corr.df, type = "n")
text(y = gen.corr.df$aa.height.fit, x = gen.corr.df$Height, labels = gen.corr.df$Genotype)
(aa.height.cor <- with(gen.corr.df, cor.test(aa.height.fit, Height, method = "pearson")))

plot(height.fit ~ aa.height.fit, gen.corr.df)
with(gen.corr.df, cor.test(height.fit, aa.height.fit, method = "pearson"))
  
# wind: trichome density
plot(trichome.fit ~ Trichome.No., gen.corr.df)
(trich.cor <- with(gen.corr.df, cor.test(trichome.fit, Trichome.No., method = "pearson")))

# ant-aphid: trichome density
plot(aa.trichome.fit ~ Trichome.No., gen.corr.df)
(aa.trich.cor <- with(gen.corr.df, cor.test(aa.trichome.fit, Trichome.No., method = "pearson")))

plot(trichome.fit ~ aa.trichome.fit, gen.corr.df) # look at correlation between wind and ant-aphid.
with(gen.corr.df, cor.test(trichome.fit, aa.trichome.fit))

# leaf C:N
plot(CN.fit ~ C_N_ratio, gen.corr.df)
(CN.cor <- with(gen.corr.df, cor.test(CN.fit, C_N_ratio, method = "pearson")))

# SLA
plot(SLA.fit ~ specific_leaf_area, gen.corr.df, type = "n")
text(y = gen.corr.df$SLA.fit, x = gen.corr.df$specific_leaf_area, labels = gen.corr.df$Genotype)
(SLA.cor <- with(gen.corr.df, cor.test(SLA.fit, specific_leaf_area, method = "pearson")))

# wind: leaf water content
plot(WC.fit ~ water_content, gen.corr.df)
text(y = gen.corr.df$WC.fit, x = gen.corr.df$water_content, labels = gen.corr.df$Genotype)
(WC.cor <- with(gen.corr.df, cor.test(WC.fit, water_content, method = "pearson")))

# ant-aphid: leaf water content
plot(aa.WC.fit ~ water_content, gen.corr.df)
(aa.WC.cor <- with(gen.corr.df, cor.test(aa.WC.fit, water_content, method = "pearson")))

plot(WC.fit ~ aa.WC.fit, gen.corr.df)
with(gen.corr.df, cor.test(aa.WC.fit, WC.fit, method = "pearson"))

# wind: Caloptilia
plot(LTF.fit ~ Gracilliaridae_miner, gen.corr.df)
(LTF.cor <- with(gen.corr.df, cor.test(LTF.fit, Gracilliaridae_miner, method = "pearson")))

plot(Gracilliaridae_miner ~ Height, gen.corr.df)

# ant-aphid: Caloptilia
plot(aa.LTF.fit ~ Gracilliaridae_miner, gen.corr.df)
(aa.LTF.cor <- with(gen.corr.df, cor.test(aa.LTF.fit, Gracilliaridae_miner, method = "pearson")))

plot(LTF.fit ~ aa.LTF.fit, gen.corr.df)
with(gen.corr.df, cor.test(LTF.fit, aa.LTF.fit))

# wind: gall midges
plot(Cecid_gall.fit ~ Cecidomyiid_gall, gen.corr.df)
(gall.cor <- with(gen.corr.df, cor.test(Cecid_gall.fit, Cecidomyiid_gall, method = "pearson")))

plot(Cecid_gall.fit ~ height.fit, gen.corr.df)
plot(Cecidomyiid_gall ~ Height, gen.corr.df) # negative correlation with Height???

cg_wind.cors <- data.frame(cor.cg_wind = c(gall.cor$estimate,
                                   LTF.cor$estimate,
                                   height.cor$estimate,
                                   trich.cor$estimate,
                                   WC.cor$estimate,
                                   SLA.cor$estimate,
                                   CN.cor$estimate),
                           Response = c("Cecid_gall","Grac_mine","Height","leaf_trichome.density","leaf_WC","SLA","leaf_C_N"))

cg_aa.cors <- data.frame(cor.cg_aa = c(aa.LTF.cor$estimate,
                                   aa.height.cor$estimate,
                                   aa.trich.cor$estimate,
                                   aa.WC.cor$estimate),
                           Response = c("Grac_mine","Height","leaf_trichome.density","leaf_WC"))

H2.cors <- left_join(herit.df, cg_wind.cors) %>%
  left_join(., cg_aa.cors)

plot(cor.cg_wind ~ common_garden, H2.cors)
summary(lm(cor.cg_wind ~ common_garden, H2.cors))

plot(cor.cg_aa ~ common_garden, H2.cors)
summary(lm(cor.cg_aa ~ common_garden, H2.cors))
