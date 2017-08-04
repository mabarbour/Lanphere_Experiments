## Analysis of root-associated communities

## load required libraries ----
source('scripts_for_analysis/required_libraries.R')

## load require data ----

## Fungal community

f.taxa <- read.csv("final_data/fungi_taxa_table.csv") %>% tbl_df() %>% rename(OTU_ID = X)

fungal.df <- read.csv("final_data/fungal.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
f.OTUs <- colnames(select(fungal.df, -(X:fungal.rarerich)))

## Bacterial community

bacteria.df <- read.csv("final_data/bacteria.df.csv") %>% tbl_df() %>% mutate(Block = as.factor(Block), X = as.factor(X))
b.OTUs <- colnames(select(bacteria.df, -(X:bacteria.rarerich)))


## Fungal rarefied richness analysis ----
# not that I get the same results if I round everything and use Poisson that accounts for overdispersion. Residuals look a bit better with linear model...

hist(fungal.df$fungal.rarerich)
fung.rarerich.lmer <- lmer(fungal.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = fungal.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))
summary(fung.rarerich.lmer)

# no effects on fungal rarerichness
(fungal.rarerich.anova <- anova.table(fung.rarerich.lmer, type = 2, experiment = "wind"))

fungal.rarerich.up <- update(fung.rarerich.lmer, .~. -Wind.Exposure*Genotype -Genotype + Wind.Exposure + (1|Genotype))
(fungal.rarerich.R2 <- var.table(fungal.rarerich.up, experiment = "wind"))


## Bacteria rarefied richness analysis ----
bact.rarerich.lmer <- lmer(bacteria.rarerich ~ Wind.Exposure*Genotype +
                             (1|Block) + (1|Block:Wind.Exposure),
                           data = bacteria.df,
                           contrasts = list(Wind.Exposure = "contr.sum",
                                            Genotype = "contr.sum"))
summary(bact.rarerich.lmer)

# no effects on bacteria rarerichness. Significant effect of wind exposure
anova.table(bact.rarerich.lmer, type = 2, experiment = "wind")

plot(Effect("Wind.Exposure", bact.rarerich.lmer)) 
as.data.frame(Effect("Wind.Exposure", bact.rarerich.lmer))
1153.492/1085.967 # 6% higher rarefied richness on exposed plants.

## Total Fungal Community analysis: hellinger ----

f.hell <- decostand(fungal.df[ ,f.OTUs], method = "hellinger") #f.sub.OTUs[-other.OTUs]], method = "hellinger") # -key.OTUs

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

## identify key taxa driving genotype effect
f.hell.geno.rda <- rda(f.hell ~ Condition(Wind.Exposure) + Condition(Block) + Genotype, data = fungal.df)
summary(f.hell.geno.rda)$cont$importance[ ,1:2] # first 2 axes only explain 2% and 1% of the variance in fungal community composition.

plot(f.hell.geno.rda, display = "sp", type = "text") # visually identify potential key taxa

## test key taxa
fungal.df$fungal.abund <- rowSums(round(fungal.df[ ,f.OTUs]),0)
OTU.lmer <- glmer(round(OTU_128,0) ~ offset(log(fungal.abund)) + Wind.Exposure  + Genotype + (1|X) + (1|Block) + (1|Block:Wind.Exposure), data = fungal.df, family = "poisson", contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
overdisp_fun(OTU.lmer)
summary(OTU.lmer)
drop1(OTU.lmer, test = "Chisq")

# OTU_128, although significant, is not part of the fungal taxa table so I'm excluding it.
taxa.sub <- which(as.character(f.taxa$OTU_ID) %in% c("OTU_54","OTU_68","OTU_134","OTU_66","OTU_119","OTU_262","OTU_1606","OTU_113"))#,"OTU_128"))

table(f.taxa$Guild)
as.data.frame(f.taxa[taxa.sub, ])
write.csv(as.data.frame(f.taxa[taxa.sub, ]), "final_data/key_fungal_taxa.csv")

## Bacteria Community analysis: hellinger distance ----
b.hell <- decostand(bacteria.df[ ,b.OTUs], method = "hellinger")

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
  geom_text(data = sites.cap, aes(x = RDA1, y = RDA2, label = Block), color = "gray") +
  geom_polygon(color = NA, fill = "gray50", alpha = 0.5) + 
  geom_text(data = centroids.cap[11:12, ], 
            aes(x = RDA1, y = RDA2, label = Wind.Exposure), size = 5) +
  ylab("RDA 2 (10%)") + xlab("RDA 1 (14%)") + ggtitle("Bacteria")
b.ord.wind

# create multipanel plot. Run plot scripts for arthropod analysis first though #aa.rGE, w.rGE, f.rich.p, b.rich.p,
(comm.p <- plot_grid(w.arth.ord.wind, f.ord, b.ord.wind, labels = "AUTO", ncol = 3, align = 'h'))

