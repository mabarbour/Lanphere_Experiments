
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)

## GET DATA ----
sd.df <- read_csv('output_brms/lanphere_SDs.csv') %>% mutate(term_group="") %>% as.data.frame()

## FILL IN TERM GROUPS AND ORDER THEM
sd.df$term_group[which(sd.df$term=="sd_Genotype__Intercept")] <- "Genotype (G)"
sd.df$term_group[which(sd.df$term=="sd_Wind.Exposure")] <- "Wind"
sd.df$term_group[which(sd.df$term=="sd_Aphid.treatment")] <- "Aphid"
sd.df$term_group[which(sd.df$term=="sd_c.Ant.mound.dist")] <- "Ant"
sd.df$term_group[which(sd.df$term=="sd_Aphid.x.Ant")] <- "Aphid x Ant"
sd.df$term_group[which(sd.df$term=="sd_Block__Intercept")] <- "Block"
sd.df$term_group[which(sd.df$term=="sd_Plot_code__Intercept")] <- "Plot"
sd.df$term_group[which(sd.df$term=="sd_plant_ID__Intercept")] <- "Individual"
sd.df$term_group[which(sd.df$term=="sd_Genotype__Wind.Exposure1")] <- "G x Wind"
sd.df$term_group[which(sd.df$term=="sd_Genotype__Aphid.treatment1")] <- "G x Aphid"
sd.df$term_group[which(sd.df$term=="sd_Genotype__c.Ant.mound.dist")] <- "G x Ant"
sd.df$term_group[which(sd.df$term=="sd_Genotype__Aphid.treatment1.c.Ant.mound.dist")] <- "G x Aphid x Ant"

sd.df$term_group <- factor(sd.df$term_group)
levels(sd.df$term_group)
sd.df$term_group_ord <- factor(sd.df$term_group, levels=c("Plot","Block","G x Aphid x Ant","G x Ant","G x Aphid","G x Wind","Aphid x Ant","Ant","Aphid","Wind","Genotype (G)"))
levels(sd.df$term_group_ord)


## FIGURE 1: COMMUNITY RICHNESS ----
plot_richness <- filter(sd.df, Response%in%c("Arthropod Richness", "scale(Fungi Rarefied Richness)", "scale(Bacteria Rarefied Richness)")) %>% 
  droplevels() %>% 
  mutate(Experiment_Year=paste(Experiment, Year, " "))

richness_gg <- ggplot(plot_richness, aes(x=term_group_ord, y=posterior_SD, fill=Response)) +
  geom_boxplot(outlier.shape = 1) +
  coord_flip() +
  ylab("Effect size (SD)") +
  xlab("") +
  facet_wrap(Year~Experiment, ncol=1, scales="free_y") 
richness_gg


## SUPP MAT 1 ----

# Color-blind friendly palette:
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_traits <- filter(varcomp.df, Response%in%c("Trait.PC1", "Trait.PC2", "Root C:N"), term != "sd_plant_ID__Intercept") %>% 
  droplevels() %>% 
  mutate(Experiment_Year=paste(Experiment, Year, " "))

traits_gg <- ggplot(plot_traits, aes(x=term_group_ord, y=SD_mean, fill=Response)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95, color=Response), size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50, color=Response), size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", shape=21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Effect size (SD)") +
  scale_fill_manual(values=cbPalette[3:5], name="Trait response") +
  scale_color_manual(values=cbPalette[3:5], name="Trait response") +
  xlab("") +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y") 

plot_soil <- filter(varcomp.df, Response%in%c("Soil.PC1", "Soil.PC2"), term != "sd_plant_ID__Intercept") %>% 
  droplevels() %>% 
  mutate(Experiment_Year=paste(Experiment, Year, " "))

soil_gg <- ggplot(plot_soil, aes(x=term_group_ord, y=SD_mean, fill=Response)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95, color=Response), size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50, color=Response), size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", shape=21, position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Effect size (SD)") +
  scale_fill_manual(values=cbPalette[8:9], name="Soil response") + 
  scale_color_manual(values=cbPalette[8:9], name="Soil response") +
  xlab("") +
  facet_grid(~Experiment_Year, scales="free_y")

plot_grid(traits_gg, soil_gg, ncol=1, align = 'v', labels = "AUTO", rel_heights=c(1,0.25))

## TRAIT - COMMUNITY RELATIONSHIPS ----

## NEED RELEVANT TRAIT/PC DATA FOR VARIANCE EXPLAINED AND LOADINGS
## NEED REGRESSIONS FROM BRM ANALYSIS FOR REGRESSION PLOTS

plot_VarExp <- function(Resp, Exper, Yr){
  response_VarExp <- filter(varcomp.df, Response%in%Resp, Experiment%in%Exper, Year%in%Yr)
  response_VarExp$term_group_ord <- factor(response_VarExp$term_group_ord, levels=rev(levels(response_VarExp$term_group_ord)))
  
  ggplot(response_VarExp, aes(x=Response, y=VarComp_mean*100, fill=term_group_ord)) +
    geom_bar(stat = "identity", position=position_stack(reverse=T), color=NA, show.legend = F) + 
    coord_flip() +
    geom_text(aes(label=term_group_ord), stat = "identity", 
              position=position_stack(reverse=T, vjust=0.5),
              angle=45, show.legend=F, size=3, check_overlap=T) + 
    scale_fill_grey(start=0.5, end=0.9)# +
  #xlab("") + 
  #scale_x_continuous(labels=NULL, breaks=NULL)
}

## ANT-APHID 2012 -- ARTHROPOD RICHNESS & TRAIT PC1 ----
VarExp_aa.trait.PC1.12 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Ant-Aphid", Yr = "2012") + ylab("Trait PC1 explained variance (%)")+
  theme_cowplot(font_size=10)

Arrow_aa.trait.PC1.12 <- aa.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC1, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(-0.1,-0.1,-0.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC1_aa.12 <- plot(marginal_effects(trait.rich.aa.2012.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                                points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_aa.trait.PC1.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_aa.trait.PC1.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

aa.traitreg.pc1 <- plot_grid(ArthRich.TraitPC1_aa.12, VarExp_aa.trait.PC1.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 

## ANT-APHID 2012 -- ARTHROPOD RICHNESS & TRAIT PC2 ----
VarExp_aa.trait.PC2.12 <- plot_VarExp(Resp = "Trait.PC2", Exper = "Ant-Aphid", Yr = "2012") + ylab("Trait PC2 explained variance (%)")+
  theme_cowplot(font_size=10)

Arrow_aa.trait.PC2.12 <- aa.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC2, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(1.1,-0.1,1.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC2_aa.12 <- plot(marginal_effects(trait.rich.aa.2012.brm, effects="trait.PC2", probs=c(0.025,0.975)), 
                                points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC2") +
  geom_segment(data=Arrow_aa.trait.PC2.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_aa.trait.PC2.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

aa.traitreg.pc2 <- plot_grid(ArthRich.TraitPC2_aa.12, VarExp_aa.trait.PC2.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 


## WIND 2012 -- ARTRHOPOD RICHNESS & TRAIT PC1 ----
VarExp_w.trait.PC1.12 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Wind", Yr = "2012") + ylab("Trait PC1 explained variance (%)") +
  theme_cowplot(font_size=10)

Arrow_w.trait.PC1.12 <- wind.trait.PCA.2012 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), 
         xend=PC1, yend=seq(-1,-0.25,length.out=5), 
         hjust=c(-0.1,-0.1,-0.1,1.1,-0.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

ArthRich.TraitPC1_w.12 <- plot(marginal_effects(trait.rich.wind.2012.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                               points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_w.trait.PC1.12, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.trait.PC1.12, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) +
  theme_cowplot(font_size=10)

w.traitreg.pc1 <- plot_grid(ArthRich.TraitPC1_w.12, VarExp_w.trait.PC1.12, ncol=1, align = 'v', scale=c(1,0.6)) #rel_heights = c(1,0.4), 


## WIND 2013 -- ARTRHOPOD RICHNESS & TRAIT PC1 ----
VarExp_w.trait.PC1.13 <- plot_VarExp(Resp = "Trait.PC1", Exper = "Wind", Yr = "2013")+ ylab("Trait PC1 explained variance (%)")

Arrow_w.trait.PC1.13 <- wind.trait.PCA.2013 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=6), 
         xend=PC1, yend=seq(-1,-0.25,length.out=6), 
         hjust=c(-0.1,-0.1,-0.1,1.1,1.1,1.1),
         label=c("Height","Shoot count", "Shoot length", "SLA", "Water content", "C:N"))

ArthRich.TraitPC1_w.13 <- plot(marginal_effects(trait.rich.wind.2013.brm, effects="trait.PC1", probs=c(0.025,0.975)), 
                               points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=Arrow_w.trait.PC1.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.trait.PC1.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

w.traitreg.pc1.2013 <- plot_grid(ArthRich.TraitPC1_w.13, VarExp_w.trait.PC1.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- BACTERIA RAREFIED RICHNESS & SOIL PC2 ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.soil.PC2.13 <- plot_VarExp(Resp = "Soil.PC2", Exper = "Wind", Yr = "2013") + ylab("Soil PC2 explained variance (%)")

Arrow_w.soil.PC2.13 <- wind.soil.PCA %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=17), 
         xend=PC2, yend=seq(-1,-0.25,length.out=17), 
         hjust=c(1.1,1.1,-0.1,-0.1,1.1,-0.1,1.1,-0.1,1.1,1.1,1.1,-0.1,1.1,1.1,-0.1,-0.1,-0.1),
         label=c("NO_3","NH_4", "Ca","Mg","K","P","Fe","Mn","Cu","Zn","B","S","Pb","Al","Cd","% TOM","Moisture"))

BactRareRich.SoilPC2_w.13 <- plot(marginal_effects(bact.trait.rarerich.wind.2013.brm, effects="soil.PC2", probs=c(0.025,0.975)), 
                                  points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Bacteria rarefied richness") + xlab("Soil PC2") +
  geom_segment(data=Arrow_w.soil.PC2.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.soil.PC2.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

bact.reg.pc2 <- plot_grid(BactRareRich.SoilPC2_w.13, VarExp_w.soil.PC2.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- FUNGI RAREFIED RICHNESS & SOIL PC1 ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.soil.PC1.13 <- plot_VarExp(Resp = "Soil.PC1", Exper = "Wind", Yr = "2013")+ ylab("Soil PC1 explained variance (%)")

Arrow_w.soil.PC1.13 <- wind.soil.PCA %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=17), 
         xend=PC1, yend=seq(-1,-0.25,length.out=17), 
         hjust=c(1.1,1.1,rep(-0.1,15)),
         label=c("NO_3","NH_4", "Ca","Mg","K","P","Fe","Mn","Cu","Zn","B","S","Pb","Al","Cd","% TOM","Moisture"))

FungRareRich.SoilPC2_w.13 <- plot(marginal_effects(trait.rarerich.wind.2013.brm, effects="soil.PC1", probs=c(0.025,0.975)), 
                                  points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Fungi rarefied richness") + xlab("Soil PC1") +
  geom_segment(data=Arrow_w.soil.PC1.13, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=Arrow_w.soil.PC1.13, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) 

fung.reg.pc1 <- plot_grid(FungRareRich.SoilPC2_w.13, VarExp_w.soil.PC1.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))


## WIND 2013 -- FUNGI RAREFIED RICHNESS & ROOT C:N ----

## NEED TO SCALE BACK RAREFIED RICHNESS

VarExp_w.root.CN.13 <- plot_VarExp(Resp = "Root C:N", Exper = "Wind", Yr = "2013") + ylab("Root C:N explained variance (%)")

FungRareRich.RootCN_w.13 <- plot(marginal_effects(trait.rarerich.wind.2013.brm, effects="log_root_CN", probs=c(0.025,0.975)), 
                                 points=T, point_args=list(color="grey", shape=1), line_args=list(color="black"))[[1]] + 
  ylab("Fungi rarefied richness") + xlab("log(Root C:N)") 

fung.reg.rootCN <- plot_grid(FungRareRich.RootCN_w.13, VarExp_w.root.CN.13, ncol=1, align = 'v', rel_heights = c(1,0.6))#, scale=c(1,0.6))



## MULTIPANEL PLOTS ----

aa.traits.var.2012 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2"), Exper = c("Ant-Aphid"), Yr = "2012")+xlab("") + ylab("Variance explained (%)")
w.traits.var.2012 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2"), Exper = c("Wind"), Yr = "2012")+xlab("") + ylab("Variance explained (%)")

## NEED TO ADD WIND TRAITS INTO PLOT, PROBABLY CAN'T DO WITH CURRENT PLOT_VAREXP FUNCTION. MAYBE JUST SELECT THE DATA FRAME AND THEN FEED INTO THE PLOTTING FUNCTION...OR GET RID OF THE FUNCTION BECAUSE I WILL ONLY NEED IT 2X AND NOT 5X LIKE I ORIGINALLY THOUGHT...
test <- plot_grid(ArthRich.TraitPC1_aa.12, ArthRich.TraitPC2_aa.12, 
                  ArthRich.TraitPC1_w.12, aa.traits.var.2012,
                  #VarExp_aa.trait.PC1.12, VarExp_aa.trait.PC2.12, VarExp_w.trait.PC1.12, 
                  align='hv')
save_plot("test.pdf", test, base_height = 8.5, base_width = 8.5)

traits.var.2013 <- plot_VarExp(Resp=c("Trait.PC1", "Trait.PC2", "Root C:N"), Exper = "Wind", Yr = "2013")+xlab("") + ylab("Variance explained (%)") #+ xlab("Plant traits") 
soil.var <- plot_VarExp(Resp=c("Soil.PC1", "Soil.PC2"), Exper = "Wind", Yr = "2013")+xlab("")  + ylab("Variance explained (%)") #+ xlab("Soil properties")


test.wind.2013 <- plot_grid(ArthRich.TraitPC1_w.13, BactRareRich.SoilPC2_w.13, 
                            FungRareRich.RootCN_w.13, FungRareRich.SoilPC2_w.13, 
                            traits.var.2013, soil.var,
                            ncol=2, align = 'hv', labels = "AUTO")
save_plot("test.wind.2013.pdf", test.wind.2013, base_height = 11, base_width = 8.5)
#plot_grid(w.traitreg.pc1.2013, fung.reg.rootCN, bact.reg.pc2, fung.reg.pc1, ncol=2, align = 'hv', labels = "AUTO")

