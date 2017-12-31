
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)

## GET DATA ----
# testing for composition
#varcomp.df <- read_csv('output_brms/lanphere_arthropods_VarComps_2012.csv') %>% mutate(term_group="", term_order="") %>% as.data.frame() %>% rename(Response=Trait)
varcomp.df <- read_csv('output_brms/lanphere_VarComps.csv') %>% mutate(term_group="", term_order="") %>% as.data.frame()
aa.trait.PCA.2012 <- read_csv('final_data/aa.trait.PCA.2012.csv') %>% select(Trait=X1, PC1, PC2)
wind.trait.PCA.2012 <- read_csv('final_data/wind.trait.PCA.2012.csv') %>% select(Trait=X1, PC1, PC2)
wind.trait.PCA.2013 <- read_csv('final_data/wind.trait.PCA.2013.csv') %>% select(Trait=X1, PC1, PC2)
wind.soil.PCA <- read_csv('final_data/soil_PC_analysis.csv') %>% select(Soil_Property=X1, PC1=Comp.1, PC2=Comp.2)

## FILL IN TERM GROUPS AND ORDER THEM
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype__Intercept")] <- "Genotype (G)"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Species__Intercept")] <- "G x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Wind.Exposure__Intercept")] <- "Wind"
varcomp.df$term_group[which(varcomp.df$term=="sd_Wind.Exposure:Species__Intercept")] <- "Wind x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment__Intercept")] <- "Aphid"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Species__Intercept")] <- "Aphid x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Ant.mound.dist__Intercept")] <- "Ant"
varcomp.df$term_group[which(varcomp.df$term=="sd_Ant.mound.dist:Species__Intercept")] <- "Ant x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Ant.mound.dist__Intercept")] <- "Aphid x Ant"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Ant.mound.dist:Species__Intercept")] <- "Aphid x Ant x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Block__Intercept")] <- "Block"
varcomp.df$term_group[which(varcomp.df$term=="sd_Block:Species__Intercept")] <- "Block x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Plot_code__Intercept")] <- "Plot"
varcomp.df$term_group[which(varcomp.df$term=="sd_plant_ID__Intercept")] <- "Individual"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Wind.Exposure__Intercept")] <- "G x Wind"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Wind.Exposure:Species__Intercept")] <- "G x Wind x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment__Intercept")] <- "G x Aphid"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Species__Intercept")] <- "G x Aphid x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Ant.mound.dist__Intercept")] <- "G x Ant"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Ant.mound.dist:Species__Intercept")] <- "G x Ant x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Ant.mound.dist__Intercept")] <- "G x Aphid x Ant"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Ant.mound.dist:Species__Intercept")] <- "G x Aphid x Ant x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Species__Intercept")] <- "Species (Sp)"

varcomp.df$term_group <- factor(varcomp.df$term_group)
levels(varcomp.df$term_group)
varcomp.df$term_group_ord <- factor(varcomp.df$term_group, levels=c("Plot","Block x Sp","Block", "G x Aphid x Ant x Sp","G x Aphid x Ant", "G x Ant x Sp","G x Ant", "G x Aphid x Sp", "G x Aphid", "G x Wind x Sp","G x Wind","Aphid x Ant x Sp",  "Aphid x Ant", "Ant x Sp","Ant", "Aphid x Sp", "Aphid", "Wind x Sp","Wind", "G x Sp", "Genotype (G)","Species (Sp)"))
levels(varcomp.df$term_group_ord)

## FIGURE 1: COMMUNITY RICHNESS ----
plot_richness <- filter(varcomp.df, Response%in%c("Arthropod Richness", "Fungal Rarefied Richness", "Bacterial Rarefied Richness"), term != "sd_plant_ID__Intercept") %>% 
  droplevels() %>% 
  mutate(Experiment_Year=paste(Experiment, Year, " "))

# get R2 values ## IGNORING FOR NOW, TOO DIFFICULT TO PLOT ##
R2_richness <- plot_richness %>%
  group_by(Experiment, Year, Response) %>%
  summarise_at(vars(VarComp_mean), sum) %>%
  unite(Experiment_Year, Experiment, Year, sep=" ") 

richness_gg <- ggplot(plot_richness, aes(x=term_group_ord, y=SD_mean, shape=Response)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Effect size (SD)") +
  scale_shape_manual(values=c(21,22,23), name="Community response") +
  xlab("") +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y") 

## COMMUNITY COMPOSITION ----
plot_composition <- filter(varcomp.df, Response%in%c("Arthropod Composition"), term != "sd_plant_ID__Intercept") %>% #, "Fungal Composition", "Bacterial Composition"
  droplevels() %>% 
  mutate(Experiment_Year=paste(Experiment, Year, " "))

# get R2 values ## IGNORING FOR NOW, TOO DIFFICULT TO PLOT ##
R2_composition <- plot_composition %>%
  group_by(Experiment, Year, Response) %>%
  summarise_at(vars(VarComp_mean), sum) %>%
  unite(Experiment_Year, Experiment, Year, sep=" ") 

composition_gg <- ggplot(plot_composition, aes(x=term_group_ord, y=SD_mean, shape=Response)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Effect size (SD)") +
  scale_shape_manual(values=c(21,22,23), name="Community response") +
  xlab("") +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y") 
composition_gg

plot_grid(richness_gg, composition_gg) # playing around with adding the composition data too.

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

#### OLD BUT DEFINITELY USEFUL ----                           
## GENERAL PLOTTING FUNCTIONS
VarComp_plot <- function(plot_data){
  ggplot(plot_data, 
         aes(x=term, y=SD_mean, color=term, fill=term)) +
    geom_linerange(aes(ymin=lower_95, ymax=upper_95), size=0.5, show.legend=F) +
    geom_linerange(aes(ymin=lower_50, ymax=upper_50), size=2, show.legend=F) +
    geom_point(size=4, color="black", shape=21, show.legend=F) +
    #coord_flip() +
    xlab("") #
    #facet_wrap(~Experiment, ncol=1, scale="free_y")
}

## COLOR-BLIND FRIENDLY PALETTE
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
G.col <- cbPalette[4]
E.col <- cbPalette[3]
GxE.col <- cbPalette[2]
Plot.col <- cbPalette[1]
Block.col <- cbPalette[1]

## PLOTS
VarComp_plot(plot_data=filter(traits_2012, Trait=="Trait.PC1"))
VarComp_plot(plot_data=filter(traits_2012, Trait=="Trait.PC2"))

VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Arthropod Richness"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Arthropod Composition"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Fungal Rarefied Richness"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Bacterial Rarefied Richness"))


## STACKED BAR CHARTS OF VARIANCE COMPONENTS

varExp <- ggplot(filter(traits_2012, Trait=="Trait.PC1", Experiment=="Wind", Year=="2012"),
       aes(x=1, y=VarComp_mean*100, fill=term)) +
  geom_bar(stat = "identity", position=position_stack(reverse=T), color="black", show.legend = F) + coord_flip() +
  ylab("Trait PC1 explained variance (%)") + #scale_y_continuous(limits=c(0,100)) +
  #scale_fill_manual(values=c(G.col, E.col, GxE.col, Block.col, Plot.col)) +
  xlab("") + scale_x_continuous(labels=NULL, breaks=NULL) #+ theme(legend.position="top", legend.title = element_blank()))#theme(legend.justification=c(1,1), legend.position=c(1,0.65), legend.title = element_blank())

ggplot(filter(traits_2012, Trait=="Trait.PC1", Experiment=="Ant-Aphid", Year=="2012"),
       aes(x=1, y=VarComp_mean*100, fill=term)) +
  geom_bar(stat = "identity", position=position_stack(reverse=T), color="black") + coord_flip() +
  ylab("Trait PC1 explained variance (%)") + scale_y_continuous(limits=c(0,100)) +
  #scale_fill_manual(values=c(G.col, E.col, GxE.col, Block.col, Plot.col)) +
  xlab("") + scale_x_continuous(labels=NULL, breaks=NULL) + theme(legend.position="top", legend.title = element_blank())#theme(legend.justification=c(1,1), legend.position=c(1,0.65), legend.title = element_blank())

## TRAIT REGRESSION PLOTS
## NEED TO APPROPRIATELY SCALE BACK TRAIT AXES TO ORIGINAL SCALE
## NEED TO CREATE CUSTOM GROBS FOR PC LOADINGS.

## WIND ARTHROPOD RICHNESS 2012
w.rich.12 <- filter(arthr)

varComp.df <- filter(arthropods_2012, Experiment=="Wind", Year=="2012", Trait=="Arthropod Richness")
varComp.df$term <- factor(varComp.df$term,
                          levels=c("sd_Genotype__Intercept", 
                                   "sd_Wind.Exposure__Intercept",
                                   "sd_Genotype:Wind.Exposure__Intercept",
                                   "sd_Block__Intercept", "sd_Plot_code__Intercept"),
                          labels=c("G (Genotype)",
                                   "Wind",
                                   "G x Wind",
                                   "Block", "Plot"))
                             
                             
sum(varComp.df$VarComp_mean)
GxE <- VarComp_plot(plot_data=varComp.df) + ylab("Arthropod richness effect size (SD)") +
  annotate(geom="text", x=5, y=1.75, label="R^2 == 0.43", parse=T)

arrow.df <- w.trait.PC.12 %>% 
  mutate(x=0, y=seq(-1,-0.25,length.out=5), xend=PC1, yend=seq(-1,-0.25,length.out=5), hjust=c(-0.1,-0.1,-0.1,1.1,-0.1),
         label=c("Height","Shoot count", "Shoot length", "Water content", "Trichomes"))

traitReg <- plot(marginal_effects(trait.rich.wind.2012.brm, effects="sc.trait.PC1", probs=c(0.025,0.975)), 
     points=T, point_args=list(color="grey", shape=1))[[1]] + 
  ylab("Arthropod richness") + xlab("Trait PC1") +
  geom_segment(data=arrow.df, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length = unit(0.1,"cm")), inherit.aes=F) +
  geom_text(data=arrow.df, aes(x=xend, y=yend, label=label, hjust=hjust), size=3, inherit.aes=F) #+
 # theme(plot.margin = unit(c(1,1,1,1), "lines")) 
  #annotation_custom(grob=varGrob, xmin=-3, xmax=3, ymin=-5, ymax=-2)

plot_grid(GxE, traitReg, varExp, ncol=1, align = 'hv', rel_heights = c(1,1,0.4), scale=c(1,1,0.6))

#varcomp.df$term_group <- factor(varcomp.df$term_group, levels=c("G","E","G x E","E_block","E_plot","E_indiv")) # c("G","G x Sp","E","E x Sp","G x E","G x E x Sp","E_block","E_block x Sp","E_plot")
#varcomp.df$term_label <- factor(varcomp.df$term,
#                                levels=c("sd_Genotype__Intercept", 
#                                         "sd_Aphid.treatment__Intercept", "sd_Ant.mound.dist__Intercept", "sd_Aphid.treatment:Ant.mound.dist__Intercept", "sd_Wind.Exposure__Intercept",
#                                         "sd_Genotype:Aphid.treatment__Intercept", "sd_Genotype:Ant.mound.dist__Intercept", "sd_Genotype:Aphid.treatment:Ant.mound.dist__Intercept", "sd_Genotype:Wind.Exposure__Intercept",
#                                         "sd_Block__Intercept", "sd_Plot_code__Intercept"),
#                                labels=c("G (Genotype)",
#                                         "Aphid", "Ant", "Aphid x Ant", "Wind",
#                                         "G x Aphid", "G x Ant", "G x Aphid x Ant", "G x Wind",
#                                         "Block", "Plot"))

## USE SYMBOLS TO DENOTE ANT VS. APHID TREATMENTS
plot_data <- filter(varcomp.df, Response=="Arthropod Richness", Year=="2012", term != "sd_plant_ID__Intercept") %>% droplevels()
#plot_data$Symbols <- c(rep("General", 5), "Ant", "Aphid", "Ant x Aphid", rep("General",2), "Ant", "Aphid", "Ant x Aphid", "General")

plot.2012 <- ggplot(plot_data, 
                    aes(x=term_group, y=SD_mean, color=Experiment, fill=Experiment, group=interaction(term, Experiment))) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95), size=0.5, position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50), size=2, position=position_dodge(width=0.5)) +
  geom_point(size=5, position=position_dodge(width=0.5)) + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plot_data$term_group))) +
  xlab("") +
  scale_y_continuous(limits=c(0,2.1)) +
  scale_fill_manual(values=c("grey","black")) +
  scale_color_manual(values=c("grey","black")) #+
#scale_shape_manual(values = c(24,23,25,21)) 

plot_data <- filter(varcomp.df, Experiment=="Wind", Year=="2013", Response%in%c("Arthropod Richness", "Fungal Rarefied Richness", "Bacterial Rarefied Richness"), term != "sd_plant_ID__Intercept") %>% droplevels()
#plot_data$Symbols <- c(rep("General", 5), "Ant", "Aphid", "Ant x Aphid", rep("General",2), "Ant", "Aphid", "Ant x Aphid", "General")

#plot.2013 <- ggplot(plot_data, 
aes(x=term_group, y=SD_mean, color=Response, fill=Response, group=interaction(term, Response))) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95), size=0.5, position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50), size=2, position=position_dodge(width=0.5)) +
  geom_point(size=5, position=position_dodge(width=0.5)) + 
  coord_flip() +
  scale_y_continuous(limits=c(0,2.1)) +
  scale_x_discrete(limits = rev(levels(plot_data$term_group))) +
  xlab("") #+
#scale_fill_manual(values=c("grey","black")) +
#scale_color_manual(values=c("grey","black")) #+
#scale_shape_manual(values = c(24,23,25,21)) 

#plot_grid(plot.2012, plot.2013, align = 'hv', ncol=1)
