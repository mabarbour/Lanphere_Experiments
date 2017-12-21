
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)

## GET DATA ----
varcomp.df <- read_csv('output_brms/lanphere_VarComps.csv') %>% mutate(term_group="", term_order="") %>% as.data.frame()

## FILL IN TERM GROUPS AND ORDER THEM
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype__Intercept")] <- "Genotype (G)"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Species__Intercept")] <- "G x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Wind.Exposure__Intercept")] <- "Wind"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Wind.Exposure:Species__Intercept")] <- "E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment__Intercept")] <- "Aphid"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Species__Intercept")] <- "E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Ant.mound.dist__Intercept")] <- "Ant"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Ant.mound.dist:Species__Intercept")] <- "E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Ant.mound.dist__Intercept")] <- "Aphid x Ant"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Aphid.treatment:Ant.mound.dist:Species__Intercept")] <- "E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Block__Intercept")] <- "Block"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Block:Species__Intercept")] <- "E_block x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Plot_code__Intercept")] <- "Plot"
varcomp.df$term_group[which(varcomp.df$term=="sd_plant_ID__Intercept")] <- "Individual"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Wind.Exposure__Intercept")] <- "G x Wind"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Wind.Exposure:Species__Intercept")] <- "G x E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment__Intercept")] <- "G x Aphid"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Species__Intercept")] <- "G x E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Ant.mound.dist__Intercept")] <- "G x Ant"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Ant.mound.dist:Species__Intercept")] <- "G x E x Sp"
varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Ant.mound.dist__Intercept")] <- "G x Aphid x Ant"
#varcomp.df$term_group[which(varcomp.df$term=="sd_Genotype:Aphid.treatment:Ant.mound.dist:Species__Intercept")] <- "G x E x Sp"

varcomp.df$term_group <- factor(varcomp.df$term_group)
levels(varcomp.df$term_group)
varcomp.df$term_group_ord <- factor(varcomp.df$term_group, levels=c("Individual","Plot","Block","G x Aphid x Ant", "G x Ant", "G x Aphid", "G x Wind", "Aphid x Ant", "Ant", "Aphid", "Wind", "Genotype (G)"))
levels(varcomp.df$term_group_ord)

## FIGURE 1 ----
plot_data <- filter(varcomp.df, Response%in%c("Arthropod Richness", "Fungal Rarefied Richness", "Bacterial Rarefied Richness"), term != "sd_plant_ID__Intercept") %>% droplevels() %>% mutate(Experiment_Year=paste(Experiment, Year, " "))

ggplot(plot_data, aes(x=term_group_ord, y=SD_mean, shape=Response)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95), color="grey", size=0.5, position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=lower_50, ymax=upper_50), color="black", size=2, position=position_dodge(width=0.75)) +
  geom_point(size=4, color="black", fill="grey", position=position_dodge(width=0.75)) +
  coord_flip() +
  ylab("Effect size (SD)") +
  scale_shape_manual(values=c(21,22,23)) +
  xlab("") +
  facet_wrap(~Experiment_Year, ncol=1, scales="free_y") 

# get R2 values
plot_data %>%
  group_by(Experiment, Year, Response) %>%
  summarise_at(vars(VarComp_mean), sum) %>%
  unite(Experiment_Year, Experiment, Year, sep=" ") %>%
  mutate(R2_label = c(rep("R^2_arthropods",3), "R^2_bacteria", "R^2_fungi"))
  

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
