library(tidyverse)
library(cowplot)

traits_2012 <- read_csv("output_brms/lanphere_trait_VarComps_2012.csv") %>%
  mutate(term_group="")
arthropods_2012 <- read_csv("output_brms/lanphere_arthropods_VarComps_2012.csv")

traits_2012$term_group[which(traits_2012$term=="sd_Genotype__Intercept")] <- "G"
traits_2012$term_group[which(traits_2012$term=="sd_Wind.Exposure__Intercept")] <- "E"
traits_2012$term_group[which(traits_2012$term=="sd_Aphid.treatment__Intercept")] <- "E"
traits_2012$term_group[which(traits_2012$term=="sd_Ant.mound.dist__Intercept")] <- "E"
traits_2012$term_group[which(traits_2012$term=="sd_Aphid.treatment:Ant.mound.dist__Intercept")] <- "E"
traits_2012$term_group[which(traits_2012$term=="sd_Block__Intercept")] <- "E_site"
traits_2012$term_group[which(traits_2012$term=="sd_Plot_code__Intercept")] <- "E_site"
traits_2012$term_group[which(traits_2012$term=="sd_Genotype:Wind.Exposure__Intercept")] <- "GxE"
traits_2012$term_group[which(traits_2012$term=="sd_Genotype:Aphid.treatment__Intercept")] <- "GxE"
traits_2012$term_group[which(traits_2012$term=="sd_Genotype:Ant.mound.dist__Intercept")] <- "GxE"
traits_2012$term_group[which(traits_2012$term=="sd_Genotype:Aphid.treatment:Ant.mound.dist__Intercept")] <- "GxE"

traits_2012$term_group <- factor(traits_2012$term_group, levels=c("G","E","GxE","E_site"))
traits_2012$term <- factor(traits_2012$term,
                           levels=c("sd_Genotype__Intercept", 
                                    "sd_Aphid.treatment__Intercept", "sd_Ant.mound.dist__Intercept", "sd_Aphid.treatment:Ant.mound.dist__Intercept", "sd_Wind.Exposure__Intercept",
                                    "sd_Genotype:Aphid.treatment__Intercept", "sd_Genotype:Ant.mound.dist__Intercept", "sd_Genotype:Aphid.treatment:Ant.mound.dist__Intercept", "sd_Genotype:Wind.Exposure__Intercept",
                                    "sd_Block__Intercept", "sd_Plot_code__Intercept"),
                           labels=c("G (Genotype)",
                                    "Aphid", "Ant", "Aphid x Ant", "Wind",
                                    "G x Aphid", "G x Ant", "G x Aphid x Ant", "G x Wind",
                                    "Block", "Plot"))
## GENERAL PLOTTING FUNCTIONS ----
VarComp_plot <- function(plot_data){
  ggplot(plot_data, 
         aes(x=term, y=SD_mean, color=Experiment)) +
    geom_linerange(aes(ymin=lower_95, ymax=upper_95), size=0.5, color="grey") +
    geom_linerange(aes(ymin=lower_50, ymax=upper_50), size=2, show.legend=F) +
    geom_point(size=3, color="black") +
    coord_flip() +
    facet_wrap(~Experiment, ncol=1, scale="free_y")
}

## PLOTS
VarComp_plot(plot_data=filter(traits_2012, Trait=="Trait.PC1"))
VarComp_plot(plot_data=filter(traits_2012, Trait=="Trait.PC2"))

VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Arthropod Richness"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Arthropod Composition"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Fungal Rarefied Richness"))
VarComp_plot(plot_data=filter(arthropods_2012, Trait=="Bacterial Rarefied Richness"))
