## Plots for manuscript

## load required libraries ----
library(dplyr)
library(tidyr)
#library(ggplot2)
library(cowplot)

## load required data ----
rootCN.means <- read.csv("rootCN.means.csv") %>% tbl_df() %>%
  mutate(experiment = "wind")
colnames(rootCN.means)[2:3] <- c("treatment","genotype")
levels(rootCN.means$treatment) <- c("Exposed","Unexposed")

arth.means <- read.csv("arth.means.csv") %>% tbl_df()
colnames(arth.means)[2:3] <- c("treatment","genotype")

trait.means <- read.csv("trait.means.csv") %>% tbl_df() %>% bind_rows(., rootCN.means) %>% bind_rows(., arth.means)

## useful settings ----
cbPal.10 <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFCC66")
ebar.w <- 0.05
l.size <- 1.5
alp <- 0.5
p.size <- 5

w.GxE <- function(resp, y.lab){
  ggplot(filter(trait.means, treatment != "NA", experiment == "wind", response == paste(resp)), aes(x = treatment, y = fit, color = genotype)) +
    geom_line(aes(group = genotype), size = l.size, alpha = alp) +
    scale_color_manual(values = cbPal.10) +
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", color = "black", width = ebar.w, size = l.size) +
    stat_summary(fun.y = mean, geom = "point", color = "black", size = p.size) +
    ylab(y.lab) +
    xlab("") +
    theme(legend.position = "none")
}

## Wind: GxE plots for individual and community phenotypes ----

# Individual phenotypes
w.height <- w.GxE(resp = "height", y.lab = "Plant height (cm)")
w.sl <- w.GxE(resp = "shoot.length", y.lab = "Avg. shoot length (cm)")
w.sc <- w.GxE(resp = "shoot.count", y.lab = "No. of shoots")
w.td <- w.GxE(resp = "leaf_trichome.density", y.lab = "No. trichomes")
w.WC <- w.GxE(resp = "leaf_WC", y.lab = "leaf water content")
w.SLA <- w.GxE(resp = "SLA", y.lab = "Specific leaf area")
w.CN <- w.GxE(resp = "leaf_CN", y.lab = "leaf C:N")
w.rCN <- w.GxE(resp = "root_CN", y.lab = "root C:N")

plot_grid(w.height, w.sl, w.sc, w.td, w.WC, w.SLA, w.CN, w.rCN,
          labels=c('A','B','C','D','E','F','G','H'))

# community phenotypes
w.arth <- w.GxE(resp = "arthropod_abund", y.lab = "No. of arthropod individuals")
w.rich <- w.GxE(resp = "arthropod_rich", y.lab = "No. of arthropod species")
w.rrich <- w.GxE(resp = "arthropod_rarerich", y.lab = "Rarefied richness")

plot_grid(w.arth, w.rich, w.rrich, 
          labels=c('A','B','C'))


## Ant-aphid: Demonstrating treatment effects on Aphis farinosa and F. obscuripes ----

# aphid response
(aphis.G <- ggplot(filter(arth.means, response == "Aphis_farinosa"), aes(x = Genotype, y = fit)) +
  geom_point(size = 10) +
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1) +
  scale_y_log10("Aphis farinosa abundance"))

# ant response
(Fobs.GxE <- ggplot(filter(arth.means, response == "F_obscuripes"), 
                    aes(x = Aphid.treatment, y = fit)) +
  geom_line(aes(group = Genotype, color = Genotype)))

(Fobs.ExE <- ggplot(filter(arth.means, response == "F_obscuripes", Ant.mound.dist > 0), 
                   aes(x = Ant.mound.dist, y = fit)) +
  geom_line(aes(group = Aphid.treatment, color = Aphid.treatment)))



